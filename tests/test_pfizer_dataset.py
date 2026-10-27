"""Common queries retain every native cell when stored as compact vectors."""

import hashlib
import json
import sqlite3
from dataclasses import replace

import pytest
from test_pfizer_evidence import build
from test_pfizer_evidence import source as source
from test_pfizer_lineage import quantification

from orfdata import (
    AssayPair,
    QuantitativeDataset,
    ReferenceMapping,
    ScopedTarget,
    SourceAsset,
    TargetRelation,
)
from orfdata.sources.pfizer import protein_manifest
from orfdata.sources.pfizer_dataset import copy_evidence, install_projection
from orfdata.sources.pfizer_evidence import EvidenceStore, decode_scalars
from orfdata.sources.pfizer_lineage import lineage_records


@pytest.fixture
def projected(source, tmp_path):
    evidence_path, _ = build(source[0], tmp_path)
    lineage, tables, converter = quantification()
    values = [m for table, columns in converter.COLUMNS.items() for row in tables[table]
              for column in columns for m in converter.peptide_abundances(table, row, column)]
    values.sort(key=lambda m: m.measurement_id)
    reference = converter.reference
    destination = replace(reference, sequence_sha256="f" * 64)
    first = next(m for m in values if m.status == "measured")
    mapped_target = ScopedTarget(destination, "peptide", "synthetic-alias", "explicit-test-alias")
    mapping = ReferenceMapping("fixture-map", converter.dataset_id, converter.dataset_version,
                               first.target, destination, (mapped_target,), "unique", "explicit-test-alias",
                               {"scope": "Query regression only; not a biological assignment"})
    samples, assays = lineage_records(lineage, reference)
    other_assay = next(a.assay_id for a in assays if a.sample_id == first.sample_id and a.assay_id != first.assay_id)
    pair = AssayPair("fixture-pair", reference, converter.dataset_id, converter.dataset_version,
                     first.assay_id, other_assay, "same_material",
                     {"scope": "Same-pool assay query regression only"}, "concurrent")
    relation = TargetRelation("fixture-relation", converter.dataset_id, converter.dataset_version,
                              first.target, first.target, "same_target", "one_to_one",
                              {"scope": "Identical target query regression"})
    asset = SourceAsset("native.json", evidence_path.as_uri(), hashlib.sha256(evidence_path.read_bytes()).hexdigest(),
                        evidence_path.stat().st_size, "Factual test slice; see vendor/NOTICE.md")
    manifest = replace(protein_manifest("dog"), dataset_id=converter.dataset_id, version=converter.dataset_version,
                       title="Explicit compact native quantification regression slice", sources=(asset,),
                       references=(reference, destination),
                       coverage={"scope": "Native regression slice, not a complete published dataset"})
    common = dict(samples=samples, assays=assays, reference_mappings=(mapping,), pairings=(pair,), relations=(relation,))
    flat = QuantitativeDataset.create(tmp_path / "flat.sqlite", manifest, measurements=values, **common)
    vector = QuantitativeDataset.create(tmp_path / "vector.sqlite", manifest, measurements=(), **common)
    with sqlite3.connect(vector.path) as connection:
        connection.execute("PRAGMA foreign_keys=ON")
        copy_evidence(connection, evidence_path)
        result = install_projection(connection, manifest, lineage, source_asset="native.json", reference=reference)
    return QuantitativeDataset(vector.path), flat, values, first, mapped_target, result


def test_all_fields_order_and_logical_counts_match_flat_dataset(projected):
    vector, flat, values, _, _, result = projected
    assert list(vector.measurements()) == list(flat.measurements()) == values
    assert vector.counts == flat.counts
    assert result == dict(vectors=12, measurements=429, layouts=4)
    assert list(vector.samples()) == list(flat.samples())
    assert list(vector.assays()) == list(flat.assays())
    assert list(vector.reference_mappings()) == list(flat.reference_mappings())
    assert list(EvidenceStore(vector.path).records("TargetPsms")) == list(EvidenceStore(flat.manifest.sources[0].url[7:]).records("TargetPsms"))
    with sqlite3.connect(vector.path) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        with pytest.raises(sqlite3.OperationalError, match="no such table"):
            connection.execute("SELECT count(*) FROM measurement")


def test_original_object_encoding_remains_queryable_through_common_api(projected):
    vector, flat, values, _, _, _ = projected
    with sqlite3.connect(vector.path) as db:
        layouts = {index: json.loads(columns) for index, columns in db.execute("SELECT i,columns_json FROM evidence_table")}
        for identifier, table, payload in db.execute("SELECT i,table_i,attributes FROM evidence_record").fetchall():
            db.execute("UPDATE evidence_record SET attributes=? WHERE i=?",
                       (json.dumps(decode_scalars(payload, layouts[table]), sort_keys=True), identifier))
        db.execute("UPDATE evidence_metadata SET value='\"orfdata.pfizer-evidence.v1\"' WHERE name='encoding'")
    legacy = QuantitativeDataset(vector.path)
    assert list(legacy.measurements()) == list(flat.measurements()) == values
    assert list(legacy.measurements(status="measured")) == list(flat.measurements(status="measured"))


def test_native_selectors_and_mapped_paired_views_match_flat_storage(projected):
    vector, flat, _, first, mapped_target, _ = projected
    reference = first.target.reference
    queries = [dict(target=first.target), dict(reference=reference), dict(sample_id=first.sample_id, reference=reference),
               dict(assay_id=first.assay_id, reference=reference), dict(status="measured"), dict(status="missing"),
               dict(status="below_detection"), dict(unit=first.unit), dict(unit="different unit"),
               dict(scale="linear"), dict(scale="log"), dict(log_base=2), dict(quantity="protein_abundance"),
               dict(quantity="peptide_abundance"), dict(target=mapped_target),
               dict(target=first.target, sample_id=first.sample_id, status="measured", unit=first.unit)]
    for query in queries:
        assert list(vector.measurements(**query)) == list(flat.measurements(**query)), query
    assert list(vector.mapped_measurements(mapped_target, status="measured")) == list(
        flat.mapped_measurements(mapped_target, status="measured"))
    for query in (dict(allowed_bases=["same_material"]),
                  dict(allowed_bases=["same_material"], left_target=first.target, right_filters={"unit": first.unit})):
        assert list(vector.matched_measurements(**query)) == list(flat.matched_measurements(**query))
    for query in (dict(status="invalid"), dict(quantity="invalid"), dict(unit=""),
                  dict(scale="linear", log_base=2), dict(sample_id=first.sample_id)):
        with pytest.raises((ValueError, TypeError)):
            list(vector.measurements(**query))
    assert list(vector._measurements(target_keys=())) == []


def test_columnar_export_and_reader_completion_gate(projected, tmp_path):
    pq = pytest.importorskip("pyarrow.parquet")
    vector, flat, _, first, _, _ = projected
    for name, dataset in (("vector", vector), ("flat", flat)):
        dataset.export_parquet(tmp_path / (name + ".parquet"), target=first.target)
    assert pq.read_table(tmp_path / "vector.parquet").to_pylist() == pq.read_table(tmp_path / "flat.parquet").to_pylist()
    with sqlite3.connect(vector.path) as connection:
        value = json.loads(connection.execute("SELECT value FROM metadata WHERE name='measurement_storage'").fetchone()[0])
        value["state"] = "building"
        connection.execute("UPDATE metadata SET value=? WHERE name='measurement_storage'", (json.dumps(value),))
    with pytest.raises(ValueError, match="incomplete"):
        QuantitativeDataset(vector.path)
