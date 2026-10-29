"""Native paired study facts, missingness, controls and deliberately limited identity."""

import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import AssayPair, DatasetManifest, DatasetRegistry
from orfdata.sources.mammary_stroma import (
    clinical_cases,
    geo_samples,
    import_stroma,
    specimen_crosswalk,
)

DIRECTORY = Path(__file__).parent / "data/mammary_stroma"


@pytest.fixture
def study(tmp_path):
    manifest = DatasetManifest.from_dict(json.loads((DIRECTORY / "manifest.json").read_text()))
    paths = {}
    for asset in manifest.sources:
        path = tmp_path / asset.asset_id
        path.write_bytes((DIRECTORY / asset.asset_id).read_bytes())
        paths[asset.asset_id] = path
    return paths, manifest


def repin(manifest, name, path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    references = tuple(
        replace(r, source_sha256=digest)
        if r.source_sha256 == next(a.sha256 for a in manifest.sources if a.asset_id == name)
        else r
        for r in manifest.references
    )
    return replace(
        manifest,
        references=references,
        sources=tuple(
            replace(a, sha256=digest, size=path.stat().st_size) if a.asset_id == name else a
            for a in manifest.sources
        ),
    )


def test_native_assays_and_metadata_derived_specimens(study, tmp_path):
    paths, manifest = study
    dataset = import_stroma(paths, tmp_path / "dataset.sqlite", manifest)
    assert dataset.counts["sample"] == dataset.counts["assay"] == 58
    pairs = list(dataset.pairings())
    assert len(pairs) == 30
    paired = [p for p in pairs if p.right_assay_id is not None]
    missing = [p for p in pairs if p.right_assay_id is None]
    assert len(paired) == 28 and len(missing) == 2
    assert {p.evidence["rna_case"] for p in missing} == {"127"}
    assert all(
        p.evidence["protein_case"] is None and p.temporal_relation == "unknown" for p in missing
    )
    assert all(p.basis == "same_specimen" for p in pairs)
    assert all(p.reference != p.right_reference for p in pairs)
    assert all(AssayPair.from_dict(p.to_dict()) == p for p in pairs)
    assert not list(dataset.relations())
    assert not list(dataset.reference_mappings())
    samples = {(s.reference.key, s.sample_id): s for s in dataset.samples()}
    for pair in paired:
        left = samples[pair.reference.key, pair.left_assay_id]
        right = samples[pair.right_reference.key, pair.right_assay_id]
        assert left.biological_unit_id == right.biological_unit_id
        assert left.condition == right.condition
        assert left.sample_id != right.sample_id
        assert not left.pool_members and not right.pool_members
        assert pair.evidence["crosswalk_status"].startswith("metadata_derived")
    assert all(
        r.annotation_version is None and r.sequence_sha256 is None for r in manifest.references
    )
    assert {s.condition for s in samples.values()} == {
        "normal stroma in a cancer-bearing mammary specimen",
        "cancer-associated mammary stroma",
    }


def test_quantities_preserve_source_values_missing_sentinels_and_controls(study, tmp_path):
    paths, manifest = study
    dataset = import_stroma(paths, tmp_path / "native.sqlite", manifest)
    observations = list(dataset.measurements())
    expected = 3 * 30 + sum(
        sum(t["measurement_statuses"].values()) for t in manifest.coverage["tables"].values()
    )
    assert len(observations) == expected
    assert {m.unit for m in observations} == {
        "estimated_count",
        "intensity",
        "LFQ intensity",
        "iBAQ",
    }
    zeros = [m for m in observations if m.value == 0]
    assert zeros and all(m.quantity == "rna_abundance" and m.status == "measured" for m in zeros)
    nans = [m for m in observations if m.provenance["native_value"] == "NaN"]
    assert nans and all(
        m.value is None and m.status == "missing" and m.provenance["missing_reason"] == "native NaN"
        for m in nans
    )
    assert all(m.provenance["source_qc"]["Reverse"] == "+" for m in nans)
    for measurement in observations:
        native = measurement.provenance["native_value"]
        if measurement.quantity != "rna_abundance" and native in ("0", "NaN"):
            assert measurement.value is None and measurement.status == "missing"
        else:
            assert measurement.value == float(native) and measurement.status == "measured"
    mappings = list(dataset.mappings())
    proteins = [m for m in mappings if m.source_namespace.endswith("proteinGroups.txt")]
    peptides = [m for m in mappings if m.source_namespace.endswith("peptides.txt")]
    assert len(proteins) == 8 and len(peptides) == 4
    assert all(not m.targets and m.status == "unresolved" for m in proteins)
    assert any(m.status == "ambiguous" and len(m.targets) > 1 for m in peptides)
    assert all(t.kind == "protein_group" for m in peptides for t in m.targets)
    assert {m.evidence["source_partition"] for m in proteins} == {"target", "contaminant", "decoy"}
    assert all(m.evidence["unique_protein_detection"] is False for m in mappings)
    # Every original MS field survives once, independently of the observation projection.
    for name in ("proteinGroups.txt", "peptides.txt"):
        with paths[name].open() as stream:
            rows = {row["id"]: row for row in csv.DictReader(stream, delimiter="\t")}
        for mapping in mappings:
            if mapping.provenance["source_asset"] == name:
                assert (
                    mapping.to_dict()["provenance"]["native_record"]
                    == rows[mapping.source_identifier]
                )


def test_metadata_crosswalk_rejects_ambiguity_and_compartment_conflict(study):
    paths, manifest = study
    clinical, rna = clinical_cases(paths["paper.xml"]), geo_samples(paths["geo.soft.gz"])
    assert specimen_crosswalk(clinical, rna) == dict(manifest.processing["specimen_crosswalk"])
    rna["139_Norm"]["traits"] = dict(rna["135_Norm"]["traits"])
    with pytest.raises(ValueError, match="compartments disagree"):
        specimen_crosswalk(clinical, rna)
    rna["139_Tumour"]["traits"] = dict(rna["135_Tumour"]["traits"])
    with pytest.raises(ValueError, match="ambiguous"):
        specimen_crosswalk(clinical, rna)


def test_source_checksum_and_reference_are_not_silently_repaired(study, tmp_path):
    paths, manifest = study
    changed = replace(
        manifest,
        references=(
            manifest.references[0],
            replace(manifest.references[1], annotation_version="invented"),
        ),
    )
    with pytest.raises(ValueError, match="reference identities"):
        import_stroma(paths, tmp_path / "reference.sqlite", changed)
    paths["proteinGroups.txt"].write_bytes(paths["proteinGroups.txt"].read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        import_stroma(paths, tmp_path / "corrupt.sqlite", manifest)
    assert not (tmp_path / "corrupt.sqlite").exists()


@pytest.mark.parametrize("mutation", ["duplicate", "negative", "flag"])
def test_invalid_native_rows_fail_without_a_partial_dataset(study, tmp_path, mutation):
    paths, manifest = study
    path = paths["proteinGroups.txt"]
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    if mutation == "duplicate":
        rows[1]["id"] = rows[0]["id"]
    elif mutation == "negative":
        rows[0][next(k for k in fields if k.startswith("Intensity "))] = "-1"
    else:
        rows[0]["Reverse"] = "unknown"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = repin(manifest, "proteinGroups.txt", path)
    output = tmp_path / "invalid.sqlite"
    with pytest.raises(ValueError):
        import_stroma(paths, output, manifest)
    assert not output.exists()


def test_registry_build_and_offline_queries_use_common_surface(study, tmp_path):
    paths, manifest = study
    registry = DatasetRegistry(
        {
            manifest.dataset_id: {
                "default_version": manifest.version,
                "versions": {
                    manifest.version: {
                        "manifest": manifest.to_dict(),
                        "assets": {},
                        "distribution": {
                            "redistributable": True,
                            "access_terms": "Factual fixture",
                        },
                    }
                },
            }
        },
        cache_root=tmp_path / "cache",
    )
    registry.build(manifest.dataset_id, source_paths=paths)
    dataset = registry.open_build(manifest.dataset_id)
    assert len(list(dataset.pairings())) == 30
    assert dataset.manifest == manifest
    assert registry.info(manifest.dataset_id).references[0].taxon_id == 9615
