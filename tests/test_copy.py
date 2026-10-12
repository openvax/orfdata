"""Unmodified CoPY source slices: calibration, censoring and cohort-only joins."""

import csv
import gzip
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import (
    DatasetManifest,
    DatasetRegistry,
    QuantitativeDataset,
    ScopedTarget,
    build_dataset,
)
from orfdata.cli import run
from orfdata.sources.copy import _quantity, _rna_rows, import_copy

FIXTURE = Path(__file__).parent / "data/copy"


@pytest.fixture(scope="module")
def source():
    manifest = DatasetManifest.from_dict(json.loads((FIXTURE / "manifest.json").read_text()))
    return manifest, {asset.asset_id: FIXTURE / asset.asset_id for asset in manifest.sources}


@pytest.fixture(scope="module")
def dataset(source, tmp_path_factory):
    manifest, paths = source
    return build_dataset(manifest, lambda i, o: import_copy(i, o, manifest),
                         cache_root=tmp_path_factory.mktemp("copy"), source_paths=paths,
                         download_options={"show_progress": False})


def test_native_protein_peptide_and_label_free_values(dataset):
    expected = json.loads((FIXTURE / "expected.json").read_text())["rows"]
    observations = {m.measurement_id: m for m in dataset.measurements()}
    for sheet, rows in expected.items():
        for line, entry in enumerate(rows, 2):
            row = entry["values"]
            columns = (("", 17),) if sheet == "ProteinQuant" else (
                tuple((f":Quant{i + 1}", i + 8) for i in range(4)) + ((":Median", 12),)
                if sheet == "PeptideQuant" else ((":Qexactive", 1), (":SAX", 2)))
            for suffix, column in columns:
                m = observations[f"{sheet}:{line}{suffix}"]
                assert tuple(row) == m.provenance["native_row"]
                native = row[column]
                if isinstance(native, (float, int)):
                    assert m.value == native and m.status == "measured"
                elif isinstance(native, str) and native.startswith("<"):
                    assert m.value is None and m.status == "below_detection"
                    assert m.detection_limit == float(native[1:])
                else:
                    assert m.value is None
                    assert m.status == ("missing" if native is None else "not_measured")
                assert m.unit == "copies/cell" and m.target.reference.taxon_id == 4932
    first = observations["ProteinQuant:2"]
    assert first.value == 3203.943997
    assert first.uncertainty == 154.1788355
    assert first.uncertainty_type == "median_absolute_deviation"
    duplicates = [m for m in observations.values() if m.target.identifier == "YNL239W"
                  and m.assay_id == "LabelFree:Qexactive"]
    assert len(duplicates) == 2 and duplicates[0].measurement_id != duplicates[1].measurement_id


def test_standard_identity_and_unresolved_assignments_remain_visible(dataset):
    mappings = list(dataset.mappings())
    mismatch = [m for m in mappings if m.evidence["standard"]["validation"] == "peptide_not_in_standard"]
    assert len(mismatch) == 4
    assert all(m.evidence["standard"]["sequence_sha256"] for m in mappings)
    assert all(m.evidence["standard"]["source_asset"] == "qconcat.fasta" for m in mappings)
    for mapping in dataset.reference_mappings():
        assert mapping.status in ("ambiguous", "unresolved") and not mapping.sequence_validated
        assert mapping.source.kind == "protein_group"
        assert all(m.measurement.target.kind == "protein_group" for m in
                   dataset.mapped_measurements(mapping.targets[0], mapping_statuses=[mapping.status]))
    reference = dataset.manifest.references[0]
    assert reference.sequence_sha256 is None  # Synthetic calibration standards are not a native proteome.
    assert reference.genome_build is None and reference.annotation_version is None


def test_rna_calibration_is_derived_and_distinct_from_replicate_fpkm(dataset):
    reference, rna_reference = dataset.manifest.references
    aggregate = list(dataset.measurements(target=ScopedTarget(reference, "gene", "YML035C", "sgd")))
    fpkm = next(m for m in aggregate if m.unit == "FPKM")
    copies = next(m for m in aggregate if m.unit == "copies/cell")
    assert fpkm.value == 40.5492 and fpkm.origin == "observed"
    assert copies.value == pytest.approx(2.1166512788981335)
    assert copies.origin == "derived" and copies.input_measurements == (fpkm.measurement_id,)
    assert copies.derivation["assumed_total_mrna_per_cell"] == 60000
    assert copies.derivation["absolute_RNA_directly_measured"] is False
    replicas = list(dataset.measurements(target=ScopedTarget(rna_reference, "transcript", "YML035C_mRNA", "CoPY:Cufflinks")))
    assert sorted(m.value for m in replicas) == [28.9895, 35.335, 48.3983, 49.1012]
    assert fpkm.value != pytest.approx(sum(m.value for m in replicas) / 4)
    lowdata = [m for m in dataset.measurements(reference=rna_reference) if m.status == "not_measured"]
    assert len(lowdata) == 4
    assert all(m.value is None and float(m.provenance["native_row"]["FPKM"]) > 0 for m in lowdata)
    assert all(m.provenance["confidence_interval"]["level"] == 0.95 for m in lowdata)
    assert len([m for m in dataset.measurements(reference=rna_reference) if m.value == 0]) == 4


def test_only_source_cohort_pairing_and_canonical_ratios(dataset):
    reference = dataset.manifest.references[0]
    assert list(dataset.matched_measurements(allowed_bases=["same_material"])) == []
    matches = list(dataset.matched_measurements(allowed_bases=["condition_aggregate"]))
    amd1 = next(m for m in matches if m.relation.left.identifier == "YML035C")
    assert amd1.pairing.evidence["individual_replicate_crosswalk"] == "unresolved"
    rna = next(m for m in amd1.left if m.unit == "copies/cell")
    protein = amd1.right[0]
    ratio = amd1.ratio(protein, rna, measurement_id="AMD1:protein-per-RNA", method="native calibrated ratio",
                       normalization={"RNA_total_mRNA_assumption": 60000, "level": "source cohort"})
    assert ratio.value == protein.value / rna.value and ratio.origin == "derived"
    group = next(m for m in matches if m.relation.relation == "gene_protein_group_context" and m.left)
    with pytest.raises(ValueError, match="group membership"):
        group.ratio(group.right[0], group.left[0], measurement_id="bad", method="bad", normalization={"bad": True})
    assert all(m.pairing.right_assay_id == "SRM:aggregate" for m in matches)
    assert all(s.biological_unit_id is None for s in dataset.samples() if s.sample_id != "cohort")
    assays = {a.assay_id: a for a in dataset.assays(reference=reference)}
    assert assays["LabelFree:SAX"].calibration["assumed_total_protein_molecules_per_cell"] == 60000000
    assert assays["LabelFree:Qexactive"].calibration["standard_fmol"] == 50
    assert assays["SRM:aggregate"].calibration["sample_cell_equivalents_per_injection"] == 200000


def test_source_build_hash_and_manifest_invariants(source, dataset, tmp_path):
    manifest, paths = source
    bad = tmp_path / "bad.xlsx"
    bad.write_bytes(b"different workbook")
    with pytest.raises(Exception, match="(?i)(hash|sha|size|integrity)"):
        build_dataset(manifest, lambda i, o: pytest.fail("Unverified input reached importer"),
                      cache_root=tmp_path, source_paths=dict(paths, **{"workbook.xlsx": bad}))
    coverage = dict(manifest.coverage, standard_sequences=92)
    with pytest.raises(ValueError, match="standard count"):
        import_copy(paths, tmp_path / "bad.sqlite", replace(manifest, coverage=coverage))
    wrong = replace(manifest.references[0], source_sha256="f" * 64)
    processing = dict(manifest.processing, reference_roles=dict(workbook=wrong.key, rna_replicates=manifest.references[1].key))
    with pytest.raises(ValueError, match="pinned sources"):
        import_copy(paths, tmp_path / "bad.sqlite", replace(manifest, references=(wrong, manifest.references[1]), processing=processing))
    assert not (tmp_path / "bad.sqlite").exists()
    assert all(m.target.reference.taxon_id == 4932 for m in dataset.measurements())


@pytest.mark.parametrize("value", [-1, float("inf"), float("nan"), True, "<no", "unknown"])
def test_malformed_quantities_rejected(value):
    with pytest.raises((TypeError, ValueError)):
        _quantity(value)


def test_rna_empty_gene_id_valid_but_duplicate_tracking_id_rejected(tmp_path):
    path = next(FIXTURE.glob("*.gz"))
    with gzip.open(path, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0]["gene_id"] = ""
    output = tmp_path / "rna.gz"
    def write(values):
        with gzip.open(output, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(values)
    write(rows)
    assert len(list(_rna_rows(output))) == 3
    write(rows + [rows[0]])
    with pytest.raises(ValueError, match="Duplicate RNA"):
        list(_rna_rows(output))


def test_cli_offline_source_rebuild(source, dataset, tmp_path, capsys):
    from test_measurements import definition

    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps(definition(dataset)))
    options = ["datasets", "--registry", str(registry), "--cache-root", str(tmp_path / "cache"),
               "build", source[0].dataset_id, "--quiet"]
    for key, path in source[1].items():
        options += ["--source-path", f"{key}={path}"]
    run(options)
    built = QuantitativeDataset(json.loads(capsys.readouterr().out)["path"])
    assert built.counts == dataset.counts and built.manifest == dataset.manifest


def test_registry_yeast_scope_and_parallel_dataset_api():
    registry = DatasetRegistry()
    yeast = registry.info("copy-yeast-abundance")
    assert all(r.taxon_id == 4932 for r in yeast.references)
    assert yeast in registry.list(species="yeast", assay="targeted_ms")
    assert yeast not in registry.list(species="human")
    assert yeast.coverage["quantified_protein_or_groups"] == 1167
    assert yeast.coverage["upper_bound_protein_or_groups"] == 668
    assert yeast.coverage["label_free_duplicate_labels"] == {"YNL239W": 2}
    assert registry.download_metadata(yeast.dataset_id)["distribution"]["redistributable"] is True
