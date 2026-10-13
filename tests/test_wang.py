"""Native healthy-tissue observations keep censoring, calibration and group scope."""

import io
import json
import zipfile
from dataclasses import replace
from pathlib import Path

import openpyxl
import pytest

from orfdata import DatasetManifest, DatasetRegistry, ScopedTarget, build_dataset
from orfdata.sources.wang import SERIES, _number, import_tissue

ROOT = Path(__file__).parent / "data/wang"


@pytest.fixture
def source():
    manifest = DatasetManifest.from_dict(json.loads((ROOT / "manifest.json").read_text()))
    return manifest, {s.asset_id: ROOT / s.asset_id for s in manifest.sources}


def build(source, tmp_path):
    manifest, inputs = source
    return build_dataset(manifest, lambda i, o: import_tissue(i, o, manifest),
                         source_paths=inputs, cache_root=tmp_path)


def test_all_real_values_keep_native_scope_and_zero_meaning(source, tmp_path):
    dataset = build(source, tmp_path)
    expected = json.loads((ROOT / "expected.json").read_text())["rows"]
    samples = {s.tissue: s.sample_id for s in dataset.samples()}
    assert samples == {"Adrenal gland": "V122", "Colon": "V269", "Duodenum": "V145"}
    measurements = {m.measurement_id: m for m in dataset.measurements()}
    assert len(measurements) == 237
    positive_subnormal = []
    for name, rows in expected.items():
        nids = SERIES[name][3]
        for line, original in enumerate(rows, 2):
            row = original["native_row"]
            for tissue, column, native in zip(source[0].processing["tissues"], original["source_columns"], row[nids:]):
                sid = samples[tissue]
                value = measurements[f"{name}:{line}:{sid}"]
                assert value.target.reference == source[0].references[0]
                assert value.sample_id == sid
                assert value.provenance["native_value"] == native
                assert value.provenance["column"] == column and value.provenance["row"] == line
                assert value.provenance["source_asset"] == SERIES[name][0]
                expected_status = "measured"
                if float(native) == 0 and name != "rna_transcripts":
                    expected_status = "below_detection" if name.startswith("rna_") else "not_measured"
                assert value.status == expected_status
                assert value.detection_limit is None  # Native normalized-unit limits are not supplied.
                assert value.value == (float(native) if expected_status == "measured" else None)
                assert value.unit == ("copies/cell" if name.endswith("copies") else "FPKM" if name.startswith("rna_") else "iBAQ")
                if isinstance(native, str) and 0 < float(native) < 1e-300:
                    positive_subnormal.append(value)
                if name.endswith("copies"):
                    assert value.origin == "derived" and value.derivation["method"]
                    if name == "rna_copies":
                        original_value, = [measurements[x] for x in value.input_measurements]
                        assert original_value.target == value.target and original_value.unit == "FPKM"
                        assert original_value.sample_id == value.sample_id
                        assert value.derivation["assumed_mrna_fraction_of_total_rna_mass"] == 0.02
                    else:
                        assert value.target.kind == "protein_group"
                        assert value.input_measurements == ()  # iBAQ is not the raw ruler input.
                        assert measurements[value.provenance["associated_iBAQ_measurement"]].target == value.target
                else:
                    assert value.origin == "observed"
    assert positive_subnormal
    assert all(m.value > 0 and m.status == "measured" for m in positive_subnormal)
    assert any(m.status == "measured" and m.value == 0 for m in measurements.values())
    assert any(m.unit == "copies/cell" and m.status == "below_detection"
               and measurements[m.input_measurements[0]].value > 0
               for m in measurements.values() if m.input_measurements)
    assert all(s.attributes["donor_id"] is None and s.biological_replicate is None for s in dataset.samples())
    assert not any(s.tissue == "Bone marrow" for s in dataset.samples())
    assert len(list(dataset.assays())) == 12  # Four processing views, three specimens.


def test_pairing_and_native_gene_ratio_never_become_group_or_orf_ratio(source, tmp_path):
    dataset = build(source, tmp_path)
    assert list(dataset.matched_measurements(allowed_bases=["same_material"])) == []
    relations = list(dataset.relations())
    gene = next(r for r in relations if r.relation == "same_target" and r.left.identifier == "ENSG00000004059")
    pair = "rna_genes:protein_genes:V122"
    match, = dataset.matched_measurements(allowed_bases=["adjacent_donor_tissue"],
                                         pair_ids=[pair], relation_ids=[gene.relation_id])
    assert len(match.left) == len(match.right) == 1
    ratio = match.ratio(match.right[0], match.left[0], measurement_id="relative-protein-per-RNA",
                        method="explicit native normalized iBAQ per FPKM",
                        normalization={"scope": "source gene-level tables; not molecular copy ratio"})
    assert ratio.quantity == "protein_rna_ratio" and ratio.origin == "derived"
    assert ratio.target.kind == "gene"
    group = next(r for r in relations if r.relation == "gene_protein_group_context"
                 and r.left.identifier == "ENSG00000004059")
    match, = dataset.matched_measurements(allowed_bases=["adjacent_donor_tissue"],
                                         pair_ids=["rna_genes:protein_groups:V122"], relation_ids=[group.relation_id])
    with pytest.raises(ValueError, match="group membership"):
        match.ratio(match.right[0], match.left[0], measurement_id="invalid", method="invalid",
                    normalization={"scope": "test"})
    assert not any(m.target.kind == "orf" for m in dataset.measurements())


def test_representatives_keep_group_ambiguity_and_native_name_differences(source, tmp_path):
    dataset = build(source, tmp_path)
    mapping = next(m for m in dataset.reference_mappings() if m.validation["native_gene_name"] == "INO80B;INO80B-WBP1")
    assert mapping.status == "unresolved" and len(mapping.targets) == 2
    assert not mapping.sequence_validated
    target = mapping.targets[0]
    assert list(dataset.measurements(target=target)) == []
    assert list(dataset.mapped_measurements(target)) == []
    views = list(dataset.mapped_measurements(target, mapping_statuses=["unresolved"], unit="copies/cell"))
    assert len(views) == 3
    assert all(v.measurement.target.kind == "protein_group" for v in views)
    assert {v.measurement.provenance["native_gene_label"] for v in views} == {"INO80B"}
    assert mapping.validation["native_gene_ids"] == views[0].measurement.provenance["native_gene_ids"]


def test_source_build_and_native_search_reference_are_offline(source, tmp_path, monkeypatch):
    import requests
    from test_measurements import definition

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Offline source fixture contacted network"))
    first = build(source, tmp_path / "build")
    assert build(source, tmp_path / "build").path == first.path
    registry = DatasetRegistry(definition(first), cache_root=tmp_path / "registry")
    assert registry.inspect(first.manifest.dataset_id).status == "missing"
    registry.download_reference(first.manifest.dataset_id, source_path=source[1]["search.fasta.gz"])
    search = registry.open_reference(first.manifest.dataset_id)
    assert len(search) == source[0].coverage["search_sequences"]
    assert search.get("ENSP00000484462.1") is not None
    assert search.get("ENSP00000484462") is None
    assert registry.inspect(first.manifest.dataset_id).status == "missing"
    reference = first.manifest.references[0]
    with pytest.raises(ValueError, match="Reference context"):
        list(first.measurements(target=ScopedTarget(replace(reference, taxon_id=10090, species="mouse"),
                                                    "gene", "ENSG00000004059", "ensembl")))


def test_real_pyensembl_gene_scopes_both_native_gene_measurements(source, tmp_path):
    import gzip

    from pyensembl import Genome

    dataset = build(source, tmp_path / "data")
    reference = dataset.manifest.references[0]
    gtf = tmp_path / "annotation.gtf"
    gtf.write_bytes(gzip.decompress(source[1]["annotation.gtf.gz"].read_bytes()))
    genome = Genome(reference_name="GRCh38", annotation_name="ensembl", annotation_version=83,
                    gtf_path_or_url=str(gtf), cache_directory_path=str(tmp_path / "genome"))
    genome.species = "human"
    genome.index()
    gene = genome.gene_by_id("ENSG00000004059")
    target = ScopedTarget.from_pyensembl(gene, kind="gene", namespace="ensembl", context=reference)
    assert len(list(dataset.measurements(target=target, unit="FPKM"))) == 3
    assert len(list(dataset.measurements(target=target, unit="iBAQ"))) == 3
    with pytest.raises(ValueError, match="annotation_version"):
        ScopedTarget.from_pyensembl(gene, kind="gene", namespace="ensembl",
                                   context=replace(reference, annotation_version="84"))


@pytest.mark.parametrize("value", (True, None, -1, float("nan"), float("inf"), "NA", "1e-999"))
def test_invalid_numeric_cells_are_never_silently_coerced(value):
    with pytest.raises(ValueError):
        _number(value)


@pytest.mark.parametrize("change", ("negative", "rna_relationship", "copy_identifier", "duplicate"))
def test_bad_source_cells_fail_atomically(source, tmp_path, change):
    manifest, inputs = source
    if change == "copy_identifier":
        asset, member, sheet, cell, value = "copies.zip", "Table_EV5.xlsx", "A. Protein copies", "B2", "wrong-gene"
    else:
        asset, member, sheet = "rna.zip", "Table_EV2.xlsx", "A. Transcripts"
        cell, value = ("E2", -1) if change == "negative" else ("B2", "wrong-protein")
    with zipfile.ZipFile(inputs[asset]) as archive:
        book = openpyxl.load_workbook(io.BytesIO(archive.read(member)))
    if change == "duplicate":
        book[sheet].append([c.value for c in book[sheet][2]])
    else:
        book[sheet][cell] = value
    raw = io.BytesIO()
    book.save(raw)
    book.close()
    altered = tmp_path / asset
    with zipfile.ZipFile(altered, "w") as archive:
        archive.writestr(member, raw.getvalue())
    output = tmp_path / "invalid.sqlite"
    with pytest.raises((ValueError, KeyError)):
        import_tissue(dict(inputs, **{asset: altered}), output, manifest)
    assert not output.exists()


def test_historical_patch_identity_is_not_dropped(source, tmp_path):
    manifest, inputs = source
    altered = replace(manifest, processing=dict(manifest.processing, source_assembly_header="GRCh38.p6"))
    with pytest.raises(ValueError, match="GTF assembly"):
        import_tissue(inputs, tmp_path / "wrong-patch.sqlite", altered)
    assert not (tmp_path / "wrong-patch.sqlite").exists()
