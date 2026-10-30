"""Native count slices retain gene, sample, strain and historical-reference scope."""

import csv
import gzip
import io
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry, ScopedTarget, build_dataset
from orfdata.cli import run
from orfdata.sources.pfizer_rna import import_rna

ROOT = Path(__file__).parent / "data/pfizer_rna"


@pytest.fixture(params=("dog", "mouse"))
def source(request):
    species = request.param
    manifest = DatasetManifest.from_dict(json.loads((ROOT / f"{species}-manifest.json").read_text()))
    expected = json.loads((ROOT / f"{species}-expected.json").read_text())
    inputs = {a.asset_id: ROOT / (f"{species}-annotation.gtf.gz" if a.asset_id == "annotation.gtf.gz" else a.asset_id)
              for a in manifest.sources}
    return manifest, expected, inputs


def build(tmp_path, source):
    manifest, _, inputs = source
    return build_dataset(manifest, lambda i, o: import_rna(i, o, manifest), source_paths=inputs,
                         cache_root=tmp_path, download_options={"show_progress": False})


def test_exact_native_counts_zeroes_and_source_identity(source, tmp_path):
    manifest, expected, _ = source
    dataset = build(tmp_path, source)
    observations = list(dataset.measurements(unit="counts", scale="linear"))
    assert len(observations) == len(expected["rows"]) == manifest.coverage["rows"]
    by_line = {o.provenance["line"]: o for o in observations}
    for native in expected["rows"]:
        observation = by_line[native["fixture_line"]]
        row = native["row"]
        assert observation.value == int(row["rcnt"])
        assert observation.status == "measured" and observation.origin == "observed"
        assert observation.target.kind == "gene" and observation.target.identifier == row["egid"]
        assert observation.target.namespace == "ensembl"
        assert observation.target.reference == manifest.references[0]
        assert observation.sample_id == observation.assay_id == row["sid"]
        assert observation.provenance["native_value"] == row["rcnt"]
        assert observation.provenance["line"] == native["fixture_line"]
    assert any(o.value == 0 for o in observations)
    assert list(dataset.measurements(unit="TPM")) == []
    assert list(dataset.pairings()) == []
    assert list(dataset.measurements(target=ScopedTarget(manifest.references[0], "gene", expected["unmeasured_gene"], "ensembl"))) == []
    assert all(s.attributes["donor_id"] is None and s.biological_replicate is None for s in dataset.samples())
    assert all(a.calibration["absolute_molecule_counts"] is False for a in dataset.assays())
    expected_table = {str(x["row"][0]): x["row"] for x in expected["sample_table"]}
    for sample in dataset.samples():
        assert sample.biological_unit_id == expected_table[sample.sample_id][-1]
        assert sample.attributes["source_sample_table"]["mapped.percent"] == expected_table[sample.sample_id][9]
    if manifest.processing["native_species"] == "mouse":
        assert {s.attributes["biological_strain"] for s in dataset.samples()} == {"c57", "cd1"}
    else:
        sample = next(s for s in dataset.samples() if s.sample_id == "1")
        assert sample.tissue == "Adipose"
        assert "tissue: AdiposeRepro" in sample.attributes["source_geo"]["characteristics_ch1"]


def test_offline_source_reuse_and_checksum_failure(source, tmp_path, monkeypatch):
    import requests

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Local source build contacted network"))
    first = build(tmp_path / "cache", source)
    assert build(tmp_path / "cache", source).path == first.path
    manifest, _, inputs = source
    corrupt = tmp_path / "bad.csv.gz"
    corrupt.write_bytes(b"bad")
    with pytest.raises(Exception, match="(?i)(hash|sha|size|integrity)"):
        build_dataset(manifest, lambda i, o: pytest.fail("unverified inputs reached importer"),
                      source_paths=dict(inputs, **{"counts.csv.gz": corrupt}), cache_root=tmp_path / "bad")


def test_species_reference_isolation_and_pyensembl_target(source, tmp_path):
    from pyensembl import Genome

    manifest, _, inputs = source
    dataset = build(tmp_path / "cache", source)
    reference, = manifest.references
    gtf = tmp_path / "annotation.gtf"
    gtf.write_bytes(gzip.decompress(inputs["annotation.gtf.gz"].read_bytes()))
    genome = Genome(reference_name=reference.genome_build, annotation_name="ensembl", annotation_version=110,
                    gtf_path_or_url=str(gtf), cache_directory_path=str(tmp_path / "genome"))
    genome.species = reference.species  # Explicit taxonomy for a local Genome.
    genome.index()
    observation = next(dataset.measurements())
    gene = genome.gene_by_id(observation.target.identifier)
    target = ScopedTarget.from_pyensembl(gene, kind="gene", namespace="ensembl", context=reference)
    assert target == observation.target
    assert len(list(dataset.measurements(target=target))) == 2
    with pytest.raises(ValueError, match="annotation_version"):
        ScopedTarget.from_pyensembl(gene, kind="gene", namespace="ensembl", context=replace(reference, annotation_version="111"))
    other_species = "mouse" if reference.taxon_id == 9615 else "dog"
    other = DatasetManifest.from_dict(json.loads((ROOT / f"{other_species}-manifest.json").read_text())).references[0]
    with pytest.raises(ValueError, match="Reference context"):
        list(dataset.measurements(target=replace(target, reference=other)))


@pytest.mark.parametrize("change,error", (("taxon", "taxon"), ("biosample", "BioSample"),
                                          ("duplicate", "Duplicate GEO"), ("assembly", "assembly")))
def test_metadata_and_annotation_disagreements_fail(source, tmp_path, change, error):
    manifest, _, inputs = source
    if change == "assembly":
        asset = "annotation.gtf.gz"
        text = gzip.decompress(inputs[asset].read_bytes()).decode().replace(
            "#!genome-build " + manifest.references[0].genome_build, "#!genome-build wrong")
    else:
        sid = "1" if manifest.processing["native_species"] == "dog" else "1345"
        asset = "GSE325525_family.soft.gz" if sid == "1" else "GSE325523_family.soft.gz"
        text = gzip.decompress(inputs[asset].read_bytes()).decode()
        if change == "taxon":
            text = text.replace("!Sample_taxid_ch1 = " + str(manifest.references[0].taxon_id), "!Sample_taxid_ch1 = 9606")
        elif change == "biosample":
            text = text.replace("BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN", "BioSample: https://www.ncbi.nlm.nih.gov/biosample/WRONG")
        else:
            text += text
    altered = tmp_path / "altered.gz"
    altered.write_bytes(gzip.compress(text.encode(), mtime=0))
    with pytest.raises(ValueError, match=error):
        import_rna(dict(inputs, **{asset: altered}), tmp_path / "bad.sqlite", manifest)
    assert not (tmp_path / "bad.sqlite").exists()


@pytest.mark.parametrize("change,error", (
    ("negative", "nonnegative integers"), ("fractional", "nonnegative integers"),
    ("nan", "nonnegative integers"), ("gene", "Unmapped/duplicate"),
    ("duplicate", "Unmapped/duplicate"), ("sample", "sample/tissue"),
    ("tissue", "sample/tissue"), ("missing", "gene count|reconciliation"),
))
def test_source_changes_fail_without_publishing_partial_output(source, tmp_path, change, error):
    manifest, _, inputs = source
    reader = csv.DictReader(io.StringIO(gzip.decompress(inputs["counts.csv.gz"].read_bytes()).decode()))
    rows = list(reader)
    index = next(i for i, row in enumerate(rows) if row["species"] == manifest.processing["native_species"])
    row = rows[index]
    if change in ("negative", "fractional", "nan"):
        row["rcnt"] = {"negative": "-1", "fractional": "1.5", "nan": "NaN"}[change]
    elif change == "gene":
        row["egid"] = "not-in-Ensembl-110"
    elif change == "duplicate":
        rows.insert(index + 1, dict(row))
    elif change in ("sample", "tissue"):
        row["sid" if change == "sample" else "tnm"] = "unknown"
    else:
        rows.pop(index)
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    altered = tmp_path / "altered.csv.gz"
    altered.write_bytes(gzip.compress(text.getvalue().encode(), mtime=0))
    output = tmp_path / "bad.sqlite"
    with pytest.raises(ValueError, match=error):
        import_rna(dict(inputs, **{"counts.csv.gz": altered}), output, manifest)
    assert not output.exists()


def test_native_rna_cli_build_and_query(source, tmp_path, capsys):
    from test_measurements import definition

    manifest, _, inputs = source
    dataset = build(tmp_path / "first", source)
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definition(dataset)))
    prefix = ["datasets", "--registry", str(config), "--cache-root", str(tmp_path / "cli")]
    overrides = [arg for name, path in inputs.items() for arg in ("--source-path", f"{name}={path}")]
    run(prefix + ["build", manifest.dataset_id, "--quiet"] + overrides + ["--json"])
    assert json.loads(capsys.readouterr().out)["counts"]["measurement"] == manifest.coverage["rows"]
    # Published-artifact queries are a distinct cache contract from local builds.
    run(prefix + ["download", manifest.dataset_id, "--source-path", str(dataset.path), "--quiet"] + ["--json"])
    capsys.readouterr()
    run(prefix + ["query", manifest.dataset_id, "--quantity", "rna_abundance", "--unit", "counts"] + ["--format", "jsonl"])
    records = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(records) == manifest.coverage["rows"]
    registry = DatasetRegistry(definition(dataset))
    assert registry.search_references(manifest.dataset_id) == ()


def test_packaged_native_release_discovery_is_offline(monkeypatch):
    import requests

    from orfdata.sources.pfizer_rna import rna_manifest

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Discovery contacted network"))
    registry = DatasetRegistry()
    for species, alias, taxon, samples, genes, zeros in (
        ("dog", "canis_lupus_familiaris", 9615, 434, 30941, 4218202),
        ("mouse", "mus_musculus", 10090, 606, 56051, 18019555),
    ):
        manifest = rna_manifest(species)
        assert rna_manifest(alias) == manifest
        assert manifest in registry.list(species=species, assay="rna_seq")
        assert manifest.references[0].taxon_id == taxon
        assert manifest.coverage["rows"] == samples * genes
        assert manifest.coverage["zeros"] == zeros
        assert manifest.coverage["samples"] == samples and manifest.coverage["genes"] == genes
        assert manifest.coverage["human_native_observations"] == 0
        assert manifest.coverage["orf_specific_abundance"] is False
        metadata = registry.download_metadata(manifest.dataset_id)
        assert 0 < metadata["distribution"]["download_size"] < metadata["distribution"]["installed_size"]
        assert metadata["assets"]["dataset.sqlite"]["decompress"] is True
    assert not any(m.dataset_id == "pfizer-atlas-human-rna"
                   for m in registry.list(species="human", assay="rna_seq"))
    with pytest.raises(ValueError, match="human counts are not deposited"):
        rna_manifest("human")
