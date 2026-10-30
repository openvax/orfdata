"""Actual three-species source slices, including zeros, NA and group ambiguity."""

import csv
import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry, QuantitativeDataset, build_dataset
from orfdata.cli import run
from orfdata.sources.pfizer import import_protein, protein_manifest, protein_sample_id

FIXTURES = Path(__file__).parent / "data/pfizer"


@pytest.fixture(params=("human", "mouse", "dog"))
def source(request):
    root = FIXTURES / request.param
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    expected = json.loads((root / "expected.json").read_text())
    inputs = {asset.asset_id: root / asset.asset_id for asset in manifest.sources}
    return manifest, inputs, expected


def build(tmp_path, source, **kwargs):
    manifest, inputs, _ = source
    return build_dataset(manifest, lambda sources, output: import_protein(sources, output, manifest),
                         cache_root=tmp_path, source_paths=inputs, **kwargs)


def test_real_native_values_and_lineage_roundtrip(tmp_path, source):
    dataset = build(tmp_path, source)
    manifest, _, expected = source
    observations = {m.measurement_id: m for m in dataset.measurements()}
    assert len(observations) == len(expected["rows"]) * 2
    for line, row in enumerate(expected["rows"], 2):
        for column, unit in (("IBAQ", "iBAQ"), ("rIBAQ", "rIBAQ")):
            observation = observations[f"row:{line}:{column}"]
            assert observation.value == (None if row[column] == "NA" else float(row[column]))
            assert observation.status == ("missing" if row[column] == "NA" else "measured")
            assert observation.target.kind == "protein_group"
            assert observation.target.identifier == row["ProteinGroup"]
            assert observation.target.reference == manifest.references[0]
            assert observation.target.reference.genome_build is None
            assert observation.sample_id == protein_sample_id(row)
            assert observation.unit == unit
            assert observation.provenance["native_Abundance"] == row["Abundance"]
            assert observation.provenance["native_sumIBAQ"] == row["sumIBAQ"]
    assert any(m.value == 0 and m.status == "measured" for m in observations.values())
    assert any(m.value is None and m.status == "missing" for m in observations.values())
    assert any(m.provenance["native_Abundance"] == "-Inf" for m in observations.values())
    assert dataset.counts["sample"] == manifest.coverage["source_counts"]["assays"]
    assert list(dataset.pairings()) == []
    for sample in dataset.samples():
        assert sample.attributes["pooling"] == "pooled"
        assert sample.pool_members == () and sample.biological_replicate is None
    for assay in dataset.assays():
        assert assay.bridge == "131C" and assay.technical_replicate == "1"
        assert assay.calibration["absolute_molecule_counts"] is False


def test_source_rebuild_is_offline_reusable_and_hash_verified(tmp_path, source):
    first = build(tmp_path, source)
    assert build(tmp_path, source).path == first.path
    corrupted = tmp_path / "wrong.tsv"
    corrupted.write_text("wrong source\n")
    manifest, paths, _ = source
    with pytest.raises(Exception, match="(?i)(hash|sha|size|integrity)"):
        build_dataset(manifest, lambda i, o: pytest.fail("unverified inputs reached importer"),
                      cache_root=tmp_path / "wrong-cache", source_paths=dict(paths, **{"protein.tsv": corrupted}))


def test_group_representative_lookup_never_allocates_individual_abundance(tmp_path, source):
    dataset = build(tmp_path, source)
    mappings = list(dataset.reference_mappings())
    assert len(mappings) == source[0].coverage["source_counts"]["groups"]
    assert all(m.status in ("unresolved", "unmapped") for m in mappings)
    found = next(m for m in mappings if m.targets)
    assert found.validation["representative_search_record"]["sequence_sha256"]
    assert not found.sequence_validated
    target = found.targets[0]
    assert list(dataset.measurements(target=target)) == []
    assert list(dataset.mapped_measurements(target)) == []
    views = list(dataset.mapped_measurements(target, mapping_statuses=["unresolved"]))
    assert views and all(v.attribution == "group_membership_context" for v in views)
    assert all(v.measurement.target.kind == "protein_group" for v in views)


def test_manifest_species_aliases_and_real_source_coverage():
    dog = protein_manifest("canis_lupus_familiaris")
    assert dog == protein_manifest("dog")
    assert dog.coverage["source_counts"]["rows"] == 313263
    assert dog.coverage["source_counts"]["groups"] == 9929
    assert dog.coverage["source_counts"]["absent_matrix_cells"] == 64039
    assert dog.coverage["source_counts"]["IBAQ:missing"] == 8
    assert protein_manifest("human").coverage["source_counts"]["rows"] == 325702
    assert protein_manifest("mouse").coverage["source_counts"]["rows"] == 615991
    with pytest.raises(ValueError, match="human, mouse and dog"):
        protein_manifest("yeast")


def test_stable_consumer_sample_id():
    row = dict(Run="Pool4_1", Channel="128_C", BioReplicate="RBI-0000137288", Condition="GI_Colon",
               TechRepMixture="1", Mixture="Pool4", **{"sample.type": "sample"})
    assert protein_sample_id(row) == "openvax:protein-sample:v1:sha256:042357dbd48da65f98474732bdc70a66111b4112ff9b295da88f386d3c9d542f"


@pytest.mark.parametrize("change,error", (("duplicate", "Duplicate"), ("representative", "Conflicting"),
                                          ("infinity", "finite"), ("negative", "nonnegative"),
                                          ("role", "bridge"), ("count", "reconciliation")))
def test_source_schema_and_coverage_changes_fail(tmp_path, change, error):
    root = FIXTURES / "dog"
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    rows = list(csv.DictReader((root / "protein.tsv").open(), delimiter="\t"))
    if change == "duplicate":
        rows.append(rows[0])
    elif change == "representative":
        other = dict(rows[0], Protein="other", Run="different-run")
        rows.append(other)
    elif change in ("infinity", "negative"):
        rows[0]["IBAQ"] = "inf" if change == "infinity" else "-1"
    elif change == "role":
        rows[0]["sample.type"] = "bridge"
    elif change == "count":
        rows.pop()
    path = tmp_path / "input.tsv"
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match=error):
        import_protein({"protein.tsv": path, "search.fasta": root / "search.fasta"},
                       tmp_path / "invalid.sqlite", manifest)
    assert not (tmp_path / "invalid.sqlite").exists()


def test_fasta_taxon_and_reference_pin_are_required(tmp_path, source):
    manifest, inputs, _ = source
    wrong = replace(manifest, references=(replace(manifest.references[0], sequence_sha256="a" * 64),))
    with pytest.raises(ValueError, match="fingerprints"):
        import_protein(inputs, tmp_path / "bad.sqlite", wrong)
    path = tmp_path / "wrong.fasta"
    path.write_text(inputs["search.fasta"].read_text().replace(f"OX={manifest.references[0].taxon_id}", "OX=999999"))
    with pytest.raises(ValueError, match="taxon"):
        import_protein(dict(inputs, **{"search.fasta": path}), tmp_path / "bad.sqlite", manifest)


def test_registry_build_and_cli_use_same_verified_source_workflow(tmp_path, source, capsys):
    from test_measurements import definition

    dataset = build(tmp_path / "first", source)
    definitions = definition(dataset)
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(definitions))
    cache = tmp_path / "cli-cache"
    options = ["datasets", "--registry", str(registry_path), "--cache-root", str(cache),
               "build", source[0].dataset_id, "--quiet"]
    for key, path in source[1].items():
        options.extend(["--source-path", f"{key}={path}"])
    run(options + ["--json"])
    result = json.loads(capsys.readouterr().out)
    rebuilt = QuantitativeDataset(result["path"])
    assert rebuilt.manifest == dataset.manifest
    registry = DatasetRegistry(definitions, cache_root=cache)
    assert registry.build(source[0].dataset_id, source_paths=source[1]).path == rebuilt.path
    assert registry.inspect(source[0].dataset_id).status == "missing"  # Builds do not overwrite installed releases.


def test_compressed_artifact_installs_verified_database_offline(tmp_path, source):
    from test_measurements import definition

    dataset = build(tmp_path / "built", source)
    definitions = definition(dataset)
    compressed = tmp_path / "dataset.sqlite.gz"
    compressed.write_bytes(gzip.compress(dataset.path.read_bytes(), mtime=0))
    artifact = definitions[dataset.manifest.dataset_id]["versions"][dataset.manifest.version]["assets"]["dataset.sqlite"]
    artifact.update(url=compressed.as_uri(), decompress=True)
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "installed")
    registry.download(dataset.manifest.dataset_id, show_progress=False,
                      source_paths={"dataset.sqlite": compressed})
    installed = registry.open(dataset.manifest.dataset_id)
    assert installed.manifest == dataset.manifest
    assert hashlib.sha256(installed.path.read_bytes()).hexdigest() == artifact["sha256"]
    assert list(installed.measurements()) == list(dataset.measurements())


def test_packaged_registry_has_parallel_source_native_releases():
    registry = DatasetRegistry()
    for species in ("human", "mouse", "dog"):
        manifest = protein_manifest(species)
        assert manifest in registry.list(species=species, assay="digest_ms")
        assert registry.info(manifest.dataset_id) == manifest
        metadata = registry.download_metadata(manifest.dataset_id)
        assert metadata["assets"]["dataset.sqlite"]["decompress"] is True
        assert metadata["distribution"]["redistributable"] is True
    assert registry.download_metadata("pfizer-atlas-dog-protein")["distribution"]["download_size"] < 100_000_000
