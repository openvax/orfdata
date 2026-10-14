"""Managed catalogues preserve the existing reference objects and source scope."""

import csv
import gzip
import hashlib
import io
import json
from dataclasses import replace
from pathlib import Path

import pytest
from datacache import FileValidationError

from orfdata import (
    CatalogManifest,
    CatalogRegistry,
    EvidenceSummary,
    ORFCatalog,
    ORFEvidence,
    ORFModel,
    ORFRecord,
    ReferenceContext,
    ScopedTarget,
    SourceAsset,
    build_catalog,
    load_catalog,
)
from orfdata.cli import run

SPECIES = ("homo_sapiens", "mus_musculus", "canis_lupus_familiaris", "saccharomyces_cerevisiae")


def setup_release(tmp_path, species="mus_musculus", *, summarize=False):
    # Native reference slices test annotation identity; no experimental support is invented.
    data = json.loads((Path(__file__).parent / "data/references" / (species + ".json")).read_text())
    fields = data["model"]
    fields["reference_context"] = reference = ReferenceContext(**fields["reference_context"])
    model = ORFModel(**fields)
    summaries = (EvidenceSummary(
        "synthetic", "1", ScopedTarget(reference, "protein", "fixture", "synthetic"),
        hashlib.sha256(model.protein_sequence.encode()).hexdigest(), "ms", "support count", None,
        "Synthetic unknown count, not experimental data", "protein_sequence",
        "https://example.invalid/fixture", "a" * 64),) if summarize else ()
    record = ORFRecord.from_model(model, ORFEvidence(summaries=summaries))
    catalog = ORFCatalog([record], dict(genome_build=record.genome_build, taxon_id=record.taxon_id,
        species=species, coordinates="1-based-inclusive", dataset_version="native-fixture",
        sources=[dict(id="native-reference-slice", scope="annotation fixture, no new experimental evidence")]))
    source = tmp_path / (species + ".json")
    source.write_text(json.dumps(catalog.to_dict()))
    asset = SourceAsset("source.json", "https://example.invalid/" + source.name,
                        hashlib.sha256(source.read_bytes()).hexdigest(), source.stat().st_size,
                        "Source-backed reference fixture; upstream provenance in tests/data/references")
    manifest = CatalogManifest(species, "1", "Native reference fixture", (reference,), (asset,),
        ("https://github.com/openvax/orfdata/tree/main/tests/data/references",),
        dict(records=1, distinct_protein_sequences=1, evidence_summaries=len(summaries),
             experimental_evidence=False),
        dict(importer="catalog-copy-v1", catalog_asset="source.json"))
    return manifest, source, catalog


def definition(manifest, artifact=None):
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest() if artifact else "a" * 64
    size = artifact.stat().st_size if artifact else 1
    return {manifest.catalog_id: dict(default_version=manifest.version, versions={manifest.version: dict(
        manifest=manifest.to_dict(), distribution=dict(access_terms="fixture", redistributable=False),
        assets={"catalog.json.gz": dict(url="https://example.invalid/catalog.json.gz", sha256=digest, size=size)})})}


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    import requests

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Unexpected network access"))


@pytest.mark.parametrize("species", SPECIES)
def test_real_species_build_download_open_and_legacy_query(tmp_path, species):
    manifest, source, original = setup_release(tmp_path, species)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    name = manifest.catalog_id
    assert registry.list(species=species) == (manifest,)
    assert registry.list(taxon_id=manifest.taxon_id + 1000000) == ()
    assert registry.inspect(name).status == registry.inspect_build(name).status == "missing"
    built = registry.build(name, source_paths={"source.json": source}, download_options={"show_progress": False})
    assert isinstance(built, ORFCatalog) and built.records == original.records
    assert not built.records[0].evidence.has_evidence("ms")
    assert built.query(gene=original.records[0].gene_id) == original.records
    inspected = registry.inspect_build(name)
    assert inspected.status == "available" and inspected.verified is False
    assert inspected.sources["source.json"]["sha256"] == manifest.sources[0].sha256
    artifact = Path(inspected.files["catalog.json.gz"].path)
    config = definition(manifest, artifact)
    installed = CatalogRegistry(config, cache_root=tmp_path / "download")
    installed.download(name, source_paths={"catalog.json.gz": artifact}, show_progress=False)
    assert installed.inspect(name).verified and installed.open(name).to_dict() == built.to_dict()
    assert load_catalog(source).records == original.records
    source.unlink()  # Neither local query path depends on original acquisition files.
    assert registry.open_build(name).records == installed.open(name).records == original.records


def test_catalogue_cli_reuses_evidence_locus_and_export_options(tmp_path, capsys):
    manifest, source, original = setup_release(tmp_path)
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definition(manifest)))
    prefix = ["catalogs", "--registry", str(config), "--cache-root", str(tmp_path / "cache")]
    run(prefix + ["list", "--species", "mouse"])
    assert json.loads(capsys.readouterr().out)[0]["catalog_id"] == manifest.catalog_id
    run(prefix + ["build", manifest.catalog_id, "--source-path", "source.json=" + str(source)])
    assert json.loads(capsys.readouterr().out)["status"] == "available"
    run(prefix + ["inspect", manifest.catalog_id, "--built"])
    assert json.loads(capsys.readouterr().out)["verified"] is False
    record = original.records[0]
    query = prefix + ["query", manifest.catalog_id, "--built"]
    run(query + ["--orf-id", record.orf_id, "--format", "fasta"])
    assert record.protein_sequence in capsys.readouterr().out
    run(query + ["--evidence", "ms", "--require-resolved"])
    assert json.loads(capsys.readouterr().out)["records"] == []
    run(query + ["--position", f"{record.contig}:{record.coding_blocks[0][0]}", "--genome-build", record.genome_build])
    assert json.loads(capsys.readouterr().out)["records"][0]["orf_id"] == record.orf_id
    with pytest.raises(SystemExit):
        run(query + ["--catalog", str(source)])


@pytest.mark.parametrize("summarize", (False, True))
def test_catalogue_cli_tsv_exports_aggregate_evidence(tmp_path, capsys, summarize):
    # Regression for #58 through the managed-catalogue query path.
    manifest, source, original = setup_release(tmp_path, summarize=summarize)
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definition(manifest)))
    prefix = ["catalogs", "--registry", str(config), "--cache-root", str(tmp_path / "cache")]
    run(prefix + ["build", manifest.catalog_id, "--source-path", "source.json=" + str(source)])
    capsys.readouterr()
    run(prefix + ["query", manifest.catalog_id, "--built", "--format", "tsv"])
    row, = csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t")
    summaries = json.loads(row["summaries_json"])
    assert summaries == [s.to_dict() for s in original.records[0].summaries]
    assert [s["value"] for s in summaries] == ([None] if summarize else [])


def test_changed_source_or_unqualified_record_cannot_publish_a_build(tmp_path):
    manifest, source, original = setup_release(tmp_path)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    changed = tmp_path / "changed.json"
    changed.write_text(source.read_text() + " ")
    with pytest.raises(FileValidationError):
        registry.build(manifest.catalog_id, source_paths={"source.json": changed})
    assert registry.inspect_build(manifest.catalog_id).status == "missing"

    def wrong_content(inputs, output):
        data = original.to_dict()
        data["metadata"]["catalog_manifest"] = manifest.to_dict()
        data["records"][0].pop("reference_context")
        Path(output).write_bytes(gzip.compress(json.dumps(data).encode(), mtime=0))

    with pytest.raises(ValueError, match="declared complete genomic reference"):
        build_catalog(manifest, wrong_content, cache_root=tmp_path / "cache", source_paths={"source.json": source})
    assert registry.inspect_build(manifest.catalog_id).status == "missing"


def test_artifact_integrity_and_whole_manifest_mismatch(tmp_path):
    manifest, source, _ = setup_release(tmp_path)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "build")
    registry.build(manifest.catalog_id, source_paths={"source.json": source})
    artifact = Path(registry.inspect_build(manifest.catalog_id).files["catalog.json.gz"].path)
    config = definition(manifest, artifact)
    altered = tmp_path / "bad.gz"
    altered.write_bytes(artifact.read_bytes() + b"bad")
    installed = CatalogRegistry(config, cache_root=tmp_path / "download")
    with pytest.raises(FileValidationError):
        installed.download(manifest.catalog_id, source_paths={"catalog.json.gz": altered})
    assert installed.inspect(manifest.catalog_id).status == "missing"
    other_manifest = replace(manifest, title="Different declared provenance")
    installed = CatalogRegistry(definition(other_manifest, artifact), cache_root=tmp_path / "other")
    installed.download(manifest.catalog_id, source_paths={"catalog.json.gz": artifact})
    with pytest.raises(ValueError, match="manifest disagrees"):
        installed.open(manifest.catalog_id)


def test_registry_and_manifest_are_immutable_and_identity_checked(tmp_path):
    manifest, _, _ = setup_release(tmp_path)
    data = definition(manifest)
    registry = CatalogRegistry(data, cache_root=tmp_path / "cache")
    data[manifest.catalog_id]["versions"]["1"]["manifest"]["title"] = "mutated"
    assert registry.info(manifest.catalog_id) == manifest
    assert CatalogManifest.from_dict(manifest.to_dict()) == manifest
    with pytest.raises(TypeError):
        manifest.coverage["records"] = 0
    for changes in (dict(catalog_id="../escape"), dict(version="a/b"), dict(coverage={"records": True}),
                    dict(coverage={"records": 1, "evidence_summaries": False}), dict(processing={}),
                    dict(references=(manifest.references[0], manifest.references[0]))):
        with pytest.raises((ValueError, TypeError)):
            replace(manifest, **changes)
    for name, version in (("unknown", None), (manifest.catalog_id, "unknown")):
        with pytest.raises(ValueError, match="Unknown"):
            registry.info(name, version)
    with pytest.raises(ValueError, match="disagree"):
        registry.list(species="mouse", taxon_id=9606)
    with pytest.raises(ValueError, match="No packaged importer"):
        CatalogRegistry(definition(replace(manifest, processing={"importer": "arbitrary.module:code"}))).build(manifest.catalog_id)


def test_force_build_preserves_readers_and_stays_separate_from_downloads(tmp_path):
    manifest, source, _ = setup_release(tmp_path)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    name = manifest.catalog_id
    registry.build(name, source_paths={"source.json": source})
    first = Path(registry.inspect_build(name).files["catalog.json.gz"].path)
    registry.build(name, source_paths={"source.json": source}, force=True)
    second = Path(registry.inspect_build(name).files["catalog.json.gz"].path)
    assert first != second and first.exists() and second.exists()
    assert registry.inspect(name).status == "missing"
    with pytest.raises(FileNotFoundError):
        registry.open(name)


def test_corrupted_local_build_fails_receipt_verification(tmp_path):
    manifest, source, _ = setup_release(tmp_path)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    registry.build(manifest.catalog_id, source_paths={"source.json": source})
    path = Path(registry.inspect_build(manifest.catalog_id).files["catalog.json.gz"].path)
    path.write_bytes(b"corrupted")
    assert registry.inspect_build(manifest.catalog_id).status != "available"
    with pytest.raises(FileValidationError):
        registry.open_build(manifest.catalog_id)


def test_same_catalogue_name_and_version_do_not_share_species_caches(tmp_path):
    paths = []
    for species in SPECIES:
        manifest, source, original = setup_release(tmp_path, species)
        manifest = replace(manifest, catalog_id="colliding-name")
        registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
        built = registry.build(manifest.catalog_id, source_paths={"source.json": source})
        assert built.taxon_id == original.taxon_id == manifest.taxon_id
        paths.append(registry.inspect_build(manifest.catalog_id).files["catalog.json.gz"].path)
    assert len(set(paths)) == 4


def test_coverage_and_unregistered_reference_are_not_accepted(tmp_path):
    manifest, source, _ = setup_release(tmp_path)
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    built = registry.build(manifest.catalog_id, source_paths={"source.json": source})
    for new_manifest in (replace(manifest, coverage={"records": 2}),
                         replace(manifest, coverage={"records": 1, "distinct_protein_sequences": 2}),
                         replace(manifest, references=(replace(manifest.references[0], annotation_version="future"),))):
        metadata = built.metadata
        metadata["catalog_manifest"] = new_manifest.to_dict()
        with pytest.raises(ValueError, match="count|reference"):
            new_manifest.validate_catalog(ORFCatalog(built.records, metadata))


def test_copy_importer_does_not_overwrite_an_existing_release_manifest(tmp_path):
    manifest, source, original = setup_release(tmp_path)
    data = original.to_dict()
    data["metadata"]["catalog_manifest"] = manifest.to_dict()
    source.write_text(json.dumps(data))
    source_asset = replace(manifest.sources[0], size=source.stat().st_size,
                           sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    manifest = replace(manifest, sources=(source_asset,))
    registry = CatalogRegistry(definition(manifest), cache_root=tmp_path / "cache")
    with pytest.raises(ValueError, match="cannot replace release provenance"):
        registry.build(manifest.catalog_id, source_paths={"source.json": source})
    assert registry.inspect_build(manifest.catalog_id).status == "missing"
