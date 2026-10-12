"""Full search-space identity is independent of reported abundance and detection."""

import hashlib
import json
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from orfdata import (
    DatasetManifest,
    DatasetRegistry,
    ScopedTarget,
    SearchMembership,
    SearchReference,
    SearchReferenceSpec,
    SearchSequence,
)
from orfdata.cli import run
from orfdata.sources.pfizer import import_protein


@pytest.fixture(params=("human", "mouse", "dog"))
def source(request, tmp_path):
    from test_measurements import definition

    root = Path(__file__).parent / "data/pfizer" / request.param
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    dataset = import_protein({a.asset_id: root / a.asset_id for a in manifest.sources}, tmp_path / "data.sqlite", manifest)
    return root, dataset, definition(dataset)


def test_discovery_download_and_offline_lookup(source, tmp_path, monkeypatch):
    root, dataset, definitions = source
    cache = tmp_path / "cache"
    registry = DatasetRegistry(definitions, cache_root=cache)
    name = dataset.manifest.dataset_id
    spec, = registry.search_references(name)
    assert SearchReferenceSpec.from_dict(spec.to_dict()) == spec
    assert not cache.exists()
    assert registry.inspect_reference(name).status == "missing"
    with pytest.raises(FileNotFoundError):
        registry.open_reference(name)
    paths = registry.download_reference(name, source_path=root / "search.fasta", show_progress=False)
    assert registry.inspect(name).status == "missing"  # No abundance artifact was downloaded.
    assert registry.inspect_reference(name).status == "available"
    reference = registry.open_reference(name)
    assert reference.path == Path(paths["search.fasta"])
    assert len(reference) == spec.expected_sequences
    assert reference.source.sha256 == spec.reference.sequence_sha256
    entry = next(reference.entries())
    assert entry.sequence_sha256 == hashlib.sha256(entry.sequence.encode()).hexdigest()
    assert entry.source_taxon_id == spec.reference.taxon_id
    assert entry.decoy is None and entry.control is None
    assert SearchSequence.from_dict(entry.to_dict()) == entry
    present = reference.membership(entry.target)
    assert present.status == "present" and present.sequence == entry
    assert SearchMembership.from_dict(present.to_dict()) == present
    missing = replace(entry.target, identifier="not-in-this-search")
    assert reference.membership(missing).status == "not_in_search"
    assert reference.membership(replace(entry.target, namespace="different")).status == "unknown"
    assert reference.membership(replace(entry.target, kind="protein_group")).status == "unknown"
    assert reference.membership(replace(entry.target, kind="orf")).status == "unknown"
    with pytest.raises(ValueError, match="different reference"):
        reference.membership(replace(entry.target, reference=replace(entry.target.reference, sequence_sha256="a" * 64)))
    with pytest.raises(FrozenInstanceError):
        reference.spec = replace(spec, namespace="changed")
    with pytest.raises(TypeError):
        reference._entries["new"] = entry
    import requests
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Offline lookup contacted network"))
    assert registry.open_reference(name).get(entry.target.identifier) == entry
    old_path = reference.path
    registry.download_reference(name, source_path=root / "search.fasta", force=True, show_progress=False)
    assert old_path.exists() and reference.get(entry.target.identifier) == entry


def test_corrupt_source_and_mismatched_scope_rejected(source, tmp_path):
    root, dataset, definitions = source
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    name = dataset.manifest.dataset_id
    spec, = registry.search_references(name)
    corrupt = tmp_path / "corrupt.fasta"
    corrupt.write_bytes(b"bad")
    with pytest.raises(Exception, match="(?i)(hash|size|sha|integrity)"):
        registry.download_reference(name, source_path=corrupt, show_progress=False)
    assert registry.inspect_reference(name).status == "missing"
    with pytest.raises(ValueError, match="SHA-256"):
        SearchReference(corrupt, spec)
    with pytest.raises(ValueError, match="sequence count"):
        SearchReference(root / "search.fasta", replace(spec, expected_sequences=999))
    wrong_reference = replace(spec.reference, taxon_id=999999, species="test_species")
    with pytest.raises(ValueError, match="taxon"):
        SearchReference(root / "search.fasta", replace(spec, reference=wrong_reference))
    with pytest.raises(ValueError, match="exact source FASTA"):
        replace(spec, reference=replace(spec.reference, sequence_sha256="a" * 64))
    with pytest.raises(ValueError, match="No declared"):
        registry.open_reference(name, reference="wrong-reference")


def test_reference_cli(source, tmp_path, capsys):
    root, dataset, definitions = source
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definitions))
    prefix = ["datasets", "--registry", str(config), "--cache-root", str(tmp_path / "cache")]
    name = dataset.manifest.dataset_id
    run(prefix + ["references", name])
    spec, = json.loads(capsys.readouterr().out)
    assert spec["source"]["sha256"] == dataset.manifest.references[0].sequence_sha256
    run(prefix + ["download-reference", name, "--source-path", str(root / "search.fasta"), "--quiet"])
    assert Path(json.loads(capsys.readouterr().out)["search.fasta"]).exists()
    identifier = next(iter(root.joinpath("search.fasta").read_text().splitlines())).split("|")[1]
    run(prefix + ["query-reference", name, "--target", identifier])
    result = json.loads(capsys.readouterr().out)
    assert result["membership"]["status"] == "present"
    assert result["source"]["source"]["license_url"]


def test_native_reference_declaration_excludes_calibration_standards():
    registry = DatasetRegistry()
    assert registry.search_references("copy-yeast-abundance") == ()
    with pytest.raises(ValueError, match="No declared native"):
        registry.open_reference("copy-yeast-abundance")
    for species, count in (("human", 20600), ("mouse", 21989), ("dog", 20649)):
        spec, = registry.search_references(f"pfizer-atlas-{species}-protein")
        assert spec.expected_sequences == count


def test_reported_and_unreported_native_search_entries(tmp_path):
    from test_measurements import definition

    root = Path(__file__).parent / "data/search_reference"
    expected = json.loads((root / "expected.json").read_text())
    native = SearchReferenceSpec.from_dict(expected["source"])
    data = (root / "dog.fasta").read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    assert digest == expected["slice_sha256"]
    reference = replace(native.reference, sequence_sha256=digest)
    asset = replace(native.source, sha256=digest, size=len(data))
    spec = replace(native, reference=reference, source=asset, expected_sequences=2)
    search = SearchReference(root / "dog.fasta", spec)
    protein_root = Path(__file__).parent / "data/pfizer/dog"
    manifest = DatasetManifest.from_dict(json.loads((protein_root / "manifest.json").read_text()))
    manifest = replace(manifest, references=(reference,),
                       sources=tuple(asset if a.asset_id == "search.fasta" else a for a in manifest.sources),
                       coverage=dict(manifest.coverage, search_sequences=2))
    dataset = import_protein({"protein.tsv": protein_root / "protein.tsv", "search.fasta": search.path},
                             tmp_path / "dataset.sqlite", manifest)
    observed = search.get(expected["reported_representative"])
    unreported = search.get(expected["no_reported_representative"])
    assert observed and unreported
    assert search.membership(unreported.target).status == "present"
    assert unreported.sequence.startswith("MAKINTQVSHPVVRTMDRDL")
    assert list(dataset.measurements(target=unreported.target)) == []
    assert list(dataset.mapped_measurements(unreported.target, mapping_statuses=["unresolved"])) == []
    views = list(dataset.mapped_measurements(observed.target, mapping_statuses=["unresolved"]))
    assert views and all(v.attribution == "group_membership_context" for v in views)
    assert search.membership(views[0].measurement.target).status == "unknown"

    # A changed reference pin cannot reuse another source's cache, even when a
    # local registry uses the same dataset name/version during fixture testing.
    first = DatasetRegistry(definition(dataset), cache_root=tmp_path / "cache")
    first.download_reference(manifest.dataset_id, source_path=search.path, show_progress=False)
    original = DatasetManifest.from_dict(json.loads((protein_root / "manifest.json").read_text()))
    other = definition(dataset)
    other[manifest.dataset_id]["versions"][manifest.version]["manifest"] = original.to_dict()
    second = DatasetRegistry(other, cache_root=tmp_path / "cache")
    assert second.inspect_reference(manifest.dataset_id).status == "missing"
    second.download_reference(manifest.dataset_id, source_path=protein_root / "search.fasta", show_progress=False)
    assert first.open_reference(manifest.dataset_id).get(unreported.target.identifier) == unreported
    assert second.open_reference(manifest.dataset_id).get(unreported.target.identifier) is None


@pytest.mark.parametrize("text,error", (
    ("MAK\n", "precedes"),
    (">\nMAK\n", "Empty"),
    (">sp|P1|TEST OX=4932\n", "empty"),
    (">sp|P1|TEST OX=4932\nMAK\n>tr|P1|TEST OX=4932\nMAK\n", "Duplicate"),
    (">sp|P1|TEST OX=4932\nM*K\n", "invalid"),
    (">sp|P1|TEST\nMAK\n", "header"),
    (">sp|P1|TEST OX=9606\nMAK\n", "taxon"),
))
def test_malformed_declared_fasta_rejected(tmp_path, text, error):
    from orfdata import ReferenceContext, SourceAsset

    path = tmp_path / "test.fasta"
    path.write_text(text)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    source = SourceAsset("test.fasta", "https://example.invalid/test.fasta", digest,
                         path.stat().st_size, "Synthetic validation fixture")
    reference = ReferenceContext(4932, "yeast", None, "uniprot", None, sequence_sha256=digest)
    spec = SearchReferenceSpec("test", "1", reference, source, "uniprot", "uniprot-fasta")
    with pytest.raises(ValueError, match=error):
        SearchReference(path, spec)


def test_arbitrary_taxon_is_scoped_without_mammal_defaults(tmp_path):
    from orfdata import ReferenceContext, SourceAsset

    path = tmp_path / "test.fasta"
    path.write_text(">tr|TEST|TEST_ODD Test sequence OX=123456789\nMAK\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    source = SourceAsset("test.fasta", "https://example.invalid/test.fasta", digest,
                         path.stat().st_size, "Synthetic validation fixture; not experimental evidence")
    reference = ReferenceContext(123456789, "synthetic_test_species", None, "uniprot", None, sequence_sha256=digest)
    spec = SearchReferenceSpec("test", "1", reference, source, "uniprot", "uniprot-fasta", 1)
    search = SearchReference(path, spec)
    target = ScopedTarget(reference, "protein", "TEST", "uniprot")
    assert search.membership(target).status == "present"
    assert search.get("TEST").gene_label is None
