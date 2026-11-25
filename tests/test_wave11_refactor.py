"""Tests for Wave 11:
- #201: deduplicate bundled uorf-evidence against orf-evidence
- #145: immutable nested Proteome.layers & evidence_detail metadata with consistent identity
- #161: unified ResourceKind protocol & registry across list/status/info/path/fetch/setup
- #162: consolidate acquisition commands onto data fetch + GDC viral inventory resource
- #163: unified --source-path grammar (ASSET=PATH and bare PATH for single-asset resources)
- #157: content-addressed SHA-256 SourceStore with re-verification, PyEnsembl seeding, and pruning
- #156: resumable downloads for compressed single-file (decompress: true) datasets and build sources
"""

from __future__ import annotations

import gzip
import hashlib
import inspect
import json
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType

import datacache.bundles
import pytest

from orfdata import (
    RESOURCE_KINDS,
    CatalogRegistry,
    DatasetRegistry,
    DataStatus,
    ProteinSequence,
    Proteome,
    ReferenceContext,
    ResourceKind,
    SourceStore,
    data_info,
    data_path,
    fetch_data,
    get_resource_kind,
    list_data,
    load_catalog,
    load_uorf_catalog,
    prune_source_store,
    registered_resource_kinds,
    resolve_source_paths,
    source_store_dir,
)
from orfdata.cli import run
from orfdata.cli_parser import UsageError, parse_source_paths

# --- #201: Deduplicated bundled uorf-evidence ---


def test_uorf_catalog_is_filtered_view_without_duplicate_data_dir():
    assert not files("orfdata").joinpath("data/uorf-evidence").exists()
    full = load_catalog()
    uorfs = load_uorf_catalog()
    assert len(uorfs.records) == 3771
    assert uorfs.records == tuple(r for r in full.records if r.biotype in ("uORF", "uoORF"))


def test_build_module_defaults_to_full_orf_evidence_catalog():
    from orfdata import build as build_mod

    assert inspect.signature(build_mod.build).parameters["all_orfs"].default is True
    assert inspect.signature(build_mod.prepare_sources).parameters["all_orfs"].default is True


# --- #145: Immutable nested Proteome.layers & identity consistency ---


def test_proteome_nested_layer_metadata_is_immutable_repro_145(tmp_path):
    """Exact regression test for #145 plus caller-input aliasing and nested lists."""
    context = ReferenceContext(10090, "mus_musculus", "GRCm39", "ensembl", "115")
    catalog_meta = {"records": 1, "tags": ["a", "b"], "nested": {"k": "v"}}
    layer_input = dict(
        role="host",
        kind="extra_orf",
        source="fixture",
        reference=context.key,
        catalog=catalog_meta,
    )
    proteome = Proteome(context, [], [layer_input])
    before = proteome.identity

    # Mutating caller-supplied nested dicts/lists after construction has no effect.
    catalog_meta["records"] = 99
    catalog_meta["tags"].append("c")
    catalog_meta["nested"]["k"] = "mutated"
    layer_input["source"] = "changed"
    assert proteome.layers[0]["catalog"]["records"] == 1
    assert proteome.layers[0]["catalog"]["tags"] == ("a", "b")
    assert proteome.layers[0]["catalog"]["nested"]["k"] == "v"
    assert proteome.layers[0]["source"] == "fixture"

    # Attempting to mutate nested dicts on proteome.layers raises TypeError.
    assert isinstance(proteome.layers[0], MappingProxyType)
    assert isinstance(proteome.layers[0]["catalog"], MappingProxyType)
    assert isinstance(proteome.layers[0]["catalog"]["nested"], MappingProxyType)
    with pytest.raises(TypeError):
        proteome.layers[0]["catalog"]["records"] = 2
    with pytest.raises(TypeError):
        proteome.layers[0]["catalog"]["nested"]["k"] = "x"
    assert isinstance(proteome.layers[0]["catalog"]["tags"], tuple)

    # save(), load(), summary() JSON serialization, and merge() preserve exact metadata and identity.
    snapshot = tmp_path / "host.json.gz"
    saved = proteome.save(snapshot)
    assert proteome.identity == before
    assert saved == before

    loaded = Proteome.load(snapshot)
    assert loaded.identity == before
    assert loaded.layers == proteome.layers
    with pytest.raises(TypeError):
        loaded.layers[0]["catalog"]["records"] = 3

    summary = proteome.summary()
    assert json.loads(json.dumps(summary))["identity"] == before
    assert proteome.merge(loaded).identity == before


def test_protein_sequence_evidence_detail_is_deep_frozen(tmp_path):
    context = ReferenceContext(10090, "mus_musculus", "GRCm39", "ensembl", "115")
    detail = {"peptides": [{"sequence": "MAAA", "datasets": ["PXD1"]}]}
    entry = ProteinSequence(
        accession="ORF1",
        sequence="MAAA",
        kind="extra_orf",
        taxon_id=10090,
        source="fixture",
        evidence=["ms"],
        evidence_detail=detail,
    )
    proteome = Proteome(context, [entry], [{"role": "host", "kind": "extra_orf", "source": "fixture"}])
    before = proteome.identity

    # Mutating caller dict does not alter frozen evidence_detail or identity.
    detail["peptides"][0]["sequence"] = "MZZZ"
    assert proteome.proteins[0].evidence_detail["peptides"][0]["sequence"] == "MAAA"
    with pytest.raises(TypeError):
        proteome.proteins[0].evidence_detail["peptides"][0]["sequence"] = "MZZZ"

    snapshot = tmp_path / "ev.json.gz"
    assert proteome.save(snapshot) == before
    loaded = Proteome.load(snapshot)
    assert loaded.identity == before
    assert loaded.proteins[0].evidence_detail["peptides"][0]["sequence"] == "MAAA"


# --- #161: Unified ResourceKind protocol and extensible registry ---


def test_every_registered_resource_kind_satisfies_protocol_and_supports_custom_kind(tmp_path, monkeypatch):
    kinds = registered_resource_kinds()
    assert set(kinds) == {"viral", "inventory", "catalog", "dataset", "reference", "genome", "chain", "records"}
    for kind in kinds:
        assert isinstance(get_resource_kind(kind), ResourceKind)

    # Registering a new kind adds it to list_data, data_info, data_path, and fetch_data automatically.
    class _ToyKind:
        kind = "toy"

        def __init__(self):
            self.installed = False
            self.target = tmp_path / "toy.txt"

        def _row(self):
            status = "available" if self.installed else "missing"
            return DataStatus(
                kind="toy",
                name="toy-item",
                version="1",
                title="Toy fixture resource",
                species=("homo_sapiens",),
                taxon_ids=(9606,),
                contexts=(),
                status=status,
                path=str(self.target),
                origin="download" if self.installed else None,
                acquire="download",
                command=None,
                default=True,
                download_size=4,
                installed_size=4,
            )

        def rows(self, *, cache_root=None, verify_files=False):
            return (self._row(),)

        def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None):
            if name != "toy-item":
                raise ValueError(f"Unknown toy {name!r}")
            return self._row()

        def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None):
            return dict(self.inspect(name, version).to_dict(), extra="toy-meta")

        def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
            return resolve_source_paths(source_paths, ("toy.txt",)), {}

        def acquire(
            self, name, version=None, *, cache_root=None, force=False,
            source_paths=None, accept_terms=False, download_options=None, options=None,
        ):
            self.target.write_text("data")
            self.installed = True
            return self.target

    toy = _ToyKind()
    monkeypatch.setitem(RESOURCE_KINDS._handlers, "toy", toy)
    assert "toy" in registered_resource_kinds()
    assert [row.name for row in list_data(kind="toy")] == ["toy-item"]
    assert data_info("toy-item")["extra"] == "toy-meta"
    fetched = fetch_data("toy-item", show_progress=False)
    assert fetched["acquired"] == "downloaded" and fetched["installed"] is True
    assert data_path("toy-item") == tmp_path / "toy.txt"


# --- #162 & #163: data fetch consolidation, GDC viral inventory, unified --source-path ---


def test_gdc_viral_inventory_is_listed_and_fetchable_with_unified_source_path(tmp_path, capsys):
    rows = list_data(kind="inventory", cache_root=tmp_path)
    assert [row.name for row in rows] == ["gdc-viruses-inventory"]
    assert rows[0].status == "missing"
    info = data_info("gdc-viruses-inventory", cache_root=tmp_path)
    assert info["source"]["filename"] == "GRCh83.d1.vd1_virus_decoy.txt"
    assert info["download_size"] == info["source"]["size"]

    # Fetch using the packaged fixture bytes via unified bare PATH --source-path.
    fixture = Path(__file__).parent / "data/viral/GRCh83.d1.vd1_virus_decoy.txt.gz"
    local_inv = tmp_path / "local-inv.txt"
    local_inv.write_bytes(gzip.decompress(fixture.read_bytes()))

    run([
        "data", "fetch", "gdc-viruses-inventory",
        "--cache-root", str(tmp_path / "cache1"),
        "--source-path", str(local_inv),
        "--json",
    ])
    fetched = json.loads(capsys.readouterr().out)
    assert fetched[0]["kind"] == "inventory"
    assert fetched[0]["status"] == "available"
    assert fetched[0]["acquired"] == "downloaded"
    assert Path(fetched[0]["path"]).is_file()

    # Also accepts ASSET=PATH syntax through both data fetch and viral inventory download.
    run([
        "viral", "inventory", "download",
        "--cache-root", str(tmp_path / "cache2"),
        "--source-path", f"GRCh83.d1.vd1_virus_decoy.txt={local_inv}",
        "--json",
    ])
    inv_out = json.loads(capsys.readouterr().out)
    assert Path(inv_out["path"]).is_file()


def test_unified_source_path_grammar_bare_and_named():
    assert parse_source_paths([]) is None
    assert parse_source_paths(["/tmp/one.gz"]) == "/tmp/one.gz"
    assert parse_source_paths(["/tmp/one.gz"], assets=("file.gz",)) == {"file.gz": "/tmp/one.gz"}
    assert parse_source_paths(["a=/tmp/a", "b=/tmp/b"]) == {"a": "/tmp/a", "b": "/tmp/b"}

    with pytest.raises(UsageError, match="Resource has 2 asset"):
        parse_source_paths(["/tmp/one.gz"], assets=("a", "b"))
    with pytest.raises(UsageError, match="Cannot combine"):
        parse_source_paths(["/tmp/one.gz", "a=/tmp/a"])
    with pytest.raises(UsageError, match="bare --source-path"):
        parse_source_paths(["/tmp/one.gz", "/tmp/two.gz"])

    assert resolve_source_paths("/tmp/only.fa", ("search.fasta",)) == {"search.fasta": "/tmp/only.fa"}
    with pytest.raises(ValueError, match="one asset"):
        resolve_source_paths("/tmp/only.fa", ("gtf", "transcripts", "proteins"))
    with pytest.raises(ValueError, match="unknown asset"):
        resolve_source_paths({"bogus": "/tmp/x"}, ("search.fasta",))


def test_data_fetch_validates_kind_specific_options_and_supports_genome_fasta(tmp_path, monkeypatch, capsys):
    # Genome-only and viral-only options are rejected for incompatible kinds before any download.
    with pytest.raises(ValueError, match="Unsupported option.*genome_fasta_path_or_url"):
        fetch_data("bundled", genome_fasta_path_or_url="/tmp/dna.fa")
    with pytest.raises(ValueError, match="Unsupported option.*prefer_reference_molecule"):
        fetch_data("bundled", prefer_reference_molecule=False)
    with pytest.raises(ValueError, match="Unsupported option.*prefer_alternate_annotation"):
        fetch_data("t2t-chm13", prefer_alternate_annotation=True, cache_root=tmp_path)

    # data fetch accepts --genome-fasta for a genome resource and installs chromosome DNA.
    from test_genomes import slice_pins

    from orfdata import genomes, load_genome

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    dna = tmp_path / "chromosomes.fa"
    dna.write_text(">MT\nATATGAAGG\n")
    argv = [
        "data", "fetch", "mouse-slice",
        "--cache-root", str(tmp_path / "gcache"),
        "--genome-fasta", str(dna),
        "--json",
    ]
    for role, path in paths.items():
        argv += ["--source-path", f"{role}={path}"]
    run(argv)
    fetched = json.loads(capsys.readouterr().out)
    assert fetched[0]["kind"] == "genome" and fetched[0]["installed"] is True
    opened = load_genome("mouse-slice", cache_root=tmp_path / "gcache", genome_fasta_path_or_url=dna)
    assert opened.sequence("MT", 1, 3) == "ATA"


def test_data_fetch_accepts_source_path_for_chains_and_references(tmp_path, monkeypatch):
    from orfdata import liftover

    chain_file = tmp_path / "toy.chain"
    chain_bytes = b"chain 1 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n"
    chain_file.write_bytes(chain_bytes)
    fake_pin = liftover.PinnedChain(
        name="A_to_B",
        species="homo_sapiens",
        source_build="A",
        target_build="B",
        url="https://example.invalid/toy.chain",
        sha256=hashlib.sha256(chain_bytes).hexdigest(),
        size=len(chain_bytes),
        provider="Fixture",
        terms="Open",
        terms_url="https://example.invalid/terms",
        restricted=False,
    )
    monkeypatch.setattr(liftover, "pinned_chains", lambda species=None: (fake_pin,))

    # Works with bare PATH and with FILENAME=PATH.
    r1 = fetch_data("A_to_B", cache_root=tmp_path / "c1", source_paths=str(chain_file), show_progress=False)
    assert r1["installed"] is True and r1["acquired"] == "downloaded"
    r2 = fetch_data(
        "A_to_B",
        cache_root=tmp_path / "c2",
        source_paths={"toy.chain": str(chain_file)},
        show_progress=False,
    )
    assert r2["installed"] is True and r2["acquired"] == "downloaded"


# --- #157: Content-addressed SHA-256 SourceStore ---


def test_source_store_reuses_by_sha256_re_verifies_on_hit_and_prunes(tmp_path, capsys):
    payload = b"verified-source-bytes\n" * 20
    digest = hashlib.sha256(payload).hexdigest()
    origin = tmp_path / "origin.tsv"
    origin.write_bytes(payload)

    cache_root = tmp_path / "cache"
    store = SourceStore(cache_root, seed_from_pyensembl=False)
    first = store.acquire(origin.as_uri(), digest, len(payload), download_options={"show_progress": False})
    assert first.is_file()
    assert str(source_store_dir(cache_root) / digest) in str(first)

    # Reuse by SHA-256 even when original file is gone and URL differs.
    origin.unlink()
    second = store.acquire("https://example.invalid/different-name.tsv", digest, len(payload))
    assert second == first and second.read_bytes() == payload

    # Tampering with stored bytes invalidates lookup() — never trusted by filename.
    first.chmod(0o644)
    first.write_bytes(b"X" * len(payload))
    assert store.lookup(digest, size=len(payload), filename="origin.tsv") is None

    # Re-acquiring repairs the tampered entry cleanly.
    origin.write_bytes(payload)
    repaired = store.acquire(origin.as_uri(), digest, len(payload), download_options={"show_progress": False})
    assert repaired.read_bytes() == payload

    # Prune dry-run reports cached sources without deleting; --yes deletes them.
    dry = prune_source_store(cache_root, apply=False)
    assert dry["dry_run"] is True and dry["removed_count"] == 1 and dry["reclaimable_bytes"] >= len(payload)
    assert store.lookup(digest, size=len(payload)) is not None

    run(["data", "prune", "--cache-root", str(cache_root), "--sources", "--json"])
    cli_dry = json.loads(capsys.readouterr().out)
    assert cli_dry["sources"]["removed_count"] == 1 and cli_dry["dry_run"] is True

    run(["data", "prune", "--cache-root", str(cache_root), "--sources", "--yes", "--json"])
    cli_yes = json.loads(capsys.readouterr().out)
    assert cli_yes["sources"]["removed_count"] == 1 and cli_yes["dry_run"] is False
    assert store.lookup(digest, size=len(payload)) is None


def test_source_store_seeds_from_pyensembl_cache_only_when_sha256_matches(tmp_path, monkeypatch):
    payload = b"ensembl-gtf-content\n" * 50
    digest = hashlib.sha256(payload).hexdigest()
    pyensembl_dir = tmp_path / "pyensembl-cache"
    pyensembl_dir.mkdir()
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(pyensembl_dir))

    # Wrong bytes with the exact same filename must NEVER be trusted.
    wrong = pyensembl_dir / "Homo_sapiens.GRCh38.101.gtf.gz"
    wrong.write_bytes(b"!" * len(payload))
    store = SourceStore(tmp_path / "orfdata-cache", seed_from_pyensembl=True)
    assert store.find_pyensembl_seed(digest, size=len(payload), filename=wrong.name) is None

    # Correct bytes (even under a different subdirectory/name) seed the store offline.
    nested = pyensembl_dir / "GRCh38" / "ensembl101"
    nested.mkdir(parents=True)
    good = nested / "matching.gtf.gz"
    good.write_bytes(payload)
    seeded = store.acquire(
        "https://example.invalid/unreachable/Homo_sapiens.GRCh38.101.gtf.gz",
        digest,
        len(payload),
        download_options={"show_progress": False},
    )
    assert seeded.read_bytes() == payload
    # ORFdata owns its own verified copy even if PyEnsembl's file is later removed.
    good.unlink()
    assert store.lookup(digest, size=len(payload)).read_bytes() == payload


def test_catalog_build_reuses_source_store_across_rebuilds(tmp_path):
    from test_catalog_registry import definition as catalog_definition
    from test_catalog_registry import setup_release

    manifest, source, original = setup_release(tmp_path)
    defs = catalog_definition(manifest)
    defs[manifest.catalog_id]["versions"]["1"]["assets"] = {}
    defs[manifest.catalog_id]["versions"]["1"]["manifest"]["sources"][0]["url"] = source.as_uri()
    registry = CatalogRegistry(defs, cache_root=tmp_path / "cache")

    built1 = registry.build(manifest.catalog_id, download_options={"show_progress": False})
    assert built1.records == original.records
    # Remove the original upstream file: forced rebuild reuses SourceStore offline.
    source.unlink()
    built2 = registry.build(manifest.catalog_id, force=True, download_options={"show_progress": False})
    assert built2.records == original.records


# --- #156: Resumable compressed single-file (decompress: true) datasets ---


def test_all_packaged_decompress_datasets_have_compressed_pins():
    data = json.loads(files("orfdata").joinpath("data/datasets.json").read_text())
    decompress_releases = []
    for name, entry in data.items():
        for version, release in entry["versions"].items():
            for asset_spec in release.get("assets", {}).values():
                if asset_spec.get("decompress"):
                    dist = release["distribution"]
                    assert isinstance(dist.get("download_sha256"), str) and len(dist["download_sha256"]) == 64
                    assert type(dist.get("download_size")) is int and dist["download_size"] > 0
                    decompress_releases.append((name, version))
    assert len(decompress_releases) == 11


def test_resumable_compressed_dataset_downloads_raw_gz_with_resume_verifies_and_decompresses(
    tmp_path, experiment, monkeypatch
):
    from test_measurements import definition as dataset_definition
    from test_measurements import write_dataset

    sqlite_ds = write_dataset(tmp_path / "orig.sqlite", experiment)
    raw_bytes = sqlite_ds.path.read_bytes()
    gz_bytes = gzip.compress(raw_bytes, mtime=0)
    gz_file = tmp_path / "orig.sqlite.gz"
    gz_file.write_bytes(gz_bytes)

    defs = dataset_definition(sqlite_ds)
    rel = defs["fixture"]["versions"]["1"]
    rel["distribution"]["download_sha256"] = hashlib.sha256(gz_bytes).hexdigest()
    rel["distribution"]["download_size"] = len(gz_bytes)
    rel["assets"]["dataset.sqlite"] = {
        "url": gz_file.as_uri(),
        "sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "size": len(raw_bytes),
        "decompress": True,
    }

    seen_calls = []
    real_install_bundle = datacache.bundles.install_bundle

    def spy_install_bundle(store_dir, assets, **kwargs):
        seen_calls.append((Path(store_dir), dict(assets), dict(kwargs.get("download_options") or {})))
        return real_install_bundle(store_dir, assets, **kwargs)

    monkeypatch.setattr(datacache.bundles, "install_bundle", spy_install_bundle)

    registry = DatasetRegistry(defs, cache_root=tmp_path / "cache")
    registry.download("fixture", "1", show_progress=False, resume="auto")

    # First call acquired the raw .gz with resume="auto" and compressed SHA-256;
    # second call decompressed the verified local .gz atomically into the dataset store.
    assert len(seen_calls) == 2
    first_dir, first_assets, first_opts = seen_calls[0]
    assert "sources" in first_dir.parts
    assert first_opts.get("resume") == "auto"
    only_asset = next(iter(first_assets.values()))
    assert only_asset["sha256"] == hashlib.sha256(gz_bytes).hexdigest()
    assert not only_asset.get("decompress")

    second_dir, second_assets, second_opts = seen_calls[1]
    assert "datasets" in second_dir.parts
    assert "resume" not in second_opts
    assert second_assets["dataset.sqlite"]["sha256"] == hashlib.sha256(raw_bytes).hexdigest()

    # Installed dataset opens and passes full integrity checks offline.
    opened = registry.open("fixture", "1", verify_files=True)
    assert list(opened.measurements()) == list(write_dataset(tmp_path / "check.sqlite", experiment).measurements())

    # Second download with resume="auto" is a no-op cache hit (no extra install_bundle calls).
    seen_calls.clear()
    registry.download("fixture", "1", show_progress=False, resume="auto")
    assert seen_calls == []


def test_resumable_compressed_dataset_rejects_corrupt_gz_before_decompression(tmp_path, experiment):
    from test_measurements import definition as dataset_definition
    from test_measurements import write_dataset

    sqlite_ds = write_dataset(tmp_path / "orig.sqlite", experiment)
    raw_bytes = sqlite_ds.path.read_bytes()
    gz_bytes = gzip.compress(raw_bytes, mtime=0)
    bad_gz = tmp_path / "bad.sqlite.gz"
    bad_gz.write_bytes(b"corrupted" + gz_bytes[9:])

    defs = dataset_definition(sqlite_ds)
    rel = defs["fixture"]["versions"]["1"]
    rel["distribution"]["download_sha256"] = hashlib.sha256(gz_bytes).hexdigest()
    rel["distribution"]["download_size"] = len(gz_bytes)
    rel["assets"]["dataset.sqlite"] = {
        "url": bad_gz.as_uri(),
        "sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "size": len(raw_bytes),
        "decompress": True,
    }
    registry = DatasetRegistry(defs, cache_root=tmp_path / "cache")
    with pytest.raises(ValueError):
        registry.download("fixture", "1", show_progress=False, resume="auto")
    assert registry.inspect("fixture", "1").status == "missing"
