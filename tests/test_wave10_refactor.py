"""Tests for Wave 10: viral DataCache split (#159), discovery sizes/species/defaults (#160),
unified status/action vocabulary (#164), and CLI/exit-status cleanups (#173)."""

import hashlib
import json
from datetime import datetime, timezone
from email.message import Message
from urllib.error import HTTPError

import pytest

from orfdata import (
    ACQUISITION_ACTIONS,
    RESOURCE_STATES,
    DatasetManifest,
    DatasetRegistry,
    ReferenceContext,
    SourceAsset,
    TermsNotAccepted,
    acquisition_action,
    catalog_status,
    data_info,
    fetch_data,
    list_data,
    normalize_status,
)
from orfdata.cli import run
from orfdata.cli_parser import ArgumentParser, share_cache_root, share_option
from orfdata.liftover import download_chain, pinned_chains
from orfdata.net import USER_AGENT, fetch, parse_retry_after
from orfdata.record_sources import (
    AUTHENTICATED_REQUEST_INTERVAL,
    REQUEST_INTERVAL,
    GenBankSource,
    ncbi_api_key,
)
from orfdata.registry import load_registry
from orfdata.viral_panel import (
    bundle_version,
    download_viral_references,
    open_viral_references,
    viral_store_path,
)
from orfdata.viral_registry import GDC_VIRUSES, ViralCollection


def _exit(argv, capsys):
    with pytest.raises(SystemExit) as exc:
        run(argv)
    return exc.value.code, capsys.readouterr().err


# --- #159: Viral DataCache download/build split, Retry-After, NCBI_API_KEY, User-Agent ---


def test_net_user_agent_does_not_name_viral_module():
    assert "viral" not in USER_AGENT.lower()
    assert USER_AGENT.startswith("orfdata/")


def test_parse_retry_after_handles_seconds_dates_and_invalid_values():
    now = datetime(2026, 11, 24, 12, 0, 0, tzinfo=timezone.utc)
    assert parse_retry_after("7", now=now) == 7.0
    assert parse_retry_after("2.5", now=now) == 2.5
    assert parse_retry_after("0", now=now) == 0.0
    assert parse_retry_after("Tue, 24 Nov 2026 12:00:15 GMT", now=now) == 15.0
    assert parse_retry_after("Tue, 24 Nov 2026 11:59:50 GMT", now=now) == 0.0
    assert parse_retry_after(None, now=now) is None
    assert parse_retry_after("", now=now) is None
    assert parse_retry_after("-3", now=now) is None
    assert parse_retry_after("not-a-date", now=now) is None


def test_net_fetch_honours_retry_after_header(monkeypatch):
    from orfdata import net

    sleeps = []
    calls = []

    class _Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return b"ok"

    def fake_urlopen(request, timeout=None):
        calls.append(request.headers.get("User-agent"))
        if len(calls) == 1:
            headers = Message()
            headers["Retry-After"] = "4"
            raise HTTPError(request.full_url, 429, "Too Many Requests", headers, None)
        return _Response()

    monkeypatch.setattr(net, "urlopen", fake_urlopen)
    monkeypatch.setattr(net.time, "sleep", sleeps.append)
    assert fetch("https://example.invalid/item", timeout=5) == b"ok"
    assert sleeps == [4.0]
    assert calls == [USER_AGENT, USER_AGENT]


def test_genbank_source_supports_ncbi_api_key(monkeypatch):
    from orfdata import record_sources

    seen_params = []
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    source = GenBankSource()
    assert ncbi_api_key() is None
    assert source.request_interval == REQUEST_INTERVAL

    def fake_fetch(url, *, timeout, params=None, **kwargs):
        seen_params.append(dict(params or {}))
        return b""

    monkeypatch.setattr(record_sources, "fetch", fake_fetch)
    monkeypatch.setattr(record_sources, "parse_genbank", lambda payload: ())
    with pytest.raises(ValueError, match="NCBI returned 0 records"):
        source.records(["ACC1.1"], timeout=5)
    assert "api_key" not in seen_params[-1]

    monkeypatch.setenv("NCBI_API_KEY", "secret-key-123")
    assert ncbi_api_key() == "secret-key-123"
    assert source.request_interval == AUTHENTICATED_REQUEST_INTERVAL
    with pytest.raises(ValueError, match="NCBI returned 0 records"):
        source.records(["ACC1.1"], timeout=5)
    assert seen_params[-1]["api_key"] == "secret-key-123"


def test_published_viral_collection_downloads_via_datacache_and_content_options_select_build(
    viral_panel, tmp_path, monkeypatch
):
    from dataclasses import replace

    from orfdata import viral_panel as vp
    from orfdata import viral_registry

    # First build a real fixture bundle so we have valid bundle files to publish as a release artifact.
    seed_cache = tmp_path / "seed-cache"
    seeded = download_viral_references(species="human", cache_root=str(seed_cache))
    artifact_dir = tmp_path / "published-files"
    artifact_dir.mkdir()
    assets = {}
    for file_path in sorted(seeded.path.rglob("*")):
        if not file_path.is_file() or file_path.name.startswith(".datacache-"):
            continue
        rel = file_path.relative_to(seeded.path).as_posix()
        dest = artifact_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        raw = file_path.read_bytes()
        dest.write_bytes(raw)
        assets[rel] = {
            "url": dest.as_uri(),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "size": len(raw),
        }
    total_size = sum(spec["size"] for spec in assets.values())
    published_collection = replace(
        GDC_VIRUSES,
        assets=assets,
        download_size=total_size,
        installed_size=total_size,
    )
    assert isinstance(published_collection, ViralCollection) and published_collection.published
    monkeypatch.setitem(viral_registry._COLLECTIONS, published_collection.name, published_collection)

    # Discovery reports download route and sizes.
    cache_root = tmp_path / "user-cache"
    info = data_info("gdc-viruses", cache_root=cache_root)
    assert info["acquire"] == "download"
    assert info["download_size"] == total_size
    assert info["installed_size"] == total_size

    # Default fetch uses DataCache VersionedDatasetRegistry download route (never calls _write_bundle).
    build_calls = []
    real_write = vp._write_bundle
    monkeypatch.setattr(vp, "_write_bundle", lambda *a, **k: build_calls.append(a) or real_write(*a, **k))

    fetched = fetch_data("gdc-viruses", cache_root=cache_root, show_progress=False)
    assert fetched["acquired"] == "downloaded"
    assert fetched["origin"] == "download"
    assert fetched["status"] == "available"
    assert build_calls == []
    dl_store = viral_store_path(cache_root, "homo_sapiens", published_collection, area="downloads")
    assert fetched["path"].startswith(str(dl_store))

    # Content-changing option (--inventory-molecules / prefer_reference_molecule=False) selects build route.
    built = download_viral_references(
        collection=published_collection,
        cache_root=cache_root,
        prefer_reference_molecule=False,
    )
    assert len(build_calls) == 1
    build_store = viral_store_path(cache_root, "homo_sapiens", published_collection, area="collections")
    assert str(built.path).startswith(str(build_store))
    assert open_viral_references(collection=published_collection, cache_root=cache_root).genome_count == 4

    # DataCache's own 'orfdata data prune' discovers and prunes stale viral generations automatically.
    download_viral_references(collection=published_collection, cache_root=cache_root, force=True)
    run(["data", "prune", "--cache-root", str(cache_root), "--yes", "--json"])


# --- #160: Species, sizes, default versions, --all-versions, and table layout ---


def test_data_list_shows_species_sizes_default_versions_and_clean_layout(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    default_rows = list_data()
    all_rows = list_data(all_versions=True)
    assert len(all_rows) > len(default_rows)
    assert all(row.default for row in default_rows)
    assert any(not row.default for row in all_rows)

    # Superseded versions are hidden by default and shown with --all-versions.
    run(["data", "list", "--kind", "catalog", "--json"])
    default_catalogs = json.loads(capsys.readouterr().out)
    assert [r["version"] for r in default_catalogs if r["name"] == "riborf-mouse"] == ["3"]
    bundled_row = next(r for r in default_catalogs if r["name"] == "bundled")
    assert "human GENCODE" in bundled_row["title"]
    assert bundled_row["status"] == "bundled" and bundled_row["installed"] is True

    run(["data", "list", "--kind", "catalog", "--all-versions", "--json"])
    all_catalogs = json.loads(capsys.readouterr().out)
    assert [r["version"] for r in all_catalogs if r["name"] == "riborf-mouse"] == ["1", "2", "3"]

    # Human table shows SPECIES, SIZE, shortened viral version, and never splits names mid-token.
    run(["data", "list", "--limit", "0"])
    table = capsys.readouterr().out
    header = table.splitlines()[0]
    assert all(col in header for col in ("KIND", "NAME", "SPECIES", "VERSION", "SIZE", "STATUS", "ACQUIRE"))
    assert "pfizer-atlas-dog-peptide-evidence:search.fasta" in table
    full_viral_version = bundle_version(GDC_VIRUSES)
    assert len(full_viral_version) > 60
    assert full_viral_version not in table
    assert "v7-" in table and "..." in table
    # Complete version remains in JSON.
    run(["data", "list", "--kind", "viral", "--json"])
    assert json.loads(capsys.readouterr().out)[0]["version"] == full_viral_version


def test_data_info_reports_download_and_installed_sizes_across_kinds(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    for name, expect_download, expect_installed in (
        ("pfizer-atlas-mouse-rna", True, True),
        ("riborf-mouse", True, True),
        ("openprot-mouse", True, False),
        ("bundled", False, True),
        ("t2t-chm13", True, False),
        ("GRCh38_to_GRCh37", True, True),
        ("gdc-viruses-records", True, True),
    ):
        info = data_info(name)
        assert (info["download_size"] is not None) is expect_download, name
        assert (info["installed_size"] is not None) is expect_installed, name
        run(["data", "info", name])
        human = capsys.readouterr().out
        if expect_download:
            assert "download size:" in human, name
        if expect_installed:
            assert "installed size:" in human, name


# --- #164: Unified status and action vocabulary ---


def test_unified_status_and_action_vocabulary():
    assert RESOURCE_STATES == ("available", "bundled", "missing", "damaged", "inaccessible")
    assert ACQUISITION_ACTIONS == (
        "downloaded",
        "built",
        "repaired",
        "replaced",
        "already installed",
        "skipped",
        "failed",
    )
    assert normalize_status("available") == "available"
    assert normalize_status("missing") == "missing"
    assert normalize_status("invalid") == "damaged"
    assert normalize_status("corrupt") == "damaged"
    assert normalize_status("damaged") == "damaged"
    assert normalize_status("inaccessible") == "inaccessible"
    assert normalize_status("available", bundled=True) == "bundled"
    with pytest.raises(ValueError, match="Unknown resource status"):
        normalize_status("broken")

    assert acquisition_action("available", acquire="download", force=False) == "already installed"
    assert acquisition_action("bundled", acquire=None, force=True) == "already installed"
    assert acquisition_action("available", acquire="download", force=True) == "replaced"
    assert acquisition_action("missing", acquire="download", force=False) == "downloaded"
    assert acquisition_action("missing", acquire="build", force=False) == "built"
    assert acquisition_action("invalid", acquire="download", force=False) == "repaired"
    assert acquisition_action("inaccessible", acquire="build", force=True) == "repaired"

    for row in catalog_status():
        assert row.status in RESOURCE_STATES


def test_setup_uses_unified_actions_for_chains(tmp_path, monkeypatch):
    from test_catalog_registry import definition as catalog_definition
    from test_catalog_registry import setup_release

    from orfdata import install, liftover
    from orfdata.catalogs import CatalogRegistry

    manifest, source, _ = setup_release(tmp_path)
    artifact = tmp_path / "catalog.json.gz"
    build_reg = CatalogRegistry(catalog_definition(manifest), cache_root=tmp_path / "src")
    build_reg.build(manifest.catalog_id, source_paths={"source.json": source},
                    download_options={"show_progress": False}).save(artifact)
    defs = catalog_definition(manifest, artifact)
    defs[manifest.catalog_id]["versions"][manifest.version]["assets"]["catalog.json.gz"]["url"] = artifact.as_uri()
    registry = CatalogRegistry(defs, cache_root=tmp_path / "cache")
    registry.download(manifest.catalog_id, show_progress=False)

    # Simulate an already-installed chain and verify setup reports "already installed" and "replaced".
    chain_state = {"status": "available"}

    class _FakeChainRegistry:
        def inspect(self, *args, **kwargs):
            return type("Insp", (), {"status": chain_state["status"]})()

    monkeypatch.setattr(liftover, "chain_registry", lambda chain, root: _FakeChainRegistry())
    monkeypatch.setattr(liftover, "download_chain", lambda *a, **k: tmp_path / "chain")
    result = install.install_species(
        "mouse", registry=registry, catalogs=[manifest.catalog_id], collections="none", accept_terms=True,
    )
    assert all(action == "already installed" for action in result["chains"].values())
    forced = install.install_species(
        "mouse", registry=registry, catalogs=[manifest.catalog_id], collections="none",
        accept_terms=True, force=True,
    )
    assert all(action == "replaced" for action in forced["chains"].values())


# --- #173: Deferred CLI and exit-status cleanups ---


def test_mismatched_reference_and_asset_exits_2(tmp_path, monkeypatch, capsys):
    from test_measurements import definition as dataset_definition
    from test_measurements import write_dataset

    from orfdata import SearchReferenceSpec, dataset_cli, datasets

    f1 = tmp_path / "ref1.fasta"
    f1.write_text(">P1\nMAAA\n")
    f2 = tmp_path / "ref2.fasta"
    f2.write_text(">P2\nMCCC\n")
    ref1 = ReferenceContext(9606, "human", None, "uniprot", "1",
                            sequence_sha256=hashlib.sha256(f1.read_bytes()).hexdigest())
    ref2 = ReferenceContext(10090, "mouse", None, "uniprot", "1",
                            sequence_sha256=hashlib.sha256(f2.read_bytes()).hexdigest())
    src1 = SourceAsset("human.fasta", f1.as_uri(), hashlib.sha256(f1.read_bytes()).hexdigest(),
                       f1.stat().st_size, "CC0", redistributable=True)
    src2 = SourceAsset("mouse.fasta", f2.as_uri(), hashlib.sha256(f2.read_bytes()).hexdigest(),
                       f2.stat().st_size, "CC0", redistributable=True)
    manifest = DatasetManifest(
        "multi-ref", "1", "Two search references", (ref1, ref2), (src1, src2),
        ("test",), ("digest_ms",), {"rows": 0}, {"software": "test", "version": "1", "parameters": {}},
    )
    sqlite_path = write_dataset(tmp_path / "multi.sqlite", (manifest, (), (), ()))
    registry = DatasetRegistry(dataset_definition(sqlite_path), cache_root=tmp_path / "cache")
    specs = (
        SearchReferenceSpec("multi-ref", "1", ref1, src1, "uniprot", "uniprot-fasta"),
        SearchReferenceSpec("multi-ref", "1", ref2, src2, "uniprot", "uniprot-fasta"),
    )
    monkeypatch.setattr(registry, "search_references", lambda name, version=None: specs)
    monkeypatch.setattr(dataset_cli, "DatasetRegistry", lambda definitions=None, **kw: registry)
    monkeypatch.setattr(datasets, "DatasetRegistry", lambda definitions=None, **kw: registry)

    for argv in (
        ["datasets", "inspect-reference", "multi-ref", "--reference", ref1.key, "--asset", "mouse.fasta"],
        ["datasets", "download-reference", "multi-ref", "--reference", ref1.key, "--asset", "mouse.fasta"],
        ["datasets", "query-reference", "multi-ref", "--target", "P1", "--reference", ref1.key, "--asset", "mouse.fasta"],
        ["data", "reference", "multi-ref", "--reference", ref1.key, "--asset", "mouse.fasta"],
    ):
        status, err = _exit(argv, capsys)
        assert status == 2, (argv, err)
        assert "do not match the same search reference" in err


def test_terms_not_accepted_is_value_error_and_exits_2(tmp_path, capsys):
    chain = next(c for c in pinned_chains() if c.restricted)
    with pytest.raises(TermsNotAccepted) as err:
        download_chain(chain.source_build, chain.target_build, species=chain.species, cache_root=tmp_path)
    assert isinstance(err.value, ValueError)
    with pytest.raises(TermsNotAccepted):
        fetch_data(chain.name, cache_root=tmp_path)
    for argv in (
        ["chains", "download", chain.source_build, chain.target_build],
        ["data", "chains", chain.source_build, chain.target_build],
        ["data", "fetch", chain.name],
    ):
        status, stderr = _exit(argv, capsys)
        assert status == 2 and "--accept-terms" in stderr, argv


def test_shared_options_propagate_to_subcommands_added_later_and_registry_works_after_subcommand(
    tmp_path, capsys
):
    parser = ArgumentParser(prog="tool")
    share_cache_root(parser)
    share_option(parser, "--registry", default=None, help="custom registry")
    commands = parser.add_subparsers(dest="command", required=True)
    late = commands.add_parser("late")
    nested = late.add_subparsers(dest="sub", required=True)
    nested.add_parser("inner")

    for argv in (
        ["--cache-root", "C", "--registry", "R", "late", "inner"],
        ["late", "--cache-root", "C", "inner", "--registry", "R"],
        ["late", "inner", "--cache-root", "C", "--registry", "R"],
    ):
        parsed = parser.parse_args(argv)
        assert (parsed.cache_root, parsed.registry) == ("C", "R"), argv

    # Real CLI accepts --registry after the subcommand too while keeping it hidden from help.
    bad_reg = tmp_path / "bad.json"
    bad_reg.write_text("{bad")
    for group in ("catalogs", "datasets"):
        status, err = _exit([group, "list", "--registry", str(bad_reg)], capsys)
        assert status == 1 and str(bad_reg) in err

    with pytest.raises(ValueError, match="Invalid dataset registry"):
        load_registry(DatasetRegistry, bad_reg)


def test_data_list_emits_single_truncation_footer(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    run(["data", "list", "--limit", "5"])
    output = capsys.readouterr().out
    assert output.count("Showing ") == 1
    assert "not shown:" in output and "--limit 0" in output
