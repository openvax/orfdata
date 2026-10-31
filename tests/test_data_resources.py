"""One acquisition contract for downloaded, source-built and bundled resources."""

import json
import os
from pathlib import Path

import pytest
from test_catalog_registry import definition as catalog_definition
from test_catalog_registry import setup_release
from test_measurements import definition as dataset_definition
from test_measurements import write_dataset

from orfdata import (
    CatalogRegistry,
    DatasetRegistry,
    NotInstalledError,
    data_info,
    data_path,
    fetch_data,
)
from orfdata.cli import run


@pytest.fixture
def published(tmp_path, experiment, monkeypatch):
    from orfdata import datasets

    original = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions = dataset_definition(original)
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    monkeypatch.setattr(datasets, "DatasetRegistry", lambda **kw: registry)
    return registry, original


def test_download_reuse_corruption_and_force_preserve_existing_generation(published):
    registry, original = published
    assert data_info("fixture")["acquire"] == "download"
    with pytest.raises(NotInstalledError, match="not installed; install it with: orfdata data fetch fixture"):
        data_path("fixture")
    installed = fetch_data("fixture", show_progress=False)
    assert installed["acquired"] == "downloaded"
    assert data_path("fixture").exists()
    assert list(registry.open("fixture").measurements()) == list(original.measurements())
    assert fetch_data("fixture")["acquired"] == "already installed"
    artifact = registry.open("fixture").path
    _flip_last_byte(artifact)
    with pytest.raises(ValueError, match="invalid"):
        data_path("fixture", verify_files=True)
    fetch_data("fixture", force=True, show_progress=False)
    assert artifact.exists() and registry.open("fixture").path != artifact
    assert list(registry.open("fixture").measurements()) == list(original.measurements())


@pytest.mark.parametrize("species", ["homo_sapiens", "mus_musculus", "canis_lupus_familiaris", "saccharomyces_cerevisiae"])
def test_pinned_source_build_uses_same_contract_for_each_species(tmp_path, monkeypatch, species):
    from orfdata import catalogs

    manifest, source, original = setup_release(tmp_path, species)
    definitions = catalog_definition(manifest)
    definitions[manifest.catalog_id]["versions"]["1"]["assets"] = {}
    registry = CatalogRegistry(definitions, cache_root=tmp_path / "cache")
    monkeypatch.setattr(catalogs, "CatalogRegistry", lambda **kw: registry)
    before = data_info(manifest.catalog_id)
    assert before["acquire"] == "build" and not before["distribution"]["redistributable"]
    result = fetch_data(manifest.catalog_id, source_paths={"source.json": source}, show_progress=False)
    assert result["acquired"] == "built" and result["installed"]
    source.unlink()
    assert fetch_data(manifest.catalog_id)["acquired"] == "already installed"
    assert registry.open_build(manifest.catalog_id).records == original.records
    assert data_path(manifest.catalog_id).exists()


@pytest.mark.parametrize("command", ["data", "downloads", "download"])
def test_cli_aliases_fetch_info_and_path(published, capsys, command):
    run([command, "info", "fixture", "--json"])
    assert json.loads(capsys.readouterr().out)["acquire"] == "download"
    run([command, "fetch", "fixture", "fixture", "--json"])
    output = json.loads(capsys.readouterr().out)
    assert len(output) == 1 and output[0]["installed"]
    run([command, "path", "fixture"])
    assert Path(capsys.readouterr().out.strip()).is_dir()


def test_batch_resolves_all_names_before_acquiring(published, capsys):
    registry, _ = published
    with pytest.raises(SystemExit):
        run(["data", "fetch", "fixture", "does-not-exist"])
    assert "Unknown data resource" in capsys.readouterr().err
    assert registry.inspect("fixture").status == "missing"
    with pytest.raises(SystemExit):
        run(["data", "fetch", "fixture", "bundled", "--version", "1"])
    assert "exactly one resource" in capsys.readouterr().err


def test_default_and_historical_catalogue_versions_are_explicit():
    assert data_info("riborf-mouse")["version"] == "3"
    old = data_info("riborf-mouse", version="1")
    assert old["version"] == "1" and old["manifest"]["coverage"]["records"] == 40089
    with pytest.raises(ValueError, match="Unknown catalogue version"):
        data_info("riborf-mouse", version="nonexistent")


def test_bundled_path_and_fetch_need_no_download(monkeypatch):
    import socket

    monkeypatch.setattr(socket, "create_connection", lambda *a, **k: pytest.fail("network"))
    path = data_path("bundled")
    assert path.name == "catalog.json.gz" and path.is_file()
    assert fetch_data("bundled")["acquired"] == "already installed"
    with pytest.raises(ValueError, match="has version"):
        data_info("bundled", version="unknown")


def test_restricted_terms_preflight_precedes_other_downloads(published, capsys):
    from orfdata import list_data

    chain = next(r for r in list_data(kind="chain") if r.restricted)
    with pytest.raises(SystemExit):
        run(["data", "fetch", "fixture", chain.name])
    assert "requires --accept-terms" in capsys.readouterr().err
    assert published[0].inspect("fixture").status == "missing"
    assert data_info(chain.name)["source"]["terms"]


def test_reference_and_genome_metadata_and_invalid_override():
    from orfdata import list_data

    reference = next(r for r in list_data(kind="reference"))
    info = data_info(reference.name)
    assert info["sources"][0]["sha256"]
    with pytest.raises(ValueError, match="Source overrides"):
        fetch_data(reference.name, source_paths={"invalid": "missing"})
    genome = next(r for r in list_data(kind="genome"))
    assert data_info(genome.name)["source"]["assets"]
    with pytest.raises(ValueError, match="Unknown data resource"):
        data_info(reference.name, kind="catalog")


def test_every_registered_resource_resolves_without_network(monkeypatch):
    import socket

    from orfdata import list_data

    monkeypatch.setattr(socket, "getaddrinfo", lambda *a, **k: pytest.fail("network"))
    for row in list_data():
        info = data_info(row.name, kind=row.kind, version=row.version)
        assert info["name"] == row.name and info["version"] == row.version
        if row.acquire is not None:
            assert row.command.startswith("orfdata data fetch ")


def test_fetch_json_separates_source_progress(monkeypatch, published, capsys):
    from orfdata import data_resources

    original = data_resources.fetch_data

    def verbose(*args, **kwargs):
        print("source processing progress")
        return original(*args, **kwargs)

    monkeypatch.setattr(data_resources, "fetch_data", verbose)
    run(["data", "fetch", "fixture", "--json"])
    output = capsys.readouterr()
    assert json.loads(output.out)[0]["installed"]
    assert "source processing progress" in output.err


def _flip_last_byte(path):
    with path.open("r+b") as handle:
        handle.seek(-1, 2)
        last = handle.read(1)
        handle.seek(-1, 2)
        handle.write(bytes([last[0] ^ 1]))


def test_status_shows_and_fetch_repairs_damaged_files(published, capsys):
    registry, original = published
    fetch_data("fixture", show_progress=False)
    artifact = registry.open("fixture").path
    _flip_last_byte(artifact)
    run(["data", "status", "--kind", "dataset", "--verify"])
    output = capsys.readouterr().out
    assert "invalid" in output and "Damaged; 'data fetch' replaces them:" in output
    assert "  orfdata data fetch fixture --kind dataset --version " in output
    run(["data", "status", "--kind", "dataset", "--verify", "--json"])
    assert [row["status"] for row in json.loads(capsys.readouterr().out)] == ["invalid"]
    repaired = fetch_data("fixture", show_progress=False)
    assert repaired["acquired"] == "repaired" and repaired["installed"]
    assert artifact.exists() and registry.open("fixture").path != artifact
    assert list(registry.open("fixture").measurements()) == list(original.measurements())
    assert fetch_data("fixture")["acquired"] == "already installed"


def test_status_omits_absent_resources(published, capsys):
    run(["data", "status", "--kind", "dataset", "--json"])
    assert json.loads(capsys.readouterr().out) == []
    run(["data", "list", "--kind", "dataset", "--json"])
    assert [row["status"] for row in json.loads(capsys.readouterr().out)] == ["missing"]


unreadable_files_are_readable = pytest.mark.skipif(
    os.name != "posix" or os.geteuid() == 0, reason="root reads files without read permission")


@unreadable_files_are_readable
def test_fetch_does_not_replace_unreadable_files_without_force(published, monkeypatch, capsys):
    from orfdata import data_resources

    registry, _ = published
    fetch_data("fixture", show_progress=False)
    artifact = registry.open("fixture").path
    artifact.chmod(0)
    try:
        assert data_info("fixture", verify_files=True)["status"] == "inaccessible"
        with pytest.raises(ValueError, match="can't be read.*command line: orfdata data fetch fixture --kind dataset .* --force$"):
            fetch_data("fixture", show_progress=False)
        run(["data", "status", "--kind", "dataset", "--verify"])
        assert "Unreadable; fix their permissions, or replace them:" in capsys.readouterr().out
        # A multi-name fetch refuses before it downloads anything.
        monkeypatch.setattr(data_resources, "fetch_data", lambda *a, **k: pytest.fail("downloaded"))
        with pytest.raises(SystemExit):
            run(["data", "fetch", "bundled", "fixture"])
        assert "can't be read" in capsys.readouterr().err
    finally:
        artifact.chmod(0o644)


def test_fetch_repairs_a_damaged_or_incomplete_viral_collection(viral_panel, tmp_path):
    from orfdata.downloads import collection_state

    root = tmp_path / "cache"
    assert fetch_data("gdc-viruses", cache_root=root)["acquired"] == "downloaded"
    path = data_path("gdc-viruses", cache_root=root)
    (path / "README.md").write_text("changed\n")
    assert data_info("gdc-viruses", cache_root=root, verify_files=True)["status"] == "corrupt"
    assert fetch_data("gdc-viruses", cache_root=root)["acquired"] == "repaired"
    assert data_info("gdc-viruses", cache_root=root, verify_files=True)["status"] == "available"
    (data_path("gdc-viruses", cache_root=root) / "metadata.json").unlink()
    assert collection_state("gdc-viruses", cache_root=root) == "corrupt"
    assert fetch_data("gdc-viruses", cache_root=root)["acquired"] == "repaired"
    assert collection_state("gdc-viruses", cache_root=root) == "available"


@unreadable_files_are_readable
def test_an_unreadable_viral_collection_is_not_replaced_without_force(viral_panel, tmp_path):
    from orfdata.downloads import collection_state

    root = tmp_path / "cache"
    fetch_data("gdc-viruses", cache_root=root)
    metadata = data_path("gdc-viruses", cache_root=root) / "metadata.json"
    metadata.chmod(0)
    try:
        assert collection_state("gdc-viruses", cache_root=root) == "inaccessible"
        with pytest.raises(ValueError, match="can't be read"):
            fetch_data("gdc-viruses", cache_root=root)
        assert fetch_data("gdc-viruses", cache_root=root, force=True)["acquired"] == "repaired"
    finally:
        if metadata.exists():
            metadata.chmod(0o644)
    assert collection_state("gdc-viruses", cache_root=root) == "available"


def test_listed_fetch_commands_name_the_cache_root(tmp_path):
    from orfdata import list_data

    row = next(r for r in list_data(kind="chain", cache_root=tmp_path))
    assert f" --cache-root {tmp_path}" in row.command
    assert list_data(kind="chain", cache_root=tmp_path, local=True) == ()


@unreadable_files_are_readable
def test_an_unreadable_records_directory_is_inaccessible_not_corrupt(viral_panel, tmp_path):
    from orfdata.downloads import collection_state

    root = tmp_path / "cache"
    fetch_data("gdc-viruses", cache_root=root)
    path = data_path("gdc-viruses", cache_root=root)
    records = path / "records"
    records.chmod(0)
    try:
        for verify in (False, True):
            assert collection_state("gdc-viruses", cache_root=root, verify_files=verify) == "inaccessible"
        with pytest.raises(ValueError, match="can't be read"):
            fetch_data("gdc-viruses", cache_root=root)
        with pytest.raises(ValueError, match="can't be read.*on the command line: orfdata data fetch gdc-viruses"):
            data_path("gdc-viruses", cache_root=root, verify_files=True)
    finally:
        records.chmod(0o755)
    assert sorted(p.name for p in path.iterdir())  # Nothing was deleted.
    assert collection_state("gdc-viruses", cache_root=root, verify_files=True) == "available"


@unreadable_files_are_readable
def test_batch_fetch_checks_verified_state_of_every_name_first(published, viral_panel, tmp_path, capsys):
    registry, _ = published
    root = tmp_path / "viral-cache"
    fetch_data("gdc-viruses", cache_root=root)
    fasta = data_path("gdc-viruses", cache_root=root) / "genomes.fasta"
    fasta.chmod(0)
    try:
        with pytest.raises(SystemExit):
            run(["data", "fetch", "fixture", "gdc-viruses", "--cache-root", str(root)])
        assert "can't be read" in capsys.readouterr().err
        assert registry.inspect("fixture").status == "missing"
    finally:
        fasta.chmod(0o644)


def test_setup_judges_whichever_catalogue_store_is_damaged():
    from types import SimpleNamespace

    from orfdata.install import _damage

    def registry(published, download, build):
        return SimpleNamespace(published=lambda *a: published,
                               inspect=lambda *a, **k: SimpleNamespace(status=download),
                               inspect_build=lambda *a, **k: SimpleNamespace(status=build))

    assert _damage(registry(True, "missing", "inaccessible"), "c", "1") == "inaccessible"
    assert _damage(registry(True, "invalid", "missing"), "c", "1") == "invalid"
    assert _damage(registry(False, "missing", "invalid"), "c", "1") == "invalid"
    assert _damage(registry(True, "available", "missing"), "c", "1") == "missing"


def test_default_list_preview_shows_every_kind(monkeypatch, tmp_path, capsys):
    from orfdata import list_data

    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    rows = list_data()
    kinds = {row.kind for row in rows}
    assert len(rows) > 25 and len(kinds) > 1
    run(["data", "list"])
    output = capsys.readouterr().out
    table = output.split("\n\n")[0].splitlines()[1:]
    assert {line.split()[0] for line in table if line[:1].strip()} == kinds
    assert f"Showing 25 of {len(rows)} entries; not shown:" in output
    run(["data", "list", "--limit", "0"])
    assert "not shown" not in capsys.readouterr().out


def test_preview_shares_rows_among_kinds():
    from types import SimpleNamespace

    from orfdata.download_cli import _preview

    rows = [SimpleNamespace(kind=kind, name=f"{kind}{index}")
            for kind, count in (("a", 6), ("b", 1), ("c", 4)) for index in range(count)]
    shown, hidden = _preview(rows, 5)
    assert [row.name for row in shown] == ["a0", "a1", "b0", "c0", "c1"]
    assert hidden == {"a": 4, "c": 2}
    assert _preview(rows, 0) == (rows, {}) and _preview(rows, 11) == (rows, {})


def test_fetch_summary_names_the_resource_and_what_happened(published, monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path / "isolated"))
    run(["data", "fetch", "fixture"])
    output = capsys.readouterr().out
    assert "name: fixture" in output and "version: " in output and "acquired: downloaded" in output
    run(["data", "fetch", "fixture"])
    assert "acquired: already installed" in capsys.readouterr().out
    # The focused command reports reuse the same way, rather than claiming another download.
    run(["data", "dataset", "fixture"])
    output = capsys.readouterr().out
    assert "acquired: already installed" in output and "title: " in output
    run(["data", "fetch", "fixture", "bundled"])
    assert "\n\nkind: catalog" in capsys.readouterr().out


def test_info_on_a_viral_collection_names_it(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    run(["data", "info", "gdc-viruses"])
    output = capsys.readouterr().out
    assert "name: gdc-viruses" in output and "version: v" in output and "command: orfdata data fetch" in output


def test_preview_rejects_a_negative_limit():
    from orfdata.download_cli import _preview

    with pytest.raises(ValueError, match="nonnegative"):
        _preview([object()], -1)


def test_viral_fetch_summary_names_the_collection(viral_panel, tmp_path, capsys):
    run(["data", "fetch", "gdc-viruses", "--cache-root", str(tmp_path)])
    output = capsys.readouterr().out
    assert "kind: viral" in output and "name: gdc-viruses" in output
    assert "version: v" in output and "acquired: downloaded" in output
