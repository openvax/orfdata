"""Downloads restart unless resuming is requested; when it is, every acquisition path asks DataCache."""

import datacache
import datacache.bundles
import pytest

from orfdata import CatalogRegistry, DatasetRegistry
from orfdata.genomes import install_genome
from orfdata.install import _fetch
from orfdata.liftover import download_chain


class Stopped(Exception):
    """Raised at the DataCache boundary, before anything is downloaded or written."""


@pytest.fixture
def boundary(monkeypatch):
    calls = []

    def stop_at(kind):
        def stop(*args, download_options=None, **kwargs):
            calls.append((kind, dict(download_options or {})))
            raise Stopped(kind)
        return stop

    monkeypatch.setattr(datacache.bundles, "install_bundle", stop_at("bundle"))
    monkeypatch.setattr(datacache, "install_archive", stop_at("archive"))
    monkeypatch.setattr(datacache, "materialize", stop_at("materialize"))
    return calls


ACQUISITIONS = {
    "published catalogue": ("bundle", lambda root, options: CatalogRegistry(cache_root=root).download(
        "ncorf-mammals-mouse", **options)),
    "published single-file dataset": ("bundle", lambda root, options: DatasetRegistry(cache_root=root).download(
        "copy-yeast-abundance", **options)),
    "published archive parts": ("archive", lambda root, options: DatasetRegistry(cache_root=root).download(
        "pfizer-atlas-dog-peptide-evidence", **options)),
    "search reference": ("bundle", lambda root, options: DatasetRegistry(cache_root=root).download_reference(
        "pfizer-atlas-dog-protein", **options)),
    "pinned genome": ("bundle", lambda root, options: install_genome(
        "t2t-chm13", cache_root=root, download_options=options)),
    "pinned chain": ("bundle", lambda root, options: download_chain(
        "GRCh38", "GRCh37", cache_root=root, download_options=options)),
    "catalogue build sources": ("materialize", lambda root, options: CatalogRegistry(cache_root=root).build(
        "openprot-mouse", download_options=options)),
    "dataset build sources": ("materialize", lambda root, options: DatasetRegistry(cache_root=root).build(
        "peptideatlas-dog-peptides", download_options=options)),
    "archived build source": ("archive", lambda root, options: DatasetRegistry(cache_root=root).build(
        "pfizer-atlas-dog-peptide-evidence", download_options=options)),
    "setup's shared source": ("bundle", lambda root, options: _fetch(
        "https://example.invalid/shared.tsv", "0" * 64, 10, root, options)),
}


@pytest.mark.parametrize("path", ACQUISITIONS)
def test_downloads_restart_unless_resume_is_requested(path, boundary, tmp_path):
    kind, acquire = ACQUISITIONS[path]
    for options in (dict(timeout=5), dict(timeout=5, resume="auto")):
        with pytest.raises(Stopped):
            acquire(tmp_path, options)
    assert [found for found, _ in boundary] == [kind, kind]
    assert "resume" not in boundary[0][1] and boundary[1][1]["resume"] == "auto"


def test_fetch_and_setup_ask_for_resume_only_when_told(boundary, tmp_path):
    from orfdata import fetch_data, list_data
    from orfdata.install import InstallError, install_species

    names = ["copy-yeast-abundance", next(r.name for r in list_data(kind="reference")), "t2t-chm13",
             next(r.name for r in list_data(kind="chain") if not r.restricted)]
    for resume in (False, True):
        boundary.clear()
        for name in names:  # A dataset, a search reference, a genome and a chain.
            with pytest.raises(Stopped):
                fetch_data(name, cache_root=tmp_path, show_progress=False, resume=resume)
        with pytest.raises(InstallError):  # Setup reports each failure after trying everything.
            install_species("yeast", cache_root=tmp_path, catalogs="all", collections="none", resume=resume)
        assert len(boundary) > len(names)
        assert all(("resume" in options) is resume for _, options in boundary)


def test_resume_needs_a_boolean_and_is_dropped_for_local_sources(boundary, tmp_path):
    from orfdata import fetch_data
    from orfdata.install import install_species

    with pytest.raises(ValueError, match="resume must be True or False"):
        fetch_data("copy-yeast-abundance", cache_root=tmp_path, resume="false")
    with pytest.raises(ValueError, match="resume must be True or False"):
        install_species("yeast", cache_root=tmp_path, resume=1)
    with pytest.raises(Stopped):
        fetch_data("copy-yeast-abundance", cache_root=tmp_path, show_progress=False, resume=True,
                   source_paths={"dataset.sqlite": str(tmp_path / "local.sqlite")})
    assert "resume" not in boundary[0][1]


def test_cli_passes_resume_through(monkeypatch, capsys):
    from orfdata import data_resources, install
    from orfdata.cli import run

    calls = []
    monkeypatch.setattr(data_resources, "check_fetch", lambda *a, **k: None)
    monkeypatch.setattr(data_resources, "fetch_data", lambda name, **options: calls.append(options) or {})
    monkeypatch.setattr(install, "install_species", lambda species, **options: calls.append(options) or {})
    run(["data", "fetch", "copy-yeast-abundance", "--resume", "--json"])
    run(["data", "fetch", "copy-yeast-abundance", "--json"])
    run(["setup", "--species", "yeast", "--resume", "--json"])
    run(["setup", "--species", "yeast", "--json"])
    capsys.readouterr()
    assert [call["resume"] for call in calls] == [True, False, True, False]
