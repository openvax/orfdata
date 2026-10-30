"""Every acquisition path asks DataCache to resume an interrupted download, unless told not to."""

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
    "published catalogue": ("bundle", lambda root: CatalogRegistry(cache_root=root).download("ncorf-mammals-mouse")),
    "published single-file dataset": (
        "bundle", lambda root: DatasetRegistry(cache_root=root).download("copy-yeast-abundance")),
    "published archive parts": (
        "archive", lambda root: DatasetRegistry(cache_root=root).download("pfizer-atlas-dog-peptide-evidence")),
    "search reference": (
        "bundle", lambda root: DatasetRegistry(cache_root=root).download_reference("pfizer-atlas-dog-protein")),
    "pinned genome": ("bundle", lambda root: install_genome("t2t-chm13", cache_root=root)),
    "pinned chain": ("bundle", lambda root: download_chain("GRCh38", "GRCh37", cache_root=root)),
    "catalogue build sources": ("materialize", lambda root: CatalogRegistry(cache_root=root).build("openprot-mouse")),
    "dataset build sources": (
        "materialize", lambda root: DatasetRegistry(cache_root=root).build("peptideatlas-dog-peptides")),
    "archived build source": (
        "archive", lambda root: DatasetRegistry(cache_root=root).build("pfizer-atlas-dog-peptide-evidence")),
    "setup's shared source": ("bundle", lambda root: _fetch(
        "https://example.invalid/shared.tsv", "0" * 64, 10, root, dict(timeout=5))),
}


@pytest.mark.parametrize("path", ACQUISITIONS)
def test_acquisition_resumes_by_default(path, boundary, tmp_path):
    kind, acquire = ACQUISITIONS[path]
    with pytest.raises(Stopped):
        acquire(tmp_path)
    assert boundary == [(kind, boundary[0][1])]
    assert boundary[0][1]["resume"] == "auto"


def test_explicit_settings_win(boundary, tmp_path):
    registry = DatasetRegistry(cache_root=tmp_path)
    with pytest.raises(Stopped):
        registry.download("copy-yeast-abundance", resume=False, timeout=7)
    with pytest.raises(Stopped):
        registry.build("peptideatlas-dog-peptides", download_options=dict(resume=True))
    with pytest.raises(Stopped):
        _fetch("https://example.invalid/shared.tsv", "0" * 64, 10, tmp_path, dict(timeout=5))
    assert [options for _, options in boundary] == [
        dict(resume=False, timeout=7), dict(resume=True), dict(resume="auto", timeout=5)]
