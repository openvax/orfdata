"""ORFdata's cache directory, store names and single-release DataCache registries."""

import hashlib
import json
from pathlib import Path

import datacache
import pytest

from orfdata import NotInstalledError
from orfdata.cache import (
    DOWNLOAD_TIMEOUT,
    content_digest,
    fetch_pinned_file,
    get_data_dir,
    installed_path,
    pinned_registry,
)
from orfdata.genomes import genome_registry, pinned_genomes
from orfdata.liftover import chain_registry, pinned_chains


@pytest.fixture
def pinned_file(tmp_path):
    """A local file and its DataCache pin, so nothing touches the network."""
    source = tmp_path / "source" / "table.tsv"
    source.parent.mkdir()
    source.write_bytes(b"gene\tvalue\nTP53\t1\n")
    return source, dict(url=source.as_uri(), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        size=source.stat().st_size)


def test_data_dir_is_the_orfdata_child_of_the_chosen_root_and_is_not_created(monkeypatch, tmp_path):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path / "orfdata-cache"))
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(tmp_path / "openvax-cache"))
    assert get_data_dir() == tmp_path / "orfdata-cache" / "orfdata"
    assert get_data_dir(tmp_path / "explicit") == tmp_path / "explicit" / "orfdata"
    monkeypatch.delenv("ORFDATA_CACHE")
    assert get_data_dir() == tmp_path / "openvax-cache" / "orfdata"
    assert not any(tmp_path.iterdir())


def test_content_digest_is_the_sha256_of_canonical_json():
    value = {"b": [1, 2], "a": {"y": None, "x": "é"}}
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    assert content_digest(value) == hashlib.sha256(canonical).hexdigest()
    assert content_digest(dict(reversed(value.items()))) == content_digest(value)
    assert content_digest({"a": 1}) != content_digest({"a": 2})
    with pytest.raises(TypeError):
        content_digest({"a": object()})


def test_pinned_registry_installs_one_release_in_its_store(tmp_path, pinned_file):
    source, pin = pinned_file
    store = tmp_path / "store"
    registry = pinned_registry("table", "7", {"table.tsv": pin}, store)
    assert registry.inspect("table", "7", verify_files=False).status == "missing"
    assert not store.exists()  # Building or inspecting a registry writes nothing.
    registry.download("table", show_progress=False)  # Its one version is the default.
    path = registry.local_path("table", "7", asset="table.tsv")
    assert path.read_bytes() == source.read_bytes()
    assert store in path.parents


def test_fetch_pinned_file_keeps_each_pin_in_its_own_store(tmp_path, pinned_file):
    source, pin = pinned_file
    path = fetch_pinned_file(pin["url"], pin["sha256"], pin["size"], tmp_path / "staging",
                             download_options=dict(show_progress=False))
    assert isinstance(path, str) and Path(path).read_bytes() == source.read_bytes()
    assert (tmp_path / "staging" / pin["sha256"]) in Path(path).parents
    with pytest.raises(datacache.FileValidationError):
        fetch_pinned_file(pin["url"], "0" * 64, pin["size"], tmp_path / "staging")


def test_installed_path_names_the_install_command_only_when_nothing_is_installed(tmp_path, pinned_file):
    _, pin = pinned_file
    registry = pinned_registry("table", "1", {"table.tsv": pin}, tmp_path / "store")
    asked = []

    def hint():
        asked.append(True)
        return "install it with: orfdata data fetch table"

    with pytest.raises(NotInstalledError, match="orfdata data fetch table"):
        installed_path(registry, "table", "1", hint, asset="table.tsv")
    with pytest.raises(NotInstalledError, match="^not here$"):
        installed_path(registry, "table", "1", "not here")
    registry.download("table", "1", show_progress=False)
    asked.clear()
    assert installed_path(registry, "table", "1", hint, asset="table.tsv").name == "table.tsv"
    assert not asked


def test_chain_and_genome_stores_are_named_by_their_pins(tmp_path):
    chain = pinned_chains()[0]
    pins = {chain.filename: dict(url=chain.url, sha256=chain.sha256, size=chain.size)}
    assert chain_registry(chain, tmp_path).inspect(chain.name, "1", verify_files=False).path == (
        str(get_data_dir(tmp_path) / "chains" / chain.name / content_digest(pins)))
    genome = pinned_genomes()[0]
    pins = {genome.assets[role]["url"].rsplit("/", 1)[-1]: dict(genome.assets[role])
            for role in ("gtf", "transcripts", "proteins")}
    assert genome_registry(genome, tmp_path).inspect(genome.name, "1", verify_files=False).path == (
        str(get_data_dir(tmp_path) / "genomes" / genome.name / content_digest(pins)))


def test_names_that_moved_here_are_still_importable_where_they_were():
    import orfdata.datasets
    import orfdata.net

    assert orfdata.datasets.get_data_dir is get_data_dir
    assert orfdata.net.DOWNLOAD_TIMEOUT == DOWNLOAD_TIMEOUT == 60
