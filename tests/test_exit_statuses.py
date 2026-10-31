"""Exit statuses: 1 for failures, 2 for command-line mistakes, 3 for missing installations."""

import json
import socket

import pytest
from test_measurements import definition as dataset_definition
from test_measurements import write_dataset

from orfdata import CatalogRegistry, DatasetRegistry, NotInstalledError
from orfdata.cli import run


def _exit(argv, capsys):
    with pytest.raises(SystemExit) as error:
        run(argv)
    output = capsys.readouterr()
    return error.value.code, output.err


@pytest.fixture
def offline(monkeypatch, tmp_path):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path / "cache"))
    monkeypatch.setattr(socket.socket, "connect", lambda *a, **k: pytest.fail("network"))


@pytest.mark.parametrize("content", [
    {"schema": "orfdata.catalog.v1"},  # No metadata.
    {"schema": "orfdata.catalog.v1", "metadata": {"coordinates": "1-based-inclusive", "sources": ["x"]},
     "records": [{"unexpected": 1}]},  # Records without their required fields.
    [],  # Not a JSON object.
])
def test_a_malformed_catalogue_file_is_a_failure_naming_the_file(content, offline, tmp_path, capsys):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(content))
    for argv in (["query", "--catalog", str(path)],
                 ["lift", "--catalog", str(path), "--to", "GRCh37", "--output", str(tmp_path / "out.json")],
                 ["proteome", "--species", "human", "--release", "115", "--nc-orf-sources", str(path)]):
        status, err = _exit(argv, capsys)
        assert status == 1, argv
        assert f"{path} is not a valid ORF catalogue" in err and "Traceback" not in err
        assert "--help' for usage" not in err


def test_missing_pyensembl_annotation_is_not_installed(offline, monkeypatch, tmp_path, capsys):
    from orfdata.proteome import _require_installed

    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "pyensembl"))
    status, err = _exit(["proteome", "--species", "human", "--release", "111"], capsys)
    assert status == 3 and "annotation is missing; install it first with: pyensembl install" in err

    class Genome:
        reference_name = "GRCh38"

        def __init__(self, status):
            self.status = status

        def inspect_data(self):
            return {"annotation": self.status}

    for state in ("missing", "incomplete", "not indexed"):
        with pytest.raises(NotInstalledError):
            _require_installed(Genome(state))
    with pytest.raises(ValueError) as damaged:
        _require_installed(Genome("invalid"))
    assert not isinstance(damaged.value, NotInstalledError)
    _require_installed(Genome("indexed"))


@pytest.mark.parametrize("argv", [
    ["query", "--species", "mouse", "--require-reviewed"],
    ["query", "--species", "mouse", "--genome-build", "GRCm39", "--region", "bad"],
    ["query", "--species", "mouse", "--ms-methods", "unknown-kind"],
    ["proteome", "--nc-orf-sources", "missing.json", "--genome", "t2t-chm13", "--species", "human"],
    ["proteome", "--nc-orf-sources", "missing.json"],
])
def test_command_line_mistakes_exit_2_before_anything_is_read(argv, offline, capsys):
    status, err = _exit(argv, capsys)
    assert status == 2 and "--help' for usage" in err and "not installed" not in err


@pytest.fixture
def installed_dataset(tmp_path, experiment, monkeypatch):
    from orfdata import dataset_cli

    original = write_dataset(tmp_path / "source.sqlite", experiment)
    registry = DatasetRegistry(dataset_definition(original), cache_root=tmp_path / "cache")
    registry.download("fixture", show_progress=False)
    monkeypatch.setattr(dataset_cli, "DatasetRegistry", lambda definitions=None, **kw: registry)
    return registry


@pytest.mark.parametrize("argv", [
    ["datasets", "query", "fixture", "--reference", "bogus"],
    ["datasets", "query", "fixture", "--target", "T", "--kind", "gene", "--namespace", "N", "--reference", "bogus"],
    ["datasets", "pairs", "fixture", "--basis", "same_specimen", "--reference", "bogus"],
    ["datasets", "pairs", "fixture", "--basis", "same_specimen", "--left-target", "T", "--left-kind", "gene",
     "--left-namespace", "N", "--left-reference", "bogus"],
    ["datasets", "inspect-reference", "fixture", "--reference", "bogus"],
    ["datasets", "query-reference", "fixture", "--target", "T", "--asset", "bogus"],
])
def test_an_unknown_dataset_reference_is_a_usage_error_listing_the_choices(argv, installed_dataset, capsys):
    status, err = _exit(argv, capsys)
    assert status == 2 and "'bogus' is not a" in err and "choose from" in err


def test_restricted_chains_need_accept_terms_in_every_command(offline, capsys):
    from orfdata import list_data
    from orfdata.liftover import pinned_chains

    chain = next(c for c in pinned_chains() if c.restricted)
    name = next(r.name for r in list_data(kind="chain") if r.restricted)
    for argv in (["chains", "download", chain.source_build, chain.target_build],
                 ["data", "chains", chain.source_build, chain.target_build],
                 ["data", "fetch", name]):
        status, err = _exit(argv, capsys)
        assert status == 2 and "--accept-terms" in err, argv


def test_install_hints_name_kind_version_cache_root_and_custom_registries(offline, tmp_path, capsys, experiment):
    root = tmp_path / "elsewhere"
    status, err = _exit(["catalogs", "--cache-root", str(root), "query", "riborf-mouse"], capsys)
    assert status == 3
    assert f"CatalogRegistry(cache_root='{root}').download('riborf-mouse', '3')" in err
    assert f"orfdata data fetch riborf-mouse --kind catalog --cache-root {root} --version 3" in err

    with pytest.raises(NotInstalledError) as built:
        CatalogRegistry(cache_root=root).open_build("riborf-mouse")  # Published, so a build is explicit.
    assert f"orfdata catalogs --cache-root {root} build riborf-mouse --version 3" in str(built.value)

    original = write_dataset(tmp_path / "source.sqlite", experiment)
    custom = DatasetRegistry(dataset_definition(original), cache_root=root)
    with pytest.raises(NotInstalledError) as missing:
        custom.open("fixture")
    message = str(missing.value)
    assert "registry.download('fixture'" in message and "orfdata datasets --registry FILE" in message
    assert "data fetch" not in message
