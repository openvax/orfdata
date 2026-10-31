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
    assert "registry.download('fixture'" in message and "the same --registry file after 'orfdata datasets'" in message
    assert f"on the command line, adding the same --registry file after 'orfdata datasets': orfdata datasets " \
           f"--cache-root {root} download fixture" in message
    assert "data fetch" not in message and "FILE" not in message


def test_install_hints_are_not_requoted(offline, tmp_path, capsys):
    root = tmp_path / "my cache"
    status, err = _exit(["catalogs", "--cache-root", str(root), "query", "riborf-mouse"], capsys)
    assert status == 3
    assert f"or on the command line: orfdata data fetch riborf-mouse --kind catalog --cache-root '{root}'" in err
    assert "'orfdata" not in err


def _gzip_cut(path, payload):
    import gzip

    data = gzip.compress(json.dumps(payload).encode())
    path.write_bytes(data[: len(data) // 2])


def test_every_catalogue_loader_names_a_malformed_or_truncated_file(offline, tmp_path, capsys):
    from orfdata.catalog import load_uorf_catalog

    truncated, corrupt, shapeless = tmp_path / "cut.json.gz", tmp_path / "bad.json.gz", tmp_path / "shape.json"
    _gzip_cut(truncated, {"schema": "orfdata.catalog.v1", "metadata": {}, "records": []})
    corrupt.write_bytes(b"\x1f\x8b" + b"not gzip at all" * 4)
    shapeless.write_text(json.dumps({"schema": "isovar.uorf_evidence.v1"}))
    for path in (truncated, corrupt, shapeless):
        for loader in (None, load_uorf_catalog):  # isovar-uorf-evidence injects load_uorf_catalog.
            with pytest.raises(SystemExit) as error:
                run(["query", "--catalog", str(path)], **({} if loader is None else {"catalog_loader": loader}))
            err = capsys.readouterr().err
            assert error.value.code == 1 and f"{path} is not a valid ORF catalogue" in err, (path, loader)


def test_species_selection_needs_the_species_catalogue_loader(offline, capsys):
    from orfdata.catalog import load_uorf_catalog

    with pytest.raises(SystemExit) as error:
        run(["query", "--species", "mouse"], catalog_loader=load_uorf_catalog)
    err = capsys.readouterr().err
    assert error.value.code == 2 and "--species is not available here" in err


@pytest.mark.parametrize("argv", [
    ["query", "canine-mammary-stroma", "--reference", "bogus"],
    ["query", "canine-mammary-stroma", "--target", "X", "--kind", "gene"],
    ["query", "canine-mammary-stroma", "--mapping-status", "unique"],
    ["query", "canine-mammary-stroma", "--mapped"],
    ["pairs", "canine-mammary-stroma", "--basis", "same_specimen", "--left-target", "X"],
    ["pairs", "canine-mammary-stroma", "--basis", "same_specimen", "--right-reference", "bogus"],
])
def test_dataset_query_mistakes_exit_2_on_an_empty_cache(argv, offline, capsys):
    status, err = _exit(["datasets", *argv], capsys)
    assert status == 2 and "not installed" not in err


def test_bad_reference_lists_the_registered_choices_without_opening(offline, capsys):
    status, err = _exit(["datasets", "query", "canine-mammary-stroma", "--reference", "bogus"], capsys)
    assert status == 2 and "choose from ref:" in err


@pytest.mark.parametrize("argv", [
    ["query", "--region", "chr1:20-10", "--genome-build", "GRCh38"],
    ["query", "--region", "chr1:0-10", "--genome-build", "GRCh38"],
    ["query", "--position", "chr1:0", "--genome-build", "GRCh38"],
    ["query", "--species", "mouse", "--min-protein-length", "0"],
    ["proteome", "--species", "human", "--release", "115", "--min-protein-length", "5"],
    ["proteome", "--species", "human", "--release", "115", "--nc-orf-sources", "ms", "--min-protein-length", "0"],
])
def test_out_of_range_values_are_usage_errors_before_loading(argv, offline, capsys):
    status, err = _exit(argv, capsys)
    assert status == 2 and "--help' for usage" in err


def test_missing_chains_for_a_supplied_catalogue_are_not_installed(tmp_path):
    from orfdata.catalog import load_catalog
    from orfdata.proteome import _onto
    from orfdata.reference import ReferenceContext

    catalog = load_catalog()  # The bundled GRCh38 catalogue.
    target = ReferenceContext(9606, "homo_sapiens", "GRCh37", "ensembl", "75")
    with pytest.raises(NotInstalledError) as error:
        _onto(catalog, target, cache_root=tmp_path / "my cache")
    message = str(error.value)
    assert "needs the chain files GRCh38_to_GRCh37" in message
    assert f"orfdata chains --cache-root '{tmp_path / 'my cache'}' download --species homo_sapiens" in message


@pytest.mark.parametrize("group, label", [("catalogs", "catalogue"), ("datasets", "dataset")])
def test_a_registry_file_with_a_json_syntax_error_is_named(group, label, offline, tmp_path, capsys):
    path = tmp_path / "registry.json"
    path.write_text("{bad")
    status, err = _exit([group, "--registry", str(path), "list"], capsys)
    assert status == 1 and f"Invalid {label} registry {path}" in err


def test_an_undownloaded_published_catalogue_with_a_build_points_to_it(tmp_path):
    from test_catalog_registry import definition as catalog_definition
    from test_catalog_registry import setup_release

    manifest, source, _ = setup_release(tmp_path)
    registry = CatalogRegistry(catalog_definition(manifest), cache_root=tmp_path / "cache")
    with pytest.raises(NotInstalledError, match="is not installed"):
        registry.open(manifest.catalog_id)
    registry.build(manifest.catalog_id, source_paths={"source.json": source})
    with pytest.raises(NotInstalledError) as error:
        registry.open(manifest.catalog_id)
    assert "a local build is installed" in str(error.value) and "--built" in str(error.value)
