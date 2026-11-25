"""One discovery surface lists every resource kind and refuses selectors it cannot honour."""

import json

import pytest

from orfdata import RESOURCE_STATES, list_data
from orfdata.cli import run
from orfdata.downloads import collection_state, default_collections
from orfdata.install import _selected_collections

KINDS = {"viral", "inventory", "catalog", "dataset", "reference", "genome", "chain", "records"}


def test_every_kind_is_listed_with_an_acquisition_command():
    rows = list_data()
    assert {row.kind for row in rows} == KINDS
    # A bundled resource needs no command; everything acquirable names one.
    assert all(row.command for row in rows if row.acquire is not None)
    assert all(row.acquire is not None or row.installed for row in rows)
    assert all(row.status in RESOURCE_STATES for row in rows)
    assert all(row.default for row in rows)
    assert len(list_data(all_versions=True)) > len(rows)
    assert [row.name for row in list_data(kind="viral")] == ["gdc-viruses"]
    assert list_data(kind="all") == rows


def test_filters_select_without_acquiring_anything():
    human = list_data(species="human")
    assert human and all("homo_sapiens" in row.species for row in human)
    assert all(row.installed is False for row in list_data(installed=False))
    assert all(row.installed for row in list_data(installed=True))
    cancer = list_data(contexts=["cancer"])
    assert [row.kind for row in cancer] == ["viral"]
    assert cancer[0].to_dict()["installed"] is cancer[0].installed


def test_unknown_selectors_are_errors_not_empty_listings():
    # A typo must not read as "nothing is available".
    with pytest.raises(ValueError, match="Unknown viral context"):
        list_data(contexts=["cancer-transcriptom"])
    with pytest.raises(ValueError, match="Unknown data kind"):
        list_data(kind="virus")
    with pytest.raises(ValueError, match="Unknown species alias"):
        list_data(species="unicorn")
    with pytest.raises(ValueError, match="context_match"):
        list_data(contexts=["cancer"], context_match="either")


def test_collections_are_declared_per_species(tmp_path):
    assert default_collections("human") == ("gdc-viruses",)
    assert default_collections("mouse") == ()
    assert collection_state("gdc-viruses", species="human", cache_root=str(tmp_path)) == "missing"
    assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "missing"
    with pytest.raises(ValueError, match="Unknown viral collection"):
        collection_state("gdc-rna", species="human", cache_root=str(tmp_path))


def test_setup_chooses_collections_like_catalogues():
    available = ("gdc-viruses",)
    assert _selected_collections(available, None, "human") == available
    assert _selected_collections(available, "all", "human") == available
    assert _selected_collections(available, ["all"], "human") == available
    assert _selected_collections(available, "none", "human") == ()
    assert _selected_collections(available, [], "human") == ()
    assert _selected_collections(available, ["gdc-viruses"], "human") == available
    assert _selected_collections((), None, "mouse") == ()
    with pytest.raises(ValueError, match="no reference collection"):
        _selected_collections(available, ["gdc-rna"], "human")


def test_cli_listings_are_readable_and_machine_readable(capsys):
    run(["download", "list", "--limit", "0"])
    output = capsys.readouterr().out
    assert "KIND" in output and "gdc-viruses" in output and "orfdata data fetch" in output

    run(["download", "list", "--kind", "viral", "--json"])
    rows = json.loads(capsys.readouterr().out)
    assert [row["name"] for row in rows] == ["gdc-viruses"]
    assert rows[0]["contexts"] == ["cancer-transcriptome-screening"]
    assert rows[0]["installed"] is False

    run(["download", "status", "--kind", "viral"])
    assert "No matching resources" in capsys.readouterr().out

    with pytest.raises(SystemExit):
        run(["download", "list", "--context", "nonsense"])
    assert "Unknown viral context" in capsys.readouterr().err


def test_setup_reports_what_each_collection_needed(monkeypatch, tmp_path):
    from orfdata import install_species, liftover
    from orfdata.resource_kinds import get_resource_kind

    viral_kind = get_resource_kind("viral")
    states = {"gdc-viruses": "missing"}
    installed = []
    real_inspect = viral_kind.inspect

    def fake_inspect(name, version=None, *, cache_root=None, verify_files=False, options=None):
        row = real_inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        from dataclasses import replace

        return replace(row, status=states[name], origin="build" if states[name] == "available" else None)

    monkeypatch.setattr(viral_kind, "inspect", fake_inspect)
    monkeypatch.setattr(
        viral_kind,
        "acquire",
        lambda name, version=None, **options: installed.append((name, options["force"])),
    )
    monkeypatch.setattr(liftover, "pinned_chains", lambda species=None: ())

    # Mock a single already-installed human catalogue so setup exercises only the viral collection.
    from orfdata import catalog as catalog_mod
    from orfdata.catalogs import CatalogRegistry

    registry = CatalogRegistry(cache_root=tmp_path)
    monkeypatch.setattr(catalog_mod, "catalog_choices", lambda reg=None: [])
    # Keep bundled human + one fake choice? Simpler: patch install_species chain list to one available chain.
    from types import SimpleNamespace

    fake_chain = SimpleNamespace(
        name="GRCh38_to_GRCh37",
        species="homo_sapiens",
        source_build="GRCh38",
        target_build="GRCh37",
        provider="Ensembl",
        terms="open",
        restricted=False,
    )
    monkeypatch.setattr(liftover, "pinned_chains", lambda species=None: (fake_chain,))
    chain_kind = get_resource_kind("chain")
    monkeypatch.setattr(
        chain_kind,
        "inspect",
        lambda name, version=None, **kw: SimpleNamespace(status="available", acquire="download"),
    )

    result = install_species("human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"])
    assert result["collections"] == {"gdc-viruses": "built"}
    assert installed == [("gdc-viruses", False)]

    states["gdc-viruses"] = "available"
    result = install_species("human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"])
    assert result["collections"] == {"gdc-viruses": "already installed"}

    # A forced replacement of a healthy collection is reported as replaced, not repaired.
    monkeypatch.setattr(chain_kind, "acquire", lambda *a, **k: None)
    result = install_species(
        "human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"], force=True,
    )
    assert result["collections"] == {"gdc-viruses": "replaced"}

    states["gdc-viruses"] = "damaged"
    result = install_species("human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"])
    assert result["collections"] == {"gdc-viruses": "repaired"}
    assert installed[-1] == ("gdc-viruses", True)

    # Unreadable files are refused, as data fetch refuses them, unless forced.
    from orfdata.install import InstallError

    states["gdc-viruses"] = "inaccessible"
    with pytest.raises(InstallError) as exc:
        install_species("human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"])
    assert exc.value.result["collections"]["gdc-viruses"].startswith("failed: gdc-viruses can't be read")
    result = install_species(
        "human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"], force=True,
    )
    assert result["collections"] == {"gdc-viruses": "repaired"}
    assert registry is not None


@pytest.mark.parametrize("error", [ImportError("missing required dependency"), RuntimeError("NCBI EFetch failed")])
def test_collection_failures_are_reported_including_missing_required_dependencies(monkeypatch, tmp_path, error):
    from dataclasses import replace
    from types import SimpleNamespace

    from orfdata import install_species, liftover
    from orfdata.install import InstallError
    from orfdata.resource_kinds import get_resource_kind

    fake_chain = SimpleNamespace(
        name="GRCh38_to_GRCh37",
        species="homo_sapiens",
        source_build="GRCh38",
        target_build="GRCh37",
        provider="Ensembl",
        terms="open",
        restricted=False,
    )
    monkeypatch.setattr(liftover, "pinned_chains", lambda species=None: (fake_chain,))
    chain_kind = get_resource_kind("chain")
    monkeypatch.setattr(
        chain_kind,
        "inspect",
        lambda name, version=None, **kw: SimpleNamespace(status="available", acquire="download"),
    )

    viral_kind = get_resource_kind("viral")
    real_inspect = viral_kind.inspect
    monkeypatch.setattr(
        viral_kind,
        "inspect",
        lambda name, version=None, *, cache_root=None, verify_files=False, options=None: replace(
            real_inspect(name, version, cache_root=cache_root, verify_files=False, options=options),
            status="missing",
            origin=None,
        ),
    )

    def refuse(*args, **options):
        raise error

    monkeypatch.setattr(viral_kind, "acquire", refuse)
    with pytest.raises(InstallError) as exc:
        install_species("human", cache_root=tmp_path, catalogs=["bundled"], collections=["gdc-viruses"])
    assert exc.value.result["collections"]["gdc-viruses"] == f"failed: {error}"
    assert issubclass(InstallError, ValueError)


def test_a_refused_request_is_not_retried(monkeypatch):
    from urllib.error import HTTPError

    from orfdata import net

    attempts = []

    def refuse(request, timeout=None):
        attempts.append(request.full_url)
        raise HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(net, "urlopen", refuse)
    monkeypatch.setattr(net.time, "sleep", lambda seconds: None)
    with pytest.raises(RuntimeError, match="was refused: HTTP 404"):
        net.fetch("https://example.invalid/record", timeout=5, description="a record")
    assert len(attempts) == 1  # A missing record will still be missing next time.


def test_a_transport_failure_is_retried(monkeypatch):
    from orfdata import net

    attempts = []

    def reset(request, timeout=None):
        attempts.append(request.full_url)
        raise ConnectionResetError("connection reset by peer")

    monkeypatch.setattr(net, "urlopen", reset)
    monkeypatch.setattr(net.time, "sleep", lambda seconds: None)
    with pytest.raises(RuntimeError, match="Could not read a record"):
        net.fetch("https://example.invalid/record", timeout=5, description="a record")
    assert len(attempts) == net.MAX_ATTEMPTS
