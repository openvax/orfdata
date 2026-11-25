"""Options shared by command groups mean the same thing wherever they're placed."""

import argparse
import json
from types import SimpleNamespace

import pytest

from orfdata.cli import run
from orfdata.cli_parser import ArgumentParser, share_cache_root


def test_cache_root_is_read_before_or_after_nested_subcommands():
    parser = ArgumentParser(prog="tool")
    groups = parser.add_subparsers(dest="command", required=True)
    leaf = groups.add_parser("leaf")
    leaf.add_argument("--cache-root", help="its own help")
    nested = groups.add_parser("group").add_subparsers(dest="inner", required=True)
    nested.add_parser("inner")
    share_cache_root(parser)
    for argv in (["--cache-root", "X", "leaf"], ["leaf", "--cache-root", "X"],
                 ["--cache-root", "X", "group", "inner"], ["group", "--cache-root", "X", "inner"],
                 ["group", "inner", "--cache-root", "X"]):
        assert parser.parse_args(argv).cache_root == "X", argv
    assert parser.parse_args(["leaf"]).cache_root is None
    assert parser.parse_args(["--cache-root", "X", "leaf", "--cache-root", "Y"]).cache_root == "Y"


@pytest.mark.parametrize("group, command", [
    ("data", ["list", "--kind", "catalog", "--json"]),
    ("catalogs", ["inspect", "riborf-mouse", "--json"]),
    ("datasets", ["inspect", "wang-human-tissue-abundance", "--json"]),
])
def test_group_commands_accept_cache_root_in_either_place(group, command, tmp_path, capsys):
    for argv in ([group, "--cache-root", str(tmp_path), *command], [group, *command, "--cache-root", str(tmp_path)]):
        run(argv)
        output = json.loads(capsys.readouterr().out)
        rows = output if isinstance(output, list) else [output]
        assert all(row["path"].startswith(str(tmp_path)) for row in rows if row.get("acquire") or "acquire" not in row)


def _capture(monkeypatch, target, result=None):
    calls = []

    def capture(*args, **kwargs):
        calls.append(kwargs)
        return result

    monkeypatch.setattr(target, capture)
    return calls


def _genome():
    return SimpleNamespace(reference_name="ref", annotation_name="fixture", annotation_version="1",
                           transcript_ids=lambda: ())


@pytest.mark.parametrize("argv, target, result", [
    (["genomes", "install", "t2t-chm13"], "orfdata.genomes.install_genome", _genome()),
    (["data", "genome", "t2t-chm13"], "orfdata.genomes.install_genome", _genome()),
    (["chains", "download", "GRCh38", "GRCh37"], "orfdata.liftover.download_chain", "chain"),
    (["data", "chains", "GRCh38", "GRCh37"], "orfdata.liftover.download_chain", "chain"),
    (["chains", "download", "--species", "human"], "orfdata.liftover.download_chain", "chain"),
])
def test_every_acquisition_command_has_the_same_timeout(argv, target, result, monkeypatch, tmp_path, capsys):
    from orfdata.resource_kinds import get_resource_kind

    calls = _capture(monkeypatch, target, result)
    # Mark inspection as available right after acquisition so fetch_data's post-check succeeds.
    kind_name = "genome" if "genome" in argv[0] or "genome" in argv[1] else "chain"
    handler = get_resource_kind(kind_name)
    real_inspect = handler.inspect
    from dataclasses import replace

    state = {"done": False}

    def fake_inspect(name, version=None, *, cache_root=None, verify_files=False, options=None):
        row = real_inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        if state["done"]:
            return replace(row, status="available", origin="download")
        return replace(row, status="missing", origin=None)

    def wrapped_capture(*args, **kwargs):
        calls.append(kwargs)
        state["done"] = True
        return result

    calls.clear()
    monkeypatch.setattr(target, wrapped_capture)
    monkeypatch.setattr(handler, "inspect", fake_inspect)

    state["done"] = False
    run([*argv, "--json"])
    first_timeout = calls[-1]["download_options"]["timeout"]
    calls.clear()
    state["done"] = False
    run([argv[0], "--cache-root", str(tmp_path), *argv[1:], "--timeout", "5", "--json"])
    capsys.readouterr()
    assert [first_timeout, calls[-1]["download_options"]["timeout"]] == [60, 5]
    assert calls[-1]["cache_root"] == str(tmp_path)


def test_dataset_downloads_no_longer_wait_indefinitely(monkeypatch, capsys):
    from orfdata.data_resources import fetch_data as real_fetch

    calls = []

    def fake_fetch(name, **options):
        calls.append(options)
        return {"kind": "dataset", "name": name, "path": "/tmp/x", "status": "available"}

    monkeypatch.setattr("orfdata.data_resources.fetch_data", fake_fetch)
    monkeypatch.setattr("orfdata.datasets.DatasetRegistry.require_published", lambda self, name, version=None: "1")
    monkeypatch.setattr(
        "orfdata.datasets.DatasetRegistry.inspect",
        lambda self, name, version=None, verify_files=True: SimpleNamespace(status="available", files={}, path="/tmp/x"),
    )
    run(["datasets", "download", "copy-yeast-abundance", "--json"])
    capsys.readouterr()
    assert calls[0]["timeout"] == 60
    assert real_fetch is not None


@pytest.mark.parametrize("argv", [
    ["viral", "--cache-root", "{root}", "inventory", "download"],
    ["viral", "inventory", "--cache-root", "{root}", "download"],
    ["viral", "inventory", "download", "--cache-root", "{root}"],
])
def test_viral_inventory_reads_cache_root_at_any_level(argv, monkeypatch, tmp_path, capsys):
    calls = _capture(monkeypatch, "orfdata.viral_cli.fetch_data", {"path": str(tmp_path / "inventory.tsv")})
    run([part.format(root=tmp_path) for part in argv] + ["--json"])
    capsys.readouterr()
    assert calls[0]["cache_root"] == str(tmp_path) and calls[0]["timeout"] == 60


def test_timeouts_must_be_positive_and_finite_or_none(capsys, monkeypatch):
    for value in ("0", "-1", "nan", "inf", "1e400", "soon"):
        with pytest.raises(SystemExit):
            run(["data", "fetch", "bundled", "--timeout", value])
        assert "positive number of seconds, or none" in capsys.readouterr().err
    calls = _capture(monkeypatch, "orfdata.data_resources.fetch_data", {"path": "chain"})
    run(["chains", "download", "GRCh38", "GRCh37", "--timeout", "none", "--json"])
    assert calls[0]["timeout"] is None


def test_cache_root_is_refused_where_nothing_reads_it(capsys):
    for argv in (["viral", "coverage", "--cache-root", "X"], ["viral", "--cache-root", "X", "contexts"],
                 ["viral", "--cache-root", "X", "subsets"]):
        with pytest.raises(SystemExit):
            run(argv)
        assert "--cache-root has no effect" in capsys.readouterr().err


def test_viral_summary_reads_cache_root_before_or_after(monkeypatch, tmp_path, capsys):
    calls = _capture(monkeypatch, "orfdata.viral_cli.viral_reference_summary", {})
    run(["viral", "--cache-root", str(tmp_path), "summary", "--json"])
    run(["viral", "summary", "--cache-root", str(tmp_path), "--json"])
    capsys.readouterr()
    assert [call["cache_root"] for call in calls] == [str(tmp_path)] * 2


def test_viral_and_dataset_downloads_share_the_default(monkeypatch, capsys):
    from orfdata.cache import DOWNLOAD_TIMEOUT
    from orfdata.viral_panel import download_viral_references

    assert DOWNLOAD_TIMEOUT == 60
    assert download_viral_references.__kwdefaults__["timeout"] == DOWNLOAD_TIMEOUT
    calls = []

    def fake_fetch(name, **options):
        calls.append(options)
        return {"kind": "viral", "name": name, "path": "p", "status": "available"}

    monkeypatch.setattr("orfdata.data_resources.fetch_data", fake_fetch)
    monkeypatch.setattr(
        "orfdata.datasets.DatasetRegistry.open_build",
        lambda self, name, version=None, verify_files=True: SimpleNamespace(
            path="p", counts={}, manifest=SimpleNamespace(dataset_id=name, version="1"),
        ),
    )
    run(["data", "viral", "--json"])
    run(["datasets", "build", "peptideatlas-dog-peptides", "--json"])
    run(["datasets", "download-reference", "pfizer-atlas-dog-protein", "--timeout", "7", "--json"])
    capsys.readouterr()
    assert [call["timeout"] for call in calls] == [60, 60, 7]


def test_suppressed_defaults_do_not_leak_into_help():
    parser = ArgumentParser(prog="tool")
    parser.add_subparsers(dest="command").add_parser("leaf")
    share_cache_root(parser)
    assert argparse.SUPPRESS not in parser.format_help()
