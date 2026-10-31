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
    (["genomes", "install", "t2t-chm13"], "orfdata.genomes_cli.install_genome", _genome()),
    (["data", "genome", "t2t-chm13"], "orfdata.genomes.install_genome", _genome()),
    (["chains", "download", "GRCh38", "GRCh37"], "orfdata.liftover_cli.download_chain", "chain"),
    (["data", "chains", "GRCh38", "GRCh37"], "orfdata.liftover.download_chain", "chain"),
    (["chains", "download"], "orfdata.liftover_cli.download_chains", dict(downloaded={}, skipped={})),
])
def test_every_acquisition_command_has_the_same_timeout(argv, target, result, monkeypatch, tmp_path, capsys):
    calls = _capture(monkeypatch, target, result)
    run([*argv, "--json"])
    run([argv[0], "--cache-root", str(tmp_path), *argv[1:], "--timeout", "5", "--json"])
    capsys.readouterr()
    assert [call["download_options"]["timeout"] for call in calls] == [60, 5]
    assert calls[1]["cache_root"] == str(tmp_path)


def test_dataset_downloads_no_longer_wait_indefinitely(monkeypatch, capsys):
    calls = []

    class Registry:
        def __init__(self, definitions=None, *, cache_root=None):
            pass

        def download(self, name, version, **options):
            calls.append(options)
            return {}

    monkeypatch.setattr("orfdata.dataset_cli.DatasetRegistry", Registry)
    run(["datasets", "download", "fixture", "--json"])
    capsys.readouterr()
    assert calls[0]["timeout"] == 60


@pytest.mark.parametrize("argv", [
    ["viral", "--cache-root", "{root}", "inventory", "download"],
    ["viral", "inventory", "--cache-root", "{root}", "download"],
    ["viral", "inventory", "download", "--cache-root", "{root}"],
])
def test_viral_inventory_reads_cache_root_at_any_level(argv, monkeypatch, tmp_path, capsys):
    calls = _capture(monkeypatch, "orfdata.viral_cli.download_viral_inventory", tmp_path / "inventory.tsv")
    run([part.format(root=tmp_path) for part in argv] + ["--json"])
    capsys.readouterr()
    assert calls[0]["cache_root"] == str(tmp_path) and calls[0]["timeout"] == 60


def test_timeouts_must_be_positive_and_finite_or_none(capsys, monkeypatch):
    for value in ("0", "-1", "nan", "inf", "1e400", "soon"):
        with pytest.raises(SystemExit):
            run(["data", "fetch", "bundled", "--timeout", value])
        assert "positive number of seconds, or none" in capsys.readouterr().err
    calls = _capture(monkeypatch, "orfdata.liftover_cli.download_chain", "chain")
    run(["chains", "download", "GRCh38", "GRCh37", "--timeout", "none", "--json"])
    assert calls[0]["download_options"]["timeout"] is None


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
    from orfdata.net import DOWNLOAD_TIMEOUT
    from orfdata.viral_panel import download_viral_references

    assert DOWNLOAD_TIMEOUT == 60
    assert download_viral_references.__kwdefaults__["timeout"] == DOWNLOAD_TIMEOUT
    viral = []

    class Bundle:
        def to_dict(self):
            return {}

    monkeypatch.setattr("orfdata.viral_panel.download_viral_references",
                        lambda **options: viral.append(options) or Bundle())
    run(["data", "viral", "--json"])
    assert viral[0]["timeout"] == 60
    calls = []

    class Registry:
        def __init__(self, definitions=None, *, cache_root=None):
            pass

        def build(self, name, version, *, download_options, **options):
            calls.append(download_options)
            return type("Built", (), dict(path="p", counts={}, manifest=type("M", (), dict(
                dataset_id=name, version="1"))))

        def download_reference(self, name, version, **options):
            calls.append(options)
            return {}

    monkeypatch.setattr("orfdata.dataset_cli.DatasetRegistry", Registry)
    run(["datasets", "build", "fixture", "--json"])
    run(["datasets", "download-reference", "fixture", "--timeout", "7", "--json"])
    capsys.readouterr()
    assert [call["timeout"] for call in calls] == [60, 7]


def test_suppressed_defaults_do_not_leak_into_help():
    parser = ArgumentParser(prog="tool")
    parser.add_subparsers(dest="command").add_parser("leaf")
    share_cache_root(parser)
    assert argparse.SUPPRESS not in parser.format_help()
