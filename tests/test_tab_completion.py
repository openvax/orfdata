"""Tab completion through argcomplete, and the cheap imports that keep each Tab press fast.

Every Tab press runs ``orfdata`` up to argcomplete, so importing the command
line and building a command's parser must not import pandas, DataCache or the
modules that run commands. Times aren't asserted, because they vary with the
machine; these tests check what is imported instead.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import orfdata
from orfdata.cli import COMMAND_ALIASES, COMMAND_TABLE, completion_parser, completion_words

# What a Tab press must not import: the scientific stack, DataCache, and ORFdata's data and command modules.
HEAVY = ("pandas", "numpy", "pyarrow", "datacache", "pyensembl", "biotite", "ahocorasick", "openpyxl", "yaml",
         "orfdata.catalog", "orfdata.catalogs", "orfdata.datasets", "orfdata.proteome", "orfdata.viral",
         "orfdata.viral_panel", "orfdata.viral_coverage", "orfdata.liftover", "orfdata.genomes",
         "orfdata.cli_output", "orfdata.download_cli", "orfdata.catalog_cli", "orfdata.dataset_cli",
         "orfdata.genomes_cli", "orfdata.liftover_cli", "orfdata.proteome_cli", "orfdata.viral_cli")
# The checkout under test, ahead of any installed copy in the subprocesses.
SOURCE = str(Path(orfdata.__file__).resolve().parent.parent)
# What the 'orfdata' console script runs.
SCRIPT = "import sys\nfrom orfdata.cli import run\nsys.exit(run())"


def run_python(code, **environ):
    """Run ``code`` in a fresh interpreter importing this checkout, and return its standard output."""
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                          env=dict(os.environ, PYTHONPATH=SOURCE, **environ)).stdout


def modules_loaded_by(code):
    """Return the names in ``sys.modules`` after running ``code`` in a fresh interpreter."""
    return set(json.loads(run_python(code + "\nimport json, sys\nprint(json.dumps(sorted(sys.modules)))")))


def completions(line, tmp_path):
    """Return what one Tab press at the end of ``line`` offers, as the shell's argcomplete hook asks for it.

    A lone completion comes with the space the shell would type after it; it is removed here.
    """
    output = tmp_path / "completions.txt"
    run_python(SCRIPT, _ARGCOMPLETE="1", COMP_LINE=line, COMP_POINT=str(len(line)), _ARGCOMPLETE_IFS="\n",
               _ARGCOMPLETE_STDOUT_FILENAME=str(output))
    return {word.rstrip(" ") for word in output.read_text().split("\n")} - {""}


def long_options(parser):
    """Return the ``--`` options of ``parser`` itself, not of its subcommands."""
    return {option for action in parser._actions for option in action.option_strings if option.startswith("--")}


@pytest.mark.parametrize("code", [
    "import orfdata",
    "import orfdata.cli",
    # What a Tab press does at most: build the guide with one command's parser under it.
    "from orfdata.cli import COMMAND_TABLE, completion_parser\n"
    "for name in COMMAND_TABLE:\n"
    "    completion_parser([name])",
])
def test_the_command_line_imports_nothing_heavy(code):
    assert sorted(modules_loaded_by(code) & set(HEAVY)) == []


def test_the_entry_module_is_marked_for_global_completion():
    with open(Path(SOURCE) / "orfdata" / "cli.py", "rb") as handle:
        assert b"# PYTHON_ARGCOMPLETE_OK" in handle.read(1024)


@pytest.mark.parametrize("line, expected", [
    ("orfdata ", set(COMMAND_TABLE) | {"--version", "--help"}),
    ("orfdata --", {"--version", "--help"}),
    ("orfdata help ", set(COMMAND_TABLE) - {"help"}),
    ("orfdata data ", {"list", "status", "info", "path", "fetch", "viral", "catalog", "dataset", "reference",
                       "genome", "chains"}),
    ("orfdata data list --kind ", {"viral", "catalog", "dataset", "reference", "genome", "chain", "records"}),
    ("orfdata downloads fetch --acc", {"--accept-terms"}),
    ("orfdata viral ", {"coverage", "contexts", "subsets", "summary", "completeness", "inventory"}),
    ("orfdata viral inventory ", {"download"}),
    ("orfdata query --ev", {"--evidence"}),
    ("orfdata catalogs query NAME --form", {"--format"}),
    ("orfdata datasets pairs NAME --left-s", {"--left-scale", "--left-status"}),
    ("orfdata proteome --isoforms ", {"all", "longest"}),
])
def test_tab_offers_commands_options_and_choices(line, expected, tmp_path):
    offered = completions(line, tmp_path)
    assert expected <= offered, expected - offered


@pytest.mark.parametrize("name", [name for name in COMMAND_TABLE if name != "help"] + sorted(COMMAND_ALIASES))
def test_tab_offers_every_option_of_each_command(name, tmp_path):
    command = COMMAND_TABLE[COMMAND_ALIASES.get(name, name)]
    expected = long_options(command.parser(prog=f"orfdata {name}"))
    assert expected and expected <= completions(f"orfdata {name} --", tmp_path)


def test_completion_words_are_those_before_the_cursor_after_the_program():
    line = "orfdata data fetch --kind cat"
    environ = dict(COMP_LINE=line, COMP_POINT=str(len(line)), _ARGCOMPLETE="1")
    assert completion_words(environ) == ["data", "fetch", "--kind"]
    environ.update(COMP_LINE="python -m orfdata.cli viral ", COMP_POINT="28", _ARGCOMPLETE="3")
    assert completion_words(environ) == ["viral"]


def test_completion_parser_attaches_only_the_typed_command():
    def subparsers(parser):
        return next(action for action in parser._actions if action.choices and "data" in action.choices).choices

    guide = subparsers(completion_parser([]))
    assert all(not parser._subparsers for parser in guide.values())
    typed = subparsers(completion_parser(["--version", "download", "fetch"]))
    assert typed["download"].prog == "orfdata download" and typed["download"]._subparsers is not None
    assert not typed["data"]._subparsers and not typed["viral"]._subparsers
