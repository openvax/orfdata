"""One table lists the top-level commands in the guide and runs them."""

import pytest

from orfdata import cli
from orfdata.cli import COMMAND_ALIASES, COMMAND_TABLE, COMMANDS, Command, command_guide, run


@pytest.fixture
def recorded(monkeypatch):
    """Replace every command's function with one that records how it was called."""
    calls = []
    for command in COMMAND_TABLE.values():
        if command.function is None:
            continue
        module, name = command.function.split(":")

        def record(argv, prog, *, recorded_as=command.function):
            calls.append((recorded_as, argv, prog))

        monkeypatch.setattr(f"{module}.{name}", record)
    return calls


def test_every_command_but_help_loads_a_function_and_builds_its_parser():
    for name, command in COMMAND_TABLE.items():
        if name == "help":
            assert command.function is None and command.parser is None
            with pytest.raises(ValueError):
                command.load()
        else:
            assert callable(command.load()), name
            assert command.parser(prog=f"orfdata {name}").prog == f"orfdata {name}"
    assert set(COMMAND_ALIASES.values()) <= set(COMMAND_TABLE)


def test_the_guide_lists_the_table_in_order(capsys):
    assert COMMANDS == {name: command.description for name, command in COMMAND_TABLE.items()}
    command_guide("orfdata").print_help()
    words = " ".join(capsys.readouterr().out.split())  # Descriptions wrap to the terminal.
    positions = [words.index(f" {name} {command.description} ") for name, command in COMMAND_TABLE.items()]
    assert positions == sorted(positions)


@pytest.mark.parametrize("name", [name for name in COMMAND_TABLE if name not in ("help", "query")])
def test_each_command_runs_its_function_with_its_own_arguments(name, recorded):
    run([name, "first", "--second"], prog="tool")
    assert recorded == [(COMMAND_TABLE[name].function, ["first", "--second"], f"tool {name}")]


@pytest.mark.parametrize("alias", sorted(COMMAND_ALIASES))
def test_an_alias_runs_its_command_under_the_name_typed(alias, recorded):
    run([alias, "list"])
    assert recorded == [(COMMAND_TABLE[COMMAND_ALIASES[alias]].function, ["list"], f"orfdata {alias}")]


def test_query_alone_takes_the_catalogue_loader(monkeypatch):
    calls = []
    monkeypatch.setattr(cli, "query", lambda argv, prog, *, catalog_loader: calls.append(
        (argv, prog, catalog_loader)))
    loader = object()
    run(["query", "--gene", "TP53"], catalog_loader=loader)
    run(["--gene", "TP53"], prog="tool", catalog_loader=loader)
    assert calls == [(["--gene", "TP53"], "orfdata query", loader), (["--gene", "TP53"], "tool query", loader)]


@pytest.mark.parametrize("argv", [["nonsense"], ["help", "help"], ["help", "nonsense"]])
def test_an_unknown_command_is_a_usage_error(argv, recorded, capsys):
    with pytest.raises(SystemExit) as stopped:
        run(argv)
    assert stopped.value.code == 2
    assert f"Unknown command {argv[-1]!r}" in capsys.readouterr().err
    assert recorded == []


def test_a_command_without_a_function_says_so():
    with pytest.raises(ValueError, match="orfdata.cli.run"):
        Command("Answered elsewhere").load()
