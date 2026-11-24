"""Subtle cyan CLI help, table header, and status styling (#195)."""

import argparse
import io
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

import orfdata
from orfdata.cli import command_guide
from orfdata.cli_help import (
    HelpArgumentParser,
    HelpFormatter,
    classify_status,
    color_enabled,
    style_help,
    style_status,
    style_table_header,
)
from orfdata.cli_output import print_table
from orfdata.cli_parser import ArgumentParser
from orfdata.command_parsers import data_parser, query_parser

ANSI_ESCAPE = re.compile(r"\033\[[0-9;]*m")
SOURCE = str(Path(orfdata.__file__).resolve().parent.parent)


class TTYBuffer(io.StringIO):
    """In-memory text buffer that reports itself as an interactive TTY."""

    def isatty(self):
        return True


def strip_ansi(text):
    """Return ``text`` with all ANSI SGR escape sequences removed."""
    return ANSI_ESCAPE.sub("", text)


def test_color_enabled_respects_tty_no_color_and_dumb_term(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    non_tty = io.StringIO()
    tty = TTYBuffer()
    assert not color_enabled(non_tty)
    assert color_enabled(tty)

    monkeypatch.setenv("NO_COLOR", "1")
    assert not color_enabled(tty)

    monkeypatch.setenv("NO_COLOR", "")
    assert color_enabled(tty)

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert not color_enabled(tty)


def test_format_help_is_plain_and_print_help_styles_only_on_tty(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    parser = HelpArgumentParser(prog="orfdata-test", description="Test CLI")
    parser.add_argument("-o", "--output", help="Write output to ``FILE``")

    plain = parser.format_help()
    assert "\033[" not in plain
    assert "``" not in plain

    non_tty = io.StringIO()
    parser.print_help(file=non_tty)
    assert non_tty.getvalue() == plain
    assert "\033[" not in non_tty.getvalue()

    tty = TTYBuffer()
    parser.print_help(file=tty)
    styled = tty.getvalue()
    assert "\033[1musage:\033[0m" in styled
    assert "\033[1;36moptions:\033[0m" in styled
    assert "\033[36m-o\033[0m" in styled
    assert "\033[36m--output\033[0m" in styled
    assert strip_ansi(styled) == plain
    # Only normal cyan (36), bold cyan (1;36), bold (1), and reset (0) are used.
    codes = set(ANSI_ESCAPE.findall(styled))
    assert codes <= {"\033[0m", "\033[1m", "\033[36m", "\033[1;36m"}

    monkeypatch.setenv("NO_COLOR", "1")
    assert style_help(plain, tty) == plain

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert style_help(plain, tty) == plain


def test_orfdata_parsers_inherit_help_styling_and_subcommands_stay_plain_in_format_help(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    for parser in (command_guide("orfdata"), query_parser(), data_parser()):
        assert isinstance(parser, ArgumentParser)
        assert isinstance(parser, HelpArgumentParser)
        assert issubclass(parser.formatter_class, HelpFormatter)

        plain = parser.format_help()
        assert "\033[" not in plain

        tty = TTYBuffer()
        parser.print_help(file=tty)
        styled = tty.getvalue()
        assert "\033[1musage:\033[0m" in styled
        assert "\033[1;36moptions:\033[0m" in styled
        assert "\033[36m--help\033[0m" in styled
        assert strip_ansi(styled) == plain


def test_help_argument_parser_disables_python_314_builtin_color(monkeypatch):
    seen_kwargs = {}

    def fake_init(self, *args, **kwargs):
        seen_kwargs.update(kwargs)

    monkeypatch.setattr(sys, "version_info", (3, 14, 0, "final", 0))
    monkeypatch.setattr(argparse.ArgumentParser, "__init__", fake_init)

    HelpArgumentParser(prog="orfdata")
    assert seen_kwargs["color"] is False
    assert seen_kwargs["formatter_class"] is HelpFormatter


def test_style_table_header_and_status_helpers(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    tty = TTYBuffer()
    non_tty = io.StringIO()

    table = "KIND     NAME\ncatalog  bundled\n"
    assert style_table_header(table, non_tty) == table
    styled_table = style_table_header(table, tty)
    assert styled_table == "\033[1;36mKIND     NAME\033[0m\ncatalog  bundled\n"
    assert style_table_header("", tty) == ""

    assert classify_status("available") == "plain"
    assert classify_status("OK") == "plain"
    assert classify_status("missing") == "dim"
    assert classify_status("not bundled") == "dim"
    assert classify_status("-") == "dim"
    assert classify_status("corrupt") == "warning"
    assert classify_status("inaccessible") == "warning"

    # Normal values stay plain (no green, no red).
    assert style_status("OK", tty) == "OK"
    assert style_status("available", tty) == "available"
    # Absent / not applicable -> dim.
    assert style_status("missing", tty) == "\033[2mmissing\033[0m"
    assert style_status("n/a", tty) == "\033[2mn/a\033[0m"
    # Warnings -> yellow.
    assert style_status("corrupt", tty) == "\033[33mcorrupt\033[0m"
    assert style_status("inaccessible", tty) == "\033[33minaccessible\033[0m"
    # Explicit overrides and non-TTY guard.
    assert style_status("custom", tty, level="dim") == "\033[2mcustom\033[0m"
    assert style_status("custom", tty, level="warning") == "\033[33mcustom\033[0m"
    assert style_status("missing", tty, level="plain") == "missing"
    assert style_status("missing", non_tty) == "missing"
    assert style_status("corrupt", non_tty) == "corrupt"

    with pytest.raises(ValueError, match="Unknown status level"):
        style_status("missing", tty, level="green")


def test_print_table_styles_header_and_status_on_tty_without_breaking_alignment(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    monkeypatch.setenv("COLUMNS", "100")

    headings = ("KIND", "NAME", "STATUS", "ACQUIRE")
    rows = [
        ("catalog", "bundled", "available", "-"),
        ("catalog", "mouse-extra", "missing", "build"),
        ("dataset", "broken-set", "corrupt", "download"),
    ]

    non_tty = io.StringIO()
    print_table(headings, rows, stream=non_tty)
    plain = non_tty.getvalue()
    assert "\033[" not in plain

    tty = TTYBuffer()
    print_table(headings, rows, stream=tty)
    styled = tty.getvalue()
    lines = styled.splitlines()
    assert lines[0].startswith("\033[1;36m") and lines[0].endswith("\033[0m")
    assert "\033[" not in lines[1]  # "available" stays plain.
    assert "\033[2mmissing\033[0m" in lines[2]
    assert "\033[33mcorrupt\033[0m" in lines[3]
    # Stripping ANSI yields byte-for-byte identical alignment to plain output.
    assert strip_ansi(styled) == plain

    monkeypatch.setenv("NO_COLOR", "1")
    tty_no_color = TTYBuffer()
    print_table(headings, rows, stream=tty_no_color)
    assert tty_no_color.getvalue() == plain


def test_subprocess_cli_help_is_plain_when_piped_and_respects_no_color():
    env = dict(os.environ, PYTHONPATH=SOURCE, TERM="xterm-256color")
    env.pop("NO_COLOR", None)
    result = subprocess.run(
        [sys.executable, "-c", "from orfdata.cli import run; run(['--help'])"],
        capture_output=True,
        text=True,
        check=True,
        env=env,
    )
    assert "\033[" not in result.stdout
    assert "usage: orfdata" in result.stdout
    assert "commands:" in result.stdout
