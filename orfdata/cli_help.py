# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Readable argparse help, table headers, and status values with optional terminal styling."""

import argparse
import os
import re
import shutil
import sys

OPTION_PATTERN = re.compile(r"(?<![\w-])(--[a-zA-Z][\w-]*|-[a-zA-Z])\b")

# Status tokens rendered dim (absent / not applicable) or yellow (warnings) on color terminals.
# Normal values, including "OK" and "available", stay plain: no green, no red.
DIM_STATUS_VALUES = frozenset({
    "-",
    "absent",
    "missing",
    "n/a",
    "none",
    "not applicable",
    "not bundled",
    "not installed",
})
WARNING_STATUS_VALUES = frozenset({
    "corrupt",
    "damaged",
    "inaccessible",
    "incomplete",
    "invalid",
    "stale",
    "warning",
})

__all__ = [
    "DIM_STATUS_VALUES",
    "HelpArgumentParser",
    "HelpFormatter",
    "OPTION_PATTERN",
    "WARNING_STATUS_VALUES",
    "classify_status",
    "color_enabled",
    "style_help",
    "style_status",
    "style_table_header",
]


def color_enabled(stream):
    """Return whether ``stream`` supports optional ANSI terminal styling.

    Styling is enabled only when ``stream`` is a TTY, ``NO_COLOR`` is unset or
    empty, and ``TERM`` is not ``"dumb"``.
    """
    return bool(
        getattr(stream, "isatty", lambda: False)()
        and not os.environ.get("NO_COLOR")
        and os.environ.get("TERM") != "dumb"
    )


def style_help(text, stream):
    """Style headings, option flags, and ``usage:`` only on color-capable terminals.

    Section headings become bold cyan, option flags become cyan, and ``usage:``
    becomes bold. Everything else stays plain. Returns ``text`` unchanged when
    :func:`color_enabled` is false.
    """
    if not color_enabled(stream):
        return text

    lines = []
    for line in text.splitlines(keepends=True):
        if line.rstrip().endswith(":") and not line[0].isspace():
            ending = "\n" if line.endswith("\n") else ""
            lines.append("\033[1;36m" + line.rstrip("\n") + "\033[0m" + ending)
        else:
            line = OPTION_PATTERN.sub(lambda match: "\033[36m" + match[0] + "\033[0m", line)
            if line.startswith("usage:"):
                line = line.replace("usage:", "\033[1musage:\033[0m", 1)
            lines.append(line)
    return "".join(lines)


def style_table_header(text, stream):
    """Style the header line of an aligned table in bold cyan on color terminals.

    Returns ``text`` unchanged when ``text`` is empty or :func:`color_enabled`
    is false.
    """
    if not text or not color_enabled(stream):
        return text
    lines = text.splitlines(keepends=True)
    ending = "\n" if lines[0].endswith("\n") else ""
    lines[0] = "\033[1;36m" + lines[0].rstrip("\n") + "\033[0m" + ending
    return "".join(lines)


def classify_status(text):
    """Classify a status label into ``"dim"``, ``"warning"``, or ``"plain"``.

    Absent / not-applicable labels return ``"dim"``, warning labels return
    ``"warning"``, and normal values (including ``"OK"`` and ``"available"``)
    return ``"plain"``.
    """
    normalized = text.strip().lower()
    if normalized in DIM_STATUS_VALUES:
        return "dim"
    if normalized in WARNING_STATUS_VALUES:
        return "warning"
    return "plain"


def style_status(text, stream, *, level=None):
    """Style a table status value on color terminals using dim or yellow only.

    Absent or not-applicable values use dim (``\\033[2m``) and warning values
    use yellow (``\\033[33m``). Normal values, including ``"OK"`` and
    ``"available"``, stay plain (no green, no red). Returns ``text`` unchanged
    when ``text`` is empty or :func:`color_enabled` is false.

    Parameters
    ----------
    text : str
        The status text to format.
    stream : file-like
        The output stream checked with :func:`color_enabled`.
    level : {None, "dim", "warning", "plain"}, optional
        Explicit styling level. When None (default), inferred with
        :func:`classify_status`.
    """
    if not text or not color_enabled(stream):
        return text
    resolved = classify_status(text) if level is None else level
    if resolved == "dim":
        return f"\033[2m{text}\033[0m"
    if resolved == "warning":
        return f"\033[33m{text}\033[0m"
    if resolved == "plain":
        return text
    raise ValueError(
        f"Unknown status level {level!r}; expected None, 'dim', 'warning', or 'plain'"
    )


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Limit reading width and separate options without adding table borders."""

    def __init__(self, prog):
        width = max(40, min(100, shutil.get_terminal_size().columns - 2))
        super().__init__(prog, max_help_position=32, width=width)

    def _split_lines(self, text, width):
        # RST literals are useful in generated docs but noisy in terminal help.
        return super()._split_lines(text.replace("``", ""), width)

    def _format_action(self, action):
        text = super()._format_action(action)
        if not text or isinstance(action, argparse._SubParsersAction._ChoicesPseudoAction):
            return text
        return text + "\n"


class HelpArgumentParser(argparse.ArgumentParser):
    """Keep ``format_help()`` plain; add color only when printing to a terminal."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", HelpFormatter)
        if sys.version_info >= (3, 14):
            # Apply the same styling on all supported Python versions.
            kwargs["color"] = False
        super().__init__(*args, **kwargs)

    def print_help(self, file=None):
        stream = sys.stdout if file is None else file
        self._print_message(style_help(self.format_help(), stream), stream)
