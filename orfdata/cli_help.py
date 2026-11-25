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

"""Readable argparse help, table headers, and status values via ``datacache.cli_help``."""

import argparse

from datacache.cli_help import (
    OPTION_PATTERN,
    color_enabled,
    style_help,
    style_table_header,
)
from datacache.cli_help import (
    HelpArgumentParser as _BaseHelpArgumentParser,
)
from datacache.cli_help import (
    HelpFormatter as _BaseHelpFormatter,
)
from datacache.cli_help import (
    classify_status as _classify_status,
)
from datacache.cli_help import (
    style_status as _style_status,
)

# ORFdata renders "missing" dim (not yet installed) alongside other absent states.
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


def classify_status(text):
    """Classify an ORFdata status label into ``"dim"``, ``"warning"``, or ``"plain"``."""
    return _classify_status(
        text,
        dim_values=DIM_STATUS_VALUES,
        warning_values=WARNING_STATUS_VALUES,
    )


def style_status(text, stream, *, level=None):
    """Style a table status value on color terminals using dim or yellow only."""
    return _style_status(
        text,
        stream=stream,
        level=level,
        dim_values=DIM_STATUS_VALUES,
        warning_values=WARNING_STATUS_VALUES,
    )


class HelpFormatter(_BaseHelpFormatter):
    """Shared formatter keeping subcommand listings compact."""

    def _format_action(self, action):
        if isinstance(action, argparse._SubParsersAction._ChoicesPseudoAction):
            return argparse.RawDescriptionHelpFormatter._format_action(self, action)
        return super()._format_action(action)


class HelpArgumentParser(_BaseHelpArgumentParser):
    """ArgumentParser using :class:`HelpFormatter` and shared cyan help styling."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", HelpFormatter)
        super().__init__(*args, **kwargs)
