"""Shared command-line output: aligned tables and machine-readable JSON.

Every ORFdata command prints a readable table or key/value block by default
and the same content as JSON with ``--json``, so listings stay comparable
across resource kinds.
"""

import json
import sys


def print_json(payload, stream=None):
    """Write one machine-readable JSON document."""
    stream = sys.stdout if stream is None else stream
    json.dump(payload, stream, indent=2, default=str)
    stream.write("\n")


def print_table(headings, rows, *, empty="No matching entries", stream=None):
    """Print left-aligned columns sized to their contents; empty fields show as ``-``."""
    stream = sys.stdout if stream is None else stream
    values = [tuple("-" if field is None or field == "" else str(field) for field in row)
              for row in rows]
    if not values:
        print(empty, file=stream)
        return
    widths = [len(heading) for heading in headings]
    for value in values:
        widths = [max(width, len(field)) for width, field in zip(widths, value, strict=True)]
    for row in [tuple(headings), *values]:
        print("  ".join(field.ljust(width)
                        for field, width in zip(row, widths, strict=True)).rstrip(), file=stream)


def print_fields(result, stream=None):
    """Print one ``key: value`` line per field, with nested values as compact JSON."""
    stream = sys.stdout if stream is None else stream
    for key, value in result.items():
        if isinstance(value, (dict, list, tuple)):
            value = json.dumps(value, sort_keys=True, default=str)
        print(f"{key}: {value}", file=stream)
