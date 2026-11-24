"""Shared human previews and explicit, complete machine-readable exports."""

import json
import shutil
import sys
import textwrap
from itertools import islice, zip_longest

import pandas as pd

from .cli_help import classify_status, style_status, style_table_header


def records(frame):
    """A data frame's rows as JSON-ready dicts, with its blanks as nulls."""
    return json.loads(frame.to_json(orient="records"))


def print_json(payload, stream=None):
    """Write one machine-readable JSON document."""
    stream = sys.stdout if stream is None else stream
    json.dump(payload, stream, indent=2, default=str)
    stream.write("\n")


def write_json(payload, path):
    """Write one machine-readable JSON document to a file.

    Parameters
    ----------
    payload : dict or list
        The document; values JSON has no type for are written as strings.
    path : str or Path
        The file to create or replace.
    """
    with open(path, "w", encoding="utf-8") as stream:
        print_json(payload, stream=stream)


def print_run_summary(summary, *, json_output=False, path=None, human=None):
    """Report the summary of a command whose result is a file or a stream of records.

    ``proteome`` and ``lift`` write a FASTA or a catalogue, so their summary is
    a second output. The complete summary goes to ``path`` when one is given,
    and to standard output as JSON when ``json_output`` is set, as ``--json``
    does in every other command. Otherwise a human summary goes to standard
    error, which leaves standard output to the records.

    Parameters
    ----------
    summary : dict
        The complete summary.
    json_output : bool
        Write the summary as one JSON document to standard output.
    path : str or Path, optional
        Also write the summary there as JSON.
    human : dict, optional
        What the human summary shows instead of ``summary``.
    """
    if path is not None:
        write_json(summary, path)
    if json_output:
        print_json(summary)
    else:
        print_fields(summary if human is None else human, stream=sys.stderr)


def _width():
    return max(40, min(120, shutil.get_terminal_size(fallback=(100, 24)).columns))


def _text(value):
    if value is None or (not isinstance(value, (dict, list, tuple)) and pd.isna(value)):
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (list, tuple)):
        return ", ".join(map(str, value)) or "-"
    return str(value) or "-"


def print_table(
    headings,
    rows,
    *,
    empty="No matching entries",
    stream=None,
    limit=25,
    total=None,
    hidden=None,
    nowrap=(),
):
    """Wrap human columns to the terminal; consume at most limit + 1 rows.

    A zero limit explicitly requests all rows. Exports bypass this renderer.

    Parameters
    ----------
    headings : sequence of str
        Column header labels.
    rows : iterable of sequence
        Row values aligned with ``headings``.
    empty : str, optional
        Message printed when ``rows`` is empty.
    stream : file-like, optional
        Output stream (defaults to ``sys.stdout``).
    limit : int, optional
        Maximum rows to render; ``0`` renders every row.
    total : int, optional
        Total unclipped row count when ``rows`` has already been pre-sliced
        (such as by kind-balanced preview selection).
    hidden : dict of (str, int), optional
        Optional breakdown of omitted rows (e.g. ``{"dataset": 4}``) rendered
        in the single shared truncation footer.
    nowrap : sequence of (int or str), optional
        Column indices or case-insensitive heading names (such as ``("NAME",)``)
        whose tokens must never split mid-token.
    """
    stream = sys.stdout if stream is None else stream
    if total is None and hasattr(rows, "__len__"):
        total = len(rows)
    values = list(rows if limit == 0 else islice(rows, limit + 1))
    more = limit > 0 and len(values) > limit
    if more:
        values = values[:limit]
    values = [tuple(_text(field) for field in row) for row in values]
    if not values:
        print(empty, file=stream)
        return
    nowrap_names = {item.strip().upper() for item in nowrap if isinstance(item, str)}
    nowrap_columns = {
        index for index, heading in enumerate(headings)
        if index in nowrap or heading.strip().upper() in nowrap_names
    }
    widths = [len(heading) for heading in headings]
    for value in values:
        widths = [max(width, len(field)) for width, field in zip(widths, value, strict=True)]
    # Keep identifiers intact when they fit; shrink prose and space-separated
    # lists first. Columns listed in ``nowrap`` keep their widest token intact
    # so names like ``dataset:search.fasta`` never split mid-token.
    minimum = []
    for index in range(len(widths)):
        widest_word = max(len(word) for row in [tuple(headings), *values] for word in row[index].split())
        minimum.append(widest_word if index in nowrap_columns else min(32, widest_word))
    available = _width() - 2 * (len(widths) - 1)
    while sum(widths) > available and max(widths) > 3:
        candidates = [
            index for index in range(len(widths))
            if widths[index] > minimum[index] and index not in nowrap_columns
        ]
        if not candidates:
            if nowrap_columns:
                break
            candidates = [index for index in range(len(widths)) if widths[index] > 3]
        if not candidates:
            break
        longest = max(candidates, key=widths.__getitem__)
        widths[longest] -= 1
    status_columns = {
        index for index, heading in enumerate(headings)
        if heading.strip().upper() == "STATUS"
    }
    for row_index, row in enumerate([tuple(headings), *values]):
        cells = [
            textwrap.wrap(
                field,
                width=width,
                break_on_hyphens=False,
                break_long_words=col_index not in nowrap_columns,
            ) or [""]
            for col_index, (field, width) in enumerate(zip(row, widths, strict=True))
        ]
        status_levels = [
            classify_status(row[col_index]) if row_index > 0 and col_index in status_columns else "plain"
            for col_index in range(len(widths))
        ]
        for line in zip_longest(*cells, fillvalue=""):
            if row_index == 0:
                rendered = "  ".join(
                    field.ljust(width) for field, width in zip(line, widths, strict=True)
                ).rstrip()
                rendered = style_table_header(rendered, stream)
            else:
                parts = []
                for col_index, (field, width) in enumerate(zip(line, widths, strict=True)):
                    padding = " " * max(0, width - len(field))
                    cell = (
                        style_status(field, stream, level=status_levels[col_index])
                        if field and status_levels[col_index] != "plain"
                        else field
                    )
                    parts.append(cell + padding)
                rendered = "  ".join(parts).rstrip()
            print(rendered, file=stream)
    if hidden:
        more_kinds = ", ".join(f"{count} {kind}" for kind, count in hidden.items())
        shown = f"{len(values)} of {total}" if total is not None else f"{len(values)}"
        print(
            f"\nShowing {shown} entries; not shown: {more_kinds}. "
            "Use --limit 0 for all rows or --json for a complete export.",
            file=stream,
        )
    elif more:
        shown = f"{limit} of {total}" if total is not None else f"the first {limit}"
        print(f"\nShowing {shown} entries. Use --limit 0 for all rows or --json for a complete export.", file=stream)


def print_fields(result, stream=None):
    """A bounded, indented summary; full nested metadata remains available via --json."""
    stream = sys.stdout if stream is None else stream
    omitted = False
    lines = []

    def visit(key, value, depth=0):
        nonlocal omitted
        prefix = "  " * depth + (str(key).replace("_", " ") if depth == 0 else str(key)) + ":"
        if isinstance(value, dict):
            if not value:
                lines.append(prefix + " none")
            elif depth >= 2:
                lines.append(prefix + f" {len(value)} fields")
                omitted = True
            else:
                lines.append(prefix)
                for name, item in list(value.items())[:12]:
                    visit(name, item, depth + 1)
                if len(value) > 12:
                    lines.append("  " * (depth + 1) + f"... {len(value) - 12} more fields")
                    omitted = True
        elif isinstance(value, (list, tuple)):
            if any(isinstance(item, (dict, list, tuple)) for item in value):
                lines.append(prefix + f" {len(value)} entries")
                for index, item in enumerate(value[:2], 1):
                    visit(index, item, depth + 1)
                omitted |= len(value) > 2
            else:
                lines.append(prefix + " " + (_text(value[:5]) if value else "none"))
                if len(value) > 5:
                    lines.append("  " * (depth + 1) + f"... {len(value) - 5} more entries")
                    omitted = True
        else:
            lines.append(prefix + " " + _text(value))

    for key, value in result.items():
        visit(key, value)
    wrapped = [part for line in lines for part in
               textwrap.wrap(line, width=_width(), subsequent_indent="  ", break_on_hyphens=False)]
    print("\n".join(wrapped[:60]), file=stream)
    if omitted or len(wrapped) > 60:
        print("\nUse --json for complete metadata.", file=stream)


def print_result(result, *, json_output=False, columns=None, limit=25, stream=None):
    """Render a management result with shared field/table conventions."""
    if json_output:
        print_json(result, stream=stream)
    elif isinstance(result, dict):
        print_fields(result, stream=stream)
    else:
        rows = list(result)
        columns = columns or (tuple(rows[0]) if rows else ())
        print_table(tuple(key.replace("_", " ").upper() for key in columns),
                    [tuple(_field(row, key) for key in columns) for row in rows], limit=limit, stream=stream)


def _field(row, path):
    for part in path.split("."):
        row = row.get(part) if isinstance(row, dict) else None
    return row


def print_record_stream(rows, *, format, columns, limit=25):
    """Stream complete JSON/JSONL, or show a bounded selection of named fields."""
    if format == "jsonl":
        for row in rows:
            print(json.dumps(row, default=str))
    elif format == "json":
        sys.stdout.write("[")
        separator = ""
        for row in rows:
            sys.stdout.write(separator + "\n" + json.dumps(row, default=str))
            separator = ","
        sys.stdout.write("\n]\n")
    else:
        references = {}

        def field(row, path):
            value = _field(row, path)
            if isinstance(value, dict) and "taxon_id" in value:
                key = json.dumps(value, sort_keys=True)
                if key not in references:
                    references[key] = (f"r{len(references) + 1}", value)
                return references[key][0]
            if isinstance(value, float):
                return f"{value:.6g}"
            if isinstance(value, str):
                if path.endswith(("sample_id", "pair_id", ".identifier")) and len(value) > 24:
                    return value[:9] + "..." + value[-12:]
                if path.endswith(("quantity", "status", "attribution")):
                    return value.replace("_", " ")
            return value

        print_table(tuple(label for label, _ in columns),
                    (tuple(field(row, path) for _, path in columns) for row in rows), limit=limit)
        if references:
            print("\nReferences (full identities in --json):")
            print_table(("REF", "SPECIES", "BUILD", "ANNOTATION", "VERSION", "SEQUENCE PIN"),
                        [(label, ref.get("species"), ref.get("genome_build"), ref.get("annotation_name"),
                          ref.get("annotation_version"), (ref.get("sequence_sha256") or "")[:12])
                         for label, ref in references.values()], limit=0)
        print("\nUse --json for complete records or --format jsonl for streaming export.")
        print("Previews abbreviate long IDs and show six significant digits; exports retain full precision.")


def manifest_summary(manifest):
    """Select useful discovery fields; complete pins and processing metadata are exports."""
    summary = {key: manifest[key] for key in ("catalog_id", "dataset_id", "version", "title", "assay_kinds")
               if key in manifest}
    summary["references"] = [" / ".join(str(ref[key]) for key in
                                         ("species", "genome_build", "annotation_name", "annotation_version")
                                         if ref.get(key) is not None)
                             for ref in manifest.get("references", ())]
    summary["source_assets"] = [source["asset_id"] for source in manifest.get("sources", ())]
    summary["coverage"] = {key: value for key, value in manifest.get("coverage", {}).items()
                           if not isinstance(value, (dict, list, tuple))}
    return summary
