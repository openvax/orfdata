"""Inspect GDC viral inventory coverage and download its pinned metadata."""

import argparse
import csv
import json
import sys
from dataclasses import asdict

from .viral_coverage import PANEL, download_viral_inventory, viral_panel_coverage


def run(argv=None, prog="orfdata viral"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    coverage = commands.add_parser("coverage", help="List all GDC panel rows and coding-reference status")
    selection = coverage.add_mutually_exclusive_group()
    selection.add_argument("--available", action="store_true", help="Show rows with bundled source CDSs")
    selection.add_argument("--missing", action="store_true", help="Show rows not yet bundled in ORFdata")
    coverage.add_argument("--virus", help="Exact GDC inventory abbreviation, e.g. HPV18")
    coverage.add_argument("--search", help="Case-insensitive name, label or source-accession search")
    coverage.add_argument("--format", choices=("table", "json", "tsv"), default="table")

    inventory = commands.add_parser("inventory", help="Manage the published GDC inventory metadata")
    inventory_commands = inventory.add_subparsers(dest="inventory_command", required=True)
    download = inventory_commands.add_parser("download", help="Download the pinned inventory via DataCache")
    download.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE or the default")
    download.add_argument("--force", action="store_true", help="Replace the cached copy from the pinned source")
    download.add_argument("--timeout", type=float, default=60, help="Per-attempt timeout in seconds")

    args = parser.parse_args(argv)
    if args.command == "coverage":
        available = True if args.available else False if args.missing else None
        try:
            rows = viral_panel_coverage(available=available, virus=args.virus, search=args.search)
            all_rows = viral_panel_coverage()
        except ValueError as error:
            parser.error(str(error))
        if args.format == "json":
            json.dump(dict(panel=PANEL, total=len(all_rows), available=sum(r.available for r in all_rows),
                           shown=len(rows), entries=[dict(asdict(r), available=r.available) for r in rows]),
                      sys.stdout, indent=2)
            sys.stdout.write("\n")
        elif args.format == "tsv":
            writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
            writer.writerow(("abbreviation", "name", "status", "coding_models", "refseq",
                             "genbank", "pave_id", "bundled_reference_id"))
            for row in rows:
                writer.writerow((row.abbreviation, row.name, row.status, row.coding_models,
                                 row.refseq, row.genbank, row.pave_id, row.bundled_reference_id))
        else:
            print(f"{PANEL}: {sum(r.available for r in all_rows)}/{len(all_rows)} "
                  f"with bundled source CDSs; showing {len(rows)}", file=sys.stdout)
            if rows:
                widths = (max(12, *(len(r.abbreviation) for r in rows)), 13, 6)
                print(f"{'Virus':<{widths[0]}}  {'Status':<{widths[1]}}  {'CDSs':>{widths[2]}}  "
                      "GDC reference / source", file=sys.stdout)
                for row in rows:
                    print(f"{row.abbreviation:<{widths[0]}}  "
                          f"{'available' if row.available else 'not bundled':<{widths[1]}}  "
                          f"{row.coding_models or '-':>{widths[2]}}  "
                          f"{row.bundled_reference_id or row.refseq or row.genbank or '-'}",
                          file=sys.stdout)
        return

    try:
        path = download_viral_inventory(cache_root=args.cache_root, force=args.force,
                                        timeout=args.timeout)
    except (ImportError, ValueError, OSError) as error:
        parser.error(str(error))
    print(path)
