"""Inspect GDC viral inventory coverage and download viral source records."""

import argparse
import csv
import json
import sys
from dataclasses import asdict

from .viral_coverage import PANEL, download_viral_inventory, viral_contexts, viral_panel_coverage
from .viral_panel import VIRAL_COLLECTION, download_viral_references, index_viral_references
from .viral_subsets import viral_subsets


def run(argv=None, prog="orfdata viral"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    coverage = commands.add_parser("coverage", help="List all GDC panel rows and coding-reference status")
    selection = coverage.add_mutually_exclusive_group()
    selection.add_argument("--available", action="store_true", help="Show rows with bundled source CDSs")
    selection.add_argument("--missing", action="store_true", help="Show rows not yet bundled in ORFdata")
    coverage.add_argument("--virus", help="Exact GDC inventory abbreviation, e.g. HPV18")
    coverage.add_argument("--search", help="Case-insensitive name, label or source-accession search")
    coverage.add_argument("--context", action="append",
                          help="Select references in this context; repeat for a union")
    coverage.add_argument("--exclude-context", action="append",
                          help="Exclude references in this context; repeat as needed")
    coverage.add_argument("--context-match", choices=("any", "all"), default="any",
                          help="Match any selected context (default) or require all")
    output = coverage.add_mutually_exclusive_group()
    output.add_argument("--format", choices=("table", "json", "tsv"), default="table")
    output.add_argument("--json", dest="format", action="store_const", const="json",
                        help="Write machine-readable JSON")

    contexts = commands.add_parser("contexts", help="List selectable viral-reference contexts")
    contexts.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    subsets = commands.add_parser("subsets", help="List named viral-reference subsets")
    subsets.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    inventory = commands.add_parser("inventory", help="Manage the published GDC inventory metadata")
    inventory_commands = inventory.add_subparsers(dest="inventory_command", required=True)
    download = inventory_commands.add_parser("download", help="Download the pinned inventory via DataCache")
    download.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE or the default")
    download.add_argument("--force", action="store_true", help="Replace the cached copy from the pinned source")
    download.add_argument("--timeout", type=float, default=60, help="Per-attempt timeout in seconds")

    download_references = commands.add_parser(
        "download", help="Compatibility alias for 'orfdata download viral'")
    download_references.add_argument("--species", default="human",
                                     help="Target analysis species (currently human)")
    download_references.add_argument(
        "--collection", default=VIRAL_COLLECTION,
        help="Reference collection (currently gdc-grch38-d1-vd1)")
    download_references.add_argument("--context", action="append",
                                     help="Build an index for this context; repeat for a union")
    download_references.add_argument("--exclude-context", action="append",
                                     help="Exclude this context from the derived index; repeat as needed")
    download_references.add_argument("--context-match", choices=("any", "all"), default="any",
                                     help="Match any selected context (default) or require all")
    download_references.add_argument("--subset", action="append",
                                     help="Build an index for this named viral subset; repeat for a union")
    download_references.add_argument("--cache-root",
                                     help="Cache directory instead of ORFDATA_CACHE or the default")
    download_references.add_argument("--force", action="store_true",
                                     help="Replace or repair the cached viral reference bundle")
    download_references.add_argument("--timeout", type=float, default=120,
                                     help="Per-attempt NCBI timeout in seconds")
    download_references.add_argument("--json", action="store_true",
                                     help="Write machine-readable JSON")

    args = parser.parse_args(argv)
    if args.command == "contexts":
        rows = viral_contexts()
        if args.json:
            json.dump([asdict(row) for row in rows], sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            headings = ("CONTEXT", "ALIASES", "REFERENCES", "DESCRIPTION")
            values = [(row.name, ",".join(row.aliases) or "-", str(row.references), row.description)
                      for row in rows]
            widths = [max(len(heading), *(len(row[index]) for row in values))
                      for index, heading in enumerate(headings)]
            print("  ".join(heading.ljust(widths[index])
                            for index, heading in enumerate(headings)))
            for row in values:
                print("  ".join(value.ljust(widths[index])
                                for index, value in enumerate(row)))
        return
    if args.command == "subsets":
        rows = viral_subsets()
        if args.json:
            json.dump(rows, sys.stdout, indent=2)
            sys.stdout.write("\n")
        else:
            headings = ("SUBSET", "REFERENCES", "DESCRIPTION")
            values = [(row["name"], str(row["references"]), row["description"]) for row in rows]
            widths = [max(len(heading), *(len(row[index]) for row in values))
                      for index, heading in enumerate(headings)]
            print("  ".join(heading.ljust(widths[index])
                            for index, heading in enumerate(headings)))
            for row in values:
                print("  ".join(value.ljust(widths[index])
                                for index, value in enumerate(row)))
        return
    if args.command == "coverage":
        available = True if args.available else False if args.missing else None
        try:
            rows = viral_panel_coverage(
                available=available, virus=args.virus, search=args.search,
                contexts=args.context, exclude_contexts=args.exclude_context,
                context_match=args.context_match)
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

    if args.command == "inventory":
        try:
            path = download_viral_inventory(cache_root=args.cache_root, force=args.force,
                                            timeout=args.timeout)
        except (ImportError, ValueError, OSError) as error:
            parser.error(str(error))
        print(path)
        return

    try:
        bundle = download_viral_references(
            species=args.species, collection=args.collection, cache_root=args.cache_root,
            force=args.force, timeout=args.timeout)
        if args.context or args.exclude_context or args.subset:
            contexts = args.context or (["cancer-transcriptome-screening"]
                                        if args.exclude_context else None)
            bundle = index_viral_references(
                contexts=contexts, exclude_contexts=args.exclude_context,
                context_match=args.context_match, subsets=args.subset,
                species=args.species, collection=args.collection,
                cache_root=args.cache_root, force=args.force, timeout=args.timeout)
    except (ImportError, RuntimeError, ValueError, OSError) as error:
        parser.error(str(error))
    if args.json:
        json.dump(bundle.to_dict(), sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(bundle.path)
