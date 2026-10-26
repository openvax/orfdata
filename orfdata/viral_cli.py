"""Inspect GDC viral coverage, contexts, subsets and installed reference layers.

Acquisition lives in 'orfdata download viral'; this command describes what is
published and what is installed.
"""

import argparse
import sys
from pathlib import Path

from .cli_output import print_json, print_table, records
from .viral_coverage import PANEL, download_viral_inventory, viral_contexts, viral_panel_coverage
from .viral_panel import viral_annotation_survey, viral_reference_summary
from .viral_subsets import viral_subsets


def _add_selection(parser):
    parser.add_argument("--context", action="append",
                        help="Select references in this context; repeat for a union")
    parser.add_argument("--exclude-context", action="append",
                        help="Exclude references in this context; repeat as needed")
    parser.add_argument("--context-match", choices=("any", "all"), default="any",
                        help="Match any selected context (default) or require all")


def run(argv=None, prog="orfdata viral"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    coverage = commands.add_parser("coverage", help="List all GDC human-screen rows and coding-reference status")
    selection = coverage.add_mutually_exclusive_group()
    selection.add_argument("--available", action="store_true", help="Show rows with bundled source CDSs")
    selection.add_argument("--missing", action="store_true", help="Show rows not yet bundled in ORFdata")
    coverage.add_argument("--virus", help="Exact GDC inventory abbreviation, e.g. HPV18")
    coverage.add_argument("--search", help="Case-insensitive name, label or source-accession search")
    _add_selection(coverage)
    output = coverage.add_mutually_exclusive_group()
    output.add_argument("--format", choices=("table", "json", "tsv"), default="table")
    output.add_argument("--json", dest="format", action="store_const", const="json",
                        help="Write machine-readable JSON")

    contexts = commands.add_parser("contexts", help="List selectable viral-reference contexts")
    contexts.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    subsets = commands.add_parser("subsets", help="List named viral-reference subsets")
    subsets.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    summary = commands.add_parser(
        "summary", help="Summarize the installed collection's genome, transcriptome and proteome")
    _add_selection(summary)
    summary.add_argument("--subset", action="append",
                         help="Count only this named viral subset; repeat for a union")
    summary.add_argument("--species", default="human", help="Target analysis species (currently human)")
    summary.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE or the default")
    summary.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    completeness = commands.add_parser(
        "completeness", help="Report how completely each installed genome is annotated")
    _add_selection(completeness)
    completeness.add_argument("--subset", action="append",
                              help="Survey only this named viral subset; repeat for a union")
    completeness.add_argument("--species", default="human",
                              help="Target analysis species (currently human)")
    completeness.add_argument("--cache-root",
                              help="Cache directory instead of ORFDATA_CACHE or the default")
    completeness.add_argument("--incomplete", action="store_true",
                              help="Show only genomes with a gap worth looking at")
    completeness.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    inventory = commands.add_parser("inventory", help="Manage the published GDC inventory metadata")
    inventory_commands = inventory.add_subparsers(dest="inventory_command", required=True)
    download = inventory_commands.add_parser("download", help="Download the pinned inventory via DataCache")
    download.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE or the default")
    download.add_argument("--force", action="store_true", help="Replace the cached copy from the pinned source")
    download.add_argument("--timeout", type=float, default=60, help="Per-attempt timeout in seconds")

    args = parser.parse_args(argv)
    if args.command == "contexts":
        rows = viral_contexts()
        if args.json:
            print_json(records(rows))
        else:
            print_table(("CONTEXT", "ALIASES", "REFERENCES", "DESCRIPTION"),
                        [(row.name, ",".join(row.aliases), row.references,
                          row.description)
                         for row in rows.itertuples(index=False)])
        return
    if args.command == "subsets":
        rows = viral_subsets()
        if args.json:
            print_json(rows)
        else:
            print_table(("SUBSET", "REFERENCES", "DESCRIPTION"),
                        [(row["name"], row["references"], row["description"]) for row in rows])
        return
    if args.command == "summary":
        try:
            result = viral_reference_summary(
                contexts=args.context, exclude_contexts=args.exclude_context,
                context_match=args.context_match, subsets=args.subset, species=args.species,
                cache_root=args.cache_root)
        except (ImportError, ValueError, OSError) as error:
            parser.error(str(error))
        _show_summary(result, args.json)
        return
    if args.command == "completeness":
        try:
            result = viral_annotation_survey(
                contexts=args.context, exclude_contexts=args.exclude_context,
                context_match=args.context_match, subsets=args.subset, species=args.species,
                cache_root=args.cache_root)
        except (ImportError, ValueError, OSError) as error:
            parser.error(str(error))
        _show_completeness(result, args.incomplete, args.json)
        return
    if args.command == "coverage":
        _show_coverage(args, parser)
        return

    try:
        path = download_viral_inventory(cache_root=args.cache_root, force=args.force,
                                        timeout=args.timeout)
    except (ImportError, ValueError, OSError) as error:
        parser.error(str(error))
    print(path)


def _show_summary(result, json_output):
    if json_output:
        print_json(result)
        return
    selection = result["selection"]
    chosen = [f"contexts={','.join(selection['contexts'])}" if selection["contexts"] else "",
              f"excluded={','.join(selection['excluded_contexts'])}"
              if selection["excluded_contexts"] else "",
              f"subsets={','.join(selection['subsets'])}" if selection["subsets"] else ""]
    print(f"{result['collection']} ({result['panel']}) for {result['target_species']} "
          f"taxon {result['target_taxon_id']}")
    print(f"version: {result['bundle_version']}")
    print(f"retrieved: {result['retrieved'] or '-'}")
    print(f"path: {result['path']}")
    print(f"selection: {', '.join(part for part in chosen if part) or 'whole collection'}")
    print()
    # The directory is on the path line above; the table names the files in it.
    print_table(("LAYER", "MOLECULE", "RECORDS", "FASTA", "TABLE"),
                [(layer["name"], layer["molecule"], layer["records"],
                  Path(layer["fasta"]).name, Path(layer["table"]).name)
                 for layer in result["layers"]])
    print(f"\nannotation: {Path(result['layers'][1]['annotation']).name}")
    print(f"{result['complete_transcript_count']} of {result['transcript_count']} coding "
          "transcripts have whole start and stop codons")
    sources = result["annotation_sources"]
    if sorted(sources) != ["genbank"]:
        print("annotation sources: " + ", ".join(f"{name} {count}"
                                                 for name, count in sorted(sources.items())))
    unannotated = result["genomes_without_transcripts"]
    if unannotated:
        print(f"{len(unannotated)} of {result['genome_count']} genomes have no annotated CDS "
              f"feature: {', '.join(unannotated)}")
    for name, reason in sorted(result["annotation_fill_skipped"].items()):
        print(f"  {name}: annotation could not be filled: {reason}")
    sources = result["genome_sources"]
    if sorted(sources) != [next(iter(sources), "genbank")]:
        print("genome sources: " + ", ".join(f"{name} {count}"
                                             for name, count in sorted(sources.items())))
    unlike = result["genomes_unlike_alternate"]
    if unlike:
        print(f"{len(unlike)} of {result['genome_count']} genomes are not the alternate "
              f"reference their inventory row names: {', '.join(unlike)}")


def _incomplete(row):
    """Whether a surveyed genome has a gap worth showing."""
    return bool(row["incomplete"] or row["unnamed_genes"] or row["shortfall"]
                or row["non_canonical_starts"] or not row["transcripts"])


def _show_completeness(result, only_incomplete, json_output):
    rows = result["genomes"]
    if only_incomplete:
        rows = [row for row in rows if _incomplete(row)]
    if json_output:
        print_json(dict(result, genomes=rows))
        return
    totals = result["totals"]
    print(f"{result['collection']} ({result['panel']}): {totals['genomes']} genomes, "
          f"{result['transcript_count']} coding transcripts")
    print(f"fully annotated: {totals['genomes_fully_annotated']}")
    print(f"with an unnamed gene: {totals['genomes_with_unnamed_genes']}")
    print(f"with an incomplete transcript: {totals['genomes_with_incomplete_transcripts']}")
    print(f"with a spliced transcript: {totals['genomes_with_spliced_transcripts']}")
    if totals["genomes_behind_alternate"]:
        print(f"behind the alternate source: {totals['genomes_behind_alternate']} genomes, "
              f"{totals['coding_features_behind_alternate']} coding features")
    starts = result["initiations"]
    print("initiation: " + ", ".join(f"{name} {count}" for name, count in sorted(starts.items())))
    for transcript, recoding in sorted(result["recoded_transcripts"].items()):
        print(f"recoded: {transcript} ({recoding})")
    print()
    print_table(("VIRUS", "SOURCE", "CDSS", "COMPLETE", "INCOMPLETE", "UNNAMED", "SPLICED",
                 "ALTERNATE", "SHORT BY"),
                [(row["abbreviation"], row["annotation_source"], row["transcripts"],
                  row["complete"], row["incomplete"], row["unnamed_genes"], row["spliced"],
                  row["alternate_features"] if row["alternate_features"] is not None else "-",
                  row["shortfall"] or "") for row in rows],
                empty="No genomes to survey")


def _show_coverage(args, parser):
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
        print_json(dict(panel=PANEL, total=len(all_rows),
                        available=int(all_rows["available"].sum()), shown=len(rows),
                        entries=records(rows)))
    elif args.format == "tsv":
        rows[["abbreviation", "name", "status", "coding_models", "refseq",
              "genbank", "pave_id", "bundled_reference_id"]].to_csv(
            sys.stdout, sep="\t", index=False, na_rep="", lineterminator="\n")
    else:
        print(f"{PANEL}: {int(all_rows['available'].sum())}/{len(all_rows)} "
              f"with bundled source CDSs; showing {len(rows)}")
        print_table(("VIRUS", "STATUS", "CDSS", "GDC REFERENCE / SOURCE"),
                    [(row.abbreviation, "available" if row.available else "not bundled",
                      row.coding_models, row.bundled_reference_id or row.refseq or row.genbank)
                     for row in rows.itertuples(index=False)], empty="No matching references")
