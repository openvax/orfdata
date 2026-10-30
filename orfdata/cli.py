"""Query and export the offline ORF evidence reference."""

import argparse
import csv
import json
import re
import sys
from contextlib import redirect_stdout

from .catalog import EVIDENCE_METHODS, UORFCatalog, load_catalog, merge_catalogs
from .cli_output import print_result, print_table
from .cli_parser import ArgumentParser, UsageError, add_json, add_limit
from .evidence import _RENAMED, MS, SYNONYMS, RiboSeq, evidence_rule

COMMANDS = {
    "data": "Discover downloadable data, inspect local files, and fetch resources",
    "query": "Find ORFs and export catalogue records",
    "catalogs": "Manage optional ORF catalogues",
    "datasets": "Manage and query quantitative RNA/protein datasets",
    "proteome": "Build a reference proteome and export FASTA",
    "viral": "Inspect viral coverage, subsets, and installed reference layers",
    "setup": "Install a species' default catalogues and reference collections",
    "genomes": "Manage pinned genome annotations",
    "chains": "List and download assembly mapping files",
    "lift": "Map an ORF catalogue to another genome assembly",
    "help": "Show help for a command",
}


def _command_guide(prog):
    parser = ArgumentParser(
        prog=prog, description="ORFdata: discover reference data, query ORFs, and build proteomes.",
        epilog="Start here:\n"
               f"  {prog} data list                 See what can be downloaded or built\n"
               f"  {prog} data status               See what is available locally\n"
               f"  {prog} data info NAME            Inspect a resource before fetching it\n"
               f"  {prog} help COMMAND              Find subcommands and options\n\n"
               "Tables and summaries are the default; use --json for machine-readable output.\n"
               "'download' and 'downloads' are aliases for 'data'.\n"
               "Documentation: https://openvax.github.io/orfdata/")
    parser.add_argument("--version", action="store_true", help="Show the installed version")
    commands = parser.add_subparsers()
    for name, description in COMMANDS.items():
        commands.add_parser(name, help=description)
    return parser


def run(argv=None, prog=None, *, catalog_loader=load_catalog):
    """Discover commands, inspect resources, and explicitly query or export records."""
    argv = list(sys.argv[1:] if argv is None else argv)
    prog = prog or "orfdata"
    if argv[:1] == ["help"]:
        argv = argv[1:] + ["--help"]
    if not argv or argv in (["--help"], ["-h"]):
        _command_guide(prog).print_help()
        return
    if argv == ["--version"]:
        from .version import __version__

        print(f"orfdata {__version__}")
        return
    if argv and argv[0] == "catalogs":
        from .catalog_cli import run as catalog_run

        return catalog_run(argv[1:], prog=(prog or "orfdata") + " catalogs")
    if argv and argv[0] == "datasets":
        from .dataset_cli import run as dataset_run

        return dataset_run(argv[1:], prog=(prog or "orfdata") + " datasets")
    if argv and argv[0] == "proteome":
        from .proteome_cli import run as proteome_run

        return proteome_run(argv[1:], prog=(prog or "orfdata") + " proteome")
    if argv and argv[0] in ("data", "download", "downloads"):
        from .download_cli import run as download_run

        return download_run(argv[1:], prog=f"{prog or 'orfdata'} {argv[0]}")
    if argv and argv[0] == "viral":
        from .viral_cli import run as viral_run

        return viral_run(argv[1:], prog=(prog or "orfdata") + " viral")
    if argv and argv[0] == "setup":
        return _setup(argv[1:], prog=(prog or "orfdata") + " setup")
    if argv and argv[0] == "genomes":
        from .genomes_cli import run as genomes_run

        return genomes_run(argv[1:], prog=(prog or "orfdata") + " genomes")
    if argv and argv[0] in ("lift", "chains"):
        from . import liftover_cli

        command = liftover_cli.lift if argv[0] == "lift" else liftover_cli.chains
        return command(argv[1:], prog=f"{prog or 'orfdata'} {argv[0]}")
    if argv[0] == "query":
        argv = argv[1:]
    elif not argv[0].startswith("-"):
        _command_guide(prog).error(f"Unknown command {argv[0]!r}")
    parser = ArgumentParser(
        description="Query an installed ORF catalogue; preview records or export them explicitly.",
        prog=prog + " query",
        epilog="Examples:\n  orfdata query --gene TP53\n"
               "  orfdata query --orf-id c7riboseqorf56 --json\n"
               "  orfdata query --species mouse --format fasta\n\n"
               "Find optional sources with 'orfdata catalogs list'.")
    _add_query_arguments(parser)
    args = parser.parse_args(argv)
    return _query(args, parser, catalog_loader)


def _setup(argv, prog):
    from .install import InstallError, install_species

    parser = ArgumentParser(
        prog=prog, description="Install a species' ORF catalogues, chains and default reference collections, so its "
                               "proteomes work offline on every assembly the chains reach. Human setup downloads "
                               "the GDC viral collection too; '--collections none' skips it.")
    parser.add_argument("--species", required=True, help="e.g. human or mouse")
    parser.add_argument("--accept-terms", action="store_true",
                        help="Accept usage terms of restricted chain files, such as UCSC's non-commercial ones")
    parser.add_argument("--force", action="store_true", help="Install again even what is already installed")
    parser.add_argument("--catalogs", metavar="NAMES", type=_catalog_names,
                        help="'all' to add the species' optional catalogues, such as human's Ribo-seq ones, or "
                             "comma-separated catalogue names; the species' default set by default")
    parser.add_argument("--collections", metavar="NAMES", type=_catalog_names,
                        help="'none' to skip the species' reference collections, such as human's GDC viral "
                             "collection, 'all', or comma-separated names; the default set by default")
    parser.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--timeout", type=float, default=60,
                        help="Seconds a download may stall before it is retried, then fails (default 60)")
    add_json(parser)
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            result = install_species(args.species, accept_terms=args.accept_terms, force=args.force,
                                     cache_root=args.cache_root, timeout=args.timeout, catalogs=args.catalogs,
                                     collections=args.collections)
    except InstallError as error:  # What did install is still reported.
        print_result(error.result, json_output=args.json)
        parser.exit(1, f"{parser.prog}: error: {error}\n")
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    print_result(result, json_output=args.json)


def _catalog_names(value):
    """``--catalogs``: "all", or comma-separated catalogue names."""
    return value if value == "all" else [name for name in value.split(",") if name]


def _add_query_arguments(parser, *, catalog_path=True):
    if catalog_path:
        parser.add_argument("--catalog", action="append",
                            help="Use a compatible local JSON/JSON.gz catalogue; repeat to merge several")
        parser.add_argument("--species", help="Use this species' catalogue, e.g. human or mouse (default human)")
    else:
        parser.set_defaults(catalog=None, species=None)
    parser.add_argument("--gene", help="Exact source gene symbol or identifier")
    parser.add_argument("--transcript", help="Exact source transcript ID, including its version when supplied")
    parser.add_argument("--orf-id", help="Exact catalogue ORF identifier")
    parser.add_argument("--evidence", type=_evidence_name(), metavar="NAME",
                        help="Evidence an ORF needs: ms, ribosome (or riboseq), proteomics, immunopeptidomics, "
                             "ms_unclassified, or all (the default)")
    _add_rule_arguments(parser)
    parser.add_argument("--biotype", help="Source ORF class, as each source labels it: uORF, lncRNA, intORF, dORF "
                                          "(GENCODE, human); uORF, uoORF, iORF, dORF, doORF, lncORF (mouse Ribo-seq; "
                                          "RibORF also truncation, extension); AltProt, Isoform (OpenProt)")
    parser.add_argument("--min-protein-length", type=int, metavar="AA",
                        help="Select ORFs with at least this many amino acids; source catalogues stay unfiltered")
    parser.add_argument("--region", help="1-based inclusive interval: chr7:66205483-66205522")
    parser.add_argument(
        "--position", help="Map a 1-based variant position: chr7:66205483 (table or JSON)"
    )
    parser.add_argument(
        "--include-initiation-context",
        action="store_true",
        help="Include spliced start flanks in locus lookup",
    )
    parser.add_argument("--ref", help="Forward-genomic SNV REF; requires --alt and --position")
    parser.add_argument("--alt", help="Forward-genomic SNV ALT; compare initiation contexts")
    parser.add_argument(
        "--genome-build", help="Required for genomic queries; selects or must match the catalogue's assembly"
    )
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--format", choices=("table", "json", "tsv", "fasta", "gtf"), default="table",
                        help="Output format (default table; exports contain all matching records)")
    output.add_argument("--json", dest="format", action="store_const", const="json",
                        help="Write complete machine-readable JSON")
    add_limit(parser)


def _evidence_name(*extra):
    """An argparse type for evidence names: those evidence_rule takes, then ``extra``.

    Names from before 0.28.0 are accepted but not listed.
    """
    names = EVIDENCE_METHODS + tuple(SYNONYMS) + ("all",) + extra

    def name(value):
        if value not in names and value not in _RENAMED:
            raise argparse.ArgumentTypeError(f"invalid choice: {value!r} (choose from {', '.join(names)})")
        return value

    return name


def _add_rule_arguments(parser):
    """Ribo-seq and MS requirements, shared by ORF queries and proteome builds."""
    parser.add_argument("--min-ribo-seq-studies", type=int, metavar="N",
                        help="Require Ribo-seq in at least N studies the source names")
    parser.add_argument("--min-ms-peptides", type=int, metavar="N",
                        help="Require at least N distinct supporting MS peptides")
    parser.add_argument("--ms-methods", metavar="KINDS",
                        help="Count only these MS peptide kinds, comma-separated: proteomics (digested cellular "
                             "proteins), immunopeptidomics (peptides eluted from MHC molecules), ms_unclassified")
    parser.add_argument("--either", action="store_true",
                        help="Accept ORFs meeting the Ribo-seq or the MS requirement, instead of both")
    parser.add_argument("--require-resolved", action="store_true",
                        help="Require peptide records or named ribosome studies, not aggregate counts")
    parser.add_argument("--require-unique", action="store_true",
                        help="Require a peptide with one protein mapping in its MS source")
    parser.add_argument("--require-reviewed", action="store_true",
                        help="Require convincing reviewed spectra (benchmark scores 4-5)")


def _rule(args, method, flag, *, default=None):
    """The evidence rule from the Ribo-seq/MS requirements, or from ``method`` (set with ``flag``).

    Without requirements, ``method`` (or ``default``) and the --require-* options
    give the rule; ``"all"`` selects every ORF.
    """
    ribo = args.min_ribo_seq_studies is not None
    ms = args.min_ms_peptides is not None or args.ms_methods is not None
    if args.either and not (ribo and ms):
        raise UsageError("--either needs both --min-ribo-seq-studies and an MS requirement")
    if not (ribo or ms):
        method = default if method is None else method
        if method in (None, "all") and (args.require_unique or args.require_reviewed or args.require_resolved):
            raise UsageError(f"--require-unique, --require-reviewed and --require-resolved need {flag} "
                             "or a Ribo-seq/MS requirement")
        try:
            return evidence_rule(None if method == "all" else method, require_unique=args.require_unique,
                                 require_reviewed=args.require_reviewed, require_resolved=args.require_resolved)
        except ValueError as error:  # The evidence name and --require-* options don't combine.
            raise UsageError(str(error)) from error
    if method is not None:
        raise UsageError(f"Use {flag} or the Ribo-seq/MS requirements, not both")
    if not ms and (args.require_unique or args.require_reviewed):
        raise UsageError("--require-unique and --require-reviewed need an MS requirement")
    try:  # Values the rules reject, such as an unknown MS kind, are command-line mistakes.
        rules = [RiboSeq(args.min_ribo_seq_studies, resolved=args.require_resolved)] if ribo else []
        if ms:
            rules.append(MS(1 if args.min_ms_peptides is None else args.min_ms_peptides,
                            unique=args.require_unique, reviewed=args.require_reviewed,
                            resolved=args.require_resolved,
                            methods=args.ms_methods.split(",") if args.ms_methods is not None else None))
    except ValueError as error:
        raise UsageError(str(error)) from error
    if len(rules) == 1:
        return rules[0]
    return rules[0] | rules[1] if args.either else rules[0] & rules[1]


def _query(args, parser, catalog_loader):
    if args.region and args.position:
        parser.error("Choose --region or --position")
    if args.position and args.format not in ("table", "json"):
        parser.error("--position requires --format json")
    if (args.region or args.position) and not args.genome_build:
        parser.error("--region and --position require --genome-build")
    if (args.ref is not None or args.alt is not None) and (
        not args.ref or not args.alt or not args.position
    ):
        parser.error("SNV annotation requires --position, --ref and --alt")
    if args.ref and (
        args.ref not in ("A", "C", "G", "T")
        or args.alt not in ("A", "C", "G", "T")
        or args.ref == args.alt
    ):
        parser.error("Provide distinct single-base forward-genomic REF and ALT")
    try:
        paths = args.catalog or [None]
        if args.species or (args.genome_build and not args.catalog and catalog_loader is load_catalog):
            # The assembly picks the catalogue too, not only checks it.
            catalogs = [catalog_loader(p, species=args.species, genome_build=args.genome_build) for p in paths]
        else:
            catalogs = [catalog_loader(p) for p in paths]
        catalog = catalogs[0] if len(catalogs) == 1 else merge_catalogs(catalogs)
        options = {
            "gene": args.gene,
            "transcript_id": args.transcript,
            "evidence": _rule(args, args.evidence, "--evidence"),
            "biotype": args.biotype,
            "min_protein_length": args.min_protein_length,
            "genome_build": args.genome_build,
            "include_initiation_context": args.include_initiation_context or bool(args.ref),
        }
        locus = args.region or args.position
        if locus:
            pattern = r"([^:]+):([0-9]+)-([0-9]+)" if args.region else r"([^:]+):([0-9]+)"
            match = re.fullmatch(pattern, locus)
            if match is None:
                raise UsageError(f"Invalid locus {locus!r}; use CONTIG:START-END for --region, CONTIG:POSITION "
                                 "for --position")
            start = int(match.group(2))
            options.update(
                contig=match.group(1),
                start=start,
                end=int(match.group(3)) if args.region else start,
            )
        records = catalog.query(**options)
        if args.orf_id:
            try:
                orf_id = catalog.get(args.orf_id).orf_id  # Any source's ID; a misspelled one is an error.
            except KeyError:
                raise ValueError(f"No ORF {args.orf_id!r} in {catalog.source}") from None
            records = tuple(r for r in records if r.orf_id == orf_id)
        initiation_annotations = None
        if args.position:
            initiation_annotations = []
            for r in records:
                annotation = r.map_initiation_position(
                    options["start"], genome_build=args.genome_build
                )
                if args.ref and annotation is not None:
                    annotation = r.annotate_initiation_snv(
                        options["start"], args.ref, args.alt, genome_build=args.genome_build
                    )
                initiation_annotations.append(annotation)
            if args.ref and not any(a is not None for a in initiation_annotations):
                raise ValueError(
                    "No mapped initiation context covers the SNV; coding consequences require full transcript reconstruction"
                )
    except (ValueError, OSError, ImportError, RuntimeError) as exc:
        parser.fail(exc)
    if args.format == "table":
        print(f"{catalog.source}: {len(records)} matching ORFs ({catalog.genome_build})")
        print_table(("ORF", "GENE", "TRANSCRIPT", "CLASS", "AA"),
                    [(r.orf_id, r.gene_name or r.gene_id, r.transcript_id, r.biotype,
                      len(r.protein_sequence)) for r in records], limit=args.limit, empty="No matching ORFs")
        if args.position:
            print("\nPosition annotations")
            positions = [r.map_genomic_position(options["start"], genome_build=args.genome_build) or {}
                         for r in records]
            print_table(("ORF", "ROLE", "AA POSITION", "REF AA", "REFERENCE MS"),
                        [(r.orf_id, a.get("role"), a.get("amino_acid_position"),
                          a.get("reference_amino_acid"), a.get("reference_ms_support"))
                         for r, a in zip(records, positions, strict=True)], limit=args.limit)
            if any(a is not None for a in initiation_annotations):
                print("\nInitiation context (sequence assessment; translation change is not measured)")
                print_table(("ORF", "CONTEXT POSITION", "REFERENCE", "ALTERNATE"),
                            [(a["orf_id"], a.get("context_position"),
                              (a.get("initiation_assessment") or {}).get("display"),
                              ((a.get("alternate_context") or {}).get("assessment") or {}).get("display"))
                             for a in initiation_annotations if a is not None], limit=args.limit)
        print("\nUse --json for evidence and provenance; --format fasta, tsv, or gtf for export.")
    elif args.format == "json":
        output = {
            "schema": "isovar.uorf_query.v1"
            if isinstance(catalog, UORFCatalog)
            else "orfdata.query.v1",
            "metadata": catalog.metadata,
            "records": [r.to_dict() for r in records],
        }
        if args.position:
            output["position_annotations"] = [
                r.map_genomic_position(options["start"], genome_build=args.genome_build)
                for r in records
            ]
            output["initiation_annotations"] = initiation_annotations
        json.dump(output, sys.stdout, indent=2)
        sys.stdout.write("\n")
    elif args.format == "gtf":
        from .gtf import orf_gtf_lines

        # One transcript per ORF on its host gene: concatenate with the assembly's
        # reference GTF to annotate variants in ncORFs with PyEnsembl-based tools.
        sys.stdout.writelines(orf_gtf_lines(records))
    elif args.format == "fasta":
        for r in records:
            reference = "" if r.reference_context is None else " reference=" + r.reference_context.key
            scope = "" if isinstance(catalog, UORFCatalog) else f" taxon={r.taxon_id}{reference}"
            sys.stdout.write(
                f">{r.orf_id} gene={r.gene_name} transcript={r.transcript_id} assembly={r.genome_build} "
                f"biotype={r.biotype}{scope} kind=extra_orf source={catalog.source} gene_id={r.gene_id}\n"
                f"{r.protein_sequence}\n"
            )
    else:
        fields = (
            "orf_id",
            "gene_id",
            "gene_name",
            "transcript_id",
            "genome_build",
            "taxon_id",
            "species",
            "translation_table",
            "source",
            "contig",
            "strand",
            "biotype",
            "protein_sequence",
            "coding_blocks_json",
            "ribosome_studies_json",
            "peptides_json",
            "initiation_context_json",
        )
        scoped = not isinstance(catalog, UORFCatalog)
        provenance = ("members", "liftover", "lifted_from", "reference_mapping")
        if scoped:
            fields += ("start_codon", "reference_context_json", "summaries_json") + tuple(
                f + "_json" for f in provenance)
        elif any(r.summaries for r in records):
            # The released upstream-only columns have no place for aggregate counts.
            parser.error("The legacy upstream-ORF TSV cannot export aggregate evidence; use --format json")
        elif any(r.members or r.lifted_from for r in records):
            parser.error("The legacy upstream-ORF TSV cannot export merged or lifted records; use --format json")
        elif any(getattr(r, "start_codon", None) for r in records):
            parser.error("The legacy upstream-ORF TSV cannot export declared start codons; use --format json")
        writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in records:
            row = r.to_dict()
            for field in ("coding_blocks", "ribosome_studies", "peptides", "initiation_context"):
                row[field + "_json"] = json.dumps(row.pop(field), separators=(",", ":"))
            if scoped:
                row["reference_context_json"] = json.dumps(row.pop("reference_context", None), separators=(",", ":"))
                # Full native records keep source, metric definition, None/0 and shared keys.
                row["summaries_json"] = json.dumps(row.pop("summaries", []), separators=(",", ":"))
                for field in provenance:  # Present only on merged or lifted records.
                    row[field + "_json"] = json.dumps(row.pop(field, None), separators=(",", ":"))
            writer.writerow(row)


if __name__ == "__main__":
    run()
