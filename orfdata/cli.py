# PYTHON_ARGCOMPLETE_OK
"""The ``orfdata`` command: its command guide, dispatch to each command, and the ``query`` and ``setup`` commands.

Importing this module is cheap, because the shell's Tab key runs it up to
argcomplete on every press: commands import pandas, DataCache and the data
modules only when they run.
"""

import csv
import json
import os
import re
import sys
from collections.abc import Callable
from contextlib import redirect_stdout
from dataclasses import dataclass
from importlib import import_module

import argcomplete

from .cli_parser import ArgumentParser, UsageError
from .command_parsers import (
    catalogs_parser,
    chains_parser,
    data_parser,
    datasets_parser,
    genomes_parser,
    lift_parser,
    proteome_parser,
    query_parser,
    setup_parser,
    viral_parser,
)
from .evidence import MS, RiboSeq, evidence_rule


@dataclass(frozen=True)
class Command:
    """One top-level ``orfdata`` command: its line in the command guide, the function that runs it, and its parser.

    Attributes
    ----------
    description : str
        What the command guide says the command does.
    function : str or None
        ``"module:name"`` of the function that runs the command, called as
        ``function(argv, prog=...)`` with the arguments after the command name.
        It is imported only when the command runs, so the guide and the other
        commands never load it. None for ``help``, which :func:`run` answers
        itself.
    parser : callable or None
        Builds the parser ``function`` parses its arguments with, called as
        ``parser(prog=...)``; tab completion reads it. None for ``help``.
    """

    description: str
    function: str | None = None
    parser: Callable[..., ArgumentParser] | None = None

    def load(self):
        """Import and return the function that runs this command.

        Returns
        -------
        callable
            Takes the command's own arguments and ``prog``.

        Raises
        ------
        ValueError
            If the command has no function of its own, as for ``help``.
        """
        if self.function is None:
            raise ValueError("This command is answered by orfdata.cli.run itself")
        module, name = self.function.split(":")
        return getattr(import_module(module), name)


# Every top-level command, in the command guide's order. Adding one here lists it and runs it.
COMMAND_TABLE = {
    "data": Command("Discover downloadable data, inspect local files, and fetch resources",
                    "orfdata.download_cli:run", data_parser),
    "query": Command("Find ORFs and export catalogue records", "orfdata.cli:query", query_parser),
    "catalogs": Command("Manage optional ORF catalogues", "orfdata.catalog_cli:run", catalogs_parser),
    "datasets": Command("Manage and query quantitative RNA/protein datasets", "orfdata.dataset_cli:run",
                        datasets_parser),
    "proteome": Command("Build a reference proteome and export FASTA", "orfdata.proteome_cli:run", proteome_parser),
    "viral": Command("Inspect viral coverage, subsets, and installed reference layers", "orfdata.viral_cli:run",
                     viral_parser),
    "setup": Command("Install a species' default catalogues and reference collections", "orfdata.cli:setup",
                     setup_parser),
    "genomes": Command("Manage pinned genome annotations", "orfdata.genomes_cli:run", genomes_parser),
    "chains": Command("List and download assembly mapping files", "orfdata.liftover_cli:chains", chains_parser),
    "lift": Command("Map an ORF catalogue to another genome assembly", "orfdata.liftover_cli:lift", lift_parser),
    "help": Command("Show help for a command"),
}
# Other names a command answers to; its help and errors name it as typed.
COMMAND_ALIASES = {"download": "data", "downloads": "data"}
# Each command's guide line by name, as before COMMAND_TABLE existed.
COMMANDS = {name: command.description for name, command in COMMAND_TABLE.items()}


def command_guide(prog="orfdata", *, parsers=None):
    """Return the parser whose help is the ``orfdata`` command guide.

    It lists every command in :data:`COMMAND_TABLE`, and parses nothing else:
    :func:`run` uses it to print the guide and to report an unknown command.

    Parameters
    ----------
    prog : str, optional
        The program name the guide shows, ``orfdata`` by default.
    parsers : dict, optional
        Parsers to put under the guide by command name or alias, in place of
        the empty ones it lists, as :func:`completion_parser` does so tab
        completion reaches a command's options. The guide's help is the same.

    Returns
    -------
    ArgumentParser
        Its ``print_help()`` prints the guide; its ``error()`` exits with
        status 2.
    """
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
    for name, command in COMMAND_TABLE.items():
        commands.add_parser(name, help=command.description)
    # Only tab completion reads this: run() answers 'help COMMAND' itself.
    commands.choices["help"].add_argument("command", nargs="?", metavar="COMMAND",
                                          choices=[name for name in COMMAND_TABLE if name != "help"])
    commands.choices.update(parsers or {})
    return parser


def completion_words(environ=None):
    """Return the words typed after the program name and before the cursor, from a tab-completion request.

    Parameters
    ----------
    environ : mapping, optional
        The environment argcomplete's shell hook sets, with ``COMP_LINE``,
        ``COMP_POINT`` and ``_ARGCOMPLETE``; :data:`os.environ` by default.

    Returns
    -------
    list of str
        The complete words, without the program name and without the word
        being completed.
    """
    environ = os.environ if environ is None else environ
    words = argcomplete.split_line(environ["COMP_LINE"], int(environ["COMP_POINT"]))[3]
    # _ARGCOMPLETE is where the program's own arguments start: 1 for 'orfdata', 3 for 'python -m orfdata.cli'.
    return words[int(environ["_ARGCOMPLETE"]):]


def completion_parser(words, prog="orfdata"):
    """Return the parser tab completion reads: the command guide, with the typed command's parser under it.

    Only the command named in ``words`` gets its parser built, so a Tab press
    builds one command's options at most, and imports none of the modules that
    run commands.

    Parameters
    ----------
    words : list of str
        The words already typed after the program name, as
        :func:`completion_words` returns them.
    prog : str, optional
        The program name, ``orfdata`` by default.

    Returns
    -------
    ArgumentParser
        Completes command names, and the options and subcommands of the
        command in ``words``, aliases included.
    """
    name = next((word for word in words if not word.startswith("-")), None)
    command = COMMAND_TABLE.get(COMMAND_ALIASES.get(name, name))
    if command is None or command.parser is None:
        return command_guide(prog)
    return command_guide(prog, parsers={name: command.parser(prog=f"{prog} {name}")})


def run(argv=None, prog=None, *, catalog_loader=None):
    """Run the ``orfdata`` command line: the command guide, ``--version``, or one command.

    ``help COMMAND`` shows that command's help. A command named in
    :data:`COMMAND_TABLE` or :data:`COMMAND_ALIASES` runs with the arguments
    after its name. Arguments that start with an option, such as
    ``--gene TP53``, are a ``query``, as they were before ``query`` was a
    command. When the shell's Tab key runs it, through argcomplete, it prints
    the completions and exits instead.

    Parameters
    ----------
    argv : list of str, optional
        The arguments after the program name; ``sys.argv[1:]`` by default.
    prog : str, optional
        The program name that help and errors show, ``orfdata`` by default.
    catalog_loader : callable, optional
        What ``query`` opens catalogues with, as in :func:`query`; other
        commands don't use it. None means :func:`~orfdata.load_catalog`.

    Returns
    -------
    object
        What the command's function returns; None for the guide and
        ``--version``.

    Raises
    ------
    SystemExit
        With status 2 for an unknown command, and however the command exits.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    prog = prog or "orfdata"
    if "_ARGCOMPLETE" in os.environ:  # A Tab press: argcomplete prints the completions and exits.
        argcomplete.autocomplete(completion_parser(completion_words(), prog))
    if argv[:1] == ["help"]:
        argv = argv[1:] + ["--help"]
    if not argv or argv in (["--help"], ["-h"]):
        command_guide(prog).print_help()
        return
    if argv == ["--version"]:
        from .version import __version__

        print(f"orfdata {__version__}")
        return
    name, arguments = argv[0], argv[1:]
    if name.startswith("-"):
        return query(argv, prog + " query", catalog_loader=catalog_loader)
    if name == "query":  # The one command that takes a catalogue loader.
        return query(arguments, prog + " query", catalog_loader=catalog_loader)
    command = COMMAND_TABLE.get(COMMAND_ALIASES.get(name, name))
    if command is None or command.function is None:
        command_guide(prog).error(f"Unknown command {name!r}")  # Exits with status 2.
    return command.load()(arguments, prog=f"{prog} {name}")


def query(argv=None, prog="orfdata query", *, catalog_loader=None):
    """Run the ``query`` command on its own, under ``prog``.

    Another tool can expose ORF queries as one of its commands with its own catalogue loader. Its help,
    errors and examples name ``prog``, and exit statuses follow :func:`run`. Any ``catalog_loader`` other
    than :func:`~orfdata.load_catalog`, which None stands for, reads one catalogue, so ``--species`` is not
    offered.
    """
    from .catalog import load_catalog

    argv = list(sys.argv[1:] if argv is None else argv)
    catalog_loader = load_catalog if catalog_loader is None else catalog_loader
    parser = query_parser(prog, species=catalog_loader is load_catalog)
    args = parser.parse_args(argv)
    return _query(args, parser, catalog_loader)


def setup(argv=None, prog="orfdata setup"):
    """Run the ``setup`` command: install a species' catalogues, chains and default reference collections.

    Parameters
    ----------
    argv : list of str, optional
        The command's arguments, such as ``["--species", "human"]``;
        ``sys.argv[1:]`` by default.
    prog : str, optional
        The program name that help and errors show.

    Raises
    ------
    SystemExit
        With status 2 when the arguments don't parse, and 1 when the species
        is unknown or anything failed to install; what did install is still
        reported first.
    """
    from .cli_output import print_result
    from .install import InstallError, install_species

    parser = setup_parser(prog)
    args = parser.parse_args(argv)
    try:
        with redirect_stdout(sys.stderr):
            result = install_species(args.species, accept_terms=args.accept_terms, force=args.force,
                                     cache_root=args.cache_root, timeout=args.timeout, catalogs=args.catalogs,
                                     collections=args.collections, resume=args.resume)
    except InstallError as error:  # What did install is still reported.
        print_result(error.result, json_output=args.json)
        parser.exit(1, f"{parser.prog}: error: {error}\n")
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    print_result(result, json_output=args.json)


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
    from .catalog import UORFCatalog, load_catalog, merge_catalogs
    from .cli_output import print_table

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
        # The whole command line is checked before any catalogue is read, so a mistake in it exits 2
        # whatever is installed.
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
            end = int(match.group(3)) if args.region else start
            if start < 1 or end < start:
                raise UsageError(f"Invalid locus {locus!r}; positions are 1-based and END must be at least START")
            options.update(contig=match.group(1), start=start, end=end)
        if args.min_protein_length is not None and args.min_protein_length < 1:
            raise UsageError("--min-protein-length must be a positive number of amino acids")
        paths = args.catalog or [None]
        if args.species or (args.genome_build and not args.catalog and catalog_loader is load_catalog):
            # The assembly picks the catalogue too, not only checks it.
            catalogs = [catalog_loader(p, species=args.species, genome_build=args.genome_build) for p in paths]
        else:
            catalogs = [catalog_loader(p) for p in paths]
        catalog = catalogs[0] if len(catalogs) == 1 else merge_catalogs(catalogs)
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
