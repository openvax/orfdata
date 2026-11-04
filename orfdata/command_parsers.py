"""The parser of every ``orfdata`` command, built without importing ORFdata's data modules.

Each command's ``run`` parses its arguments with the parser built here, and
tab completion builds the typed command's parser from here on every Tab press,
so this module imports only :mod:`argparse` and other cheap modules. The
functions that run the commands import pandas, DataCache and the rest.
"""

import argparse

from .cli_parser import (
    ArgumentParser,
    add_cache_root,
    add_json,
    add_limit,
    add_resume,
    add_timeout,
    share_cache_root,
)
from .evidence import _RENAMED, EVIDENCE_METHODS, SYNONYMS


def data_parser(prog="orfdata download"):
    """Return the parser of ``orfdata data``, which discovers and fetches every kind of resource.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``list``, ``status``, ``info``, ``path``,
        ``fetch``, ``viral``, ``catalog``, ``dataset``, ``reference``,
        ``genome`` and ``chains``.
    """
    parser = ArgumentParser(prog=prog, description="Discover and acquire ORFdata resources through one command.",
                            epilog="Start here:\n  orfdata data list --species dog\n"
                                   "  orfdata data status\n  orfdata data info NAME\n"
                                   "  orfdata data fetch NAME\n\n"
                                   "List, status, info, and path never download data.")
    commands = parser.add_subparsers(dest="command", required=True)

    listing = commands.add_parser("list", help="List every bundled, downloadable or buildable resource")
    add_discovery(listing)
    status = commands.add_parser("status", help="Show resources with local files, including damaged and unreadable ones")
    add_discovery(status, status=True)

    for verb, description in (("info", "Describe a resource's source pins and download/build route"),
                              ("path", "Print the installed path without downloading"),
                              ("fetch", "Install named resources, downloading or building as needed")):
        command = commands.add_parser(verb, help=description)
        command.add_argument("names", metavar="NAME", nargs="+" if verb == "fetch" else None)
        command.add_argument("--kind", choices=("viral", "catalog", "dataset", "reference", "genome", "chain", "records"))
        command.add_argument("--version", help="Exact resource version; defaults to its registered version")
        add_json(command, help="Write machine-readable JSON")
        if verb == "fetch":
            add_force_timeout(command)
            command.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH")
            command.add_argument("--quiet", action="store_true")
            command.add_argument("--accept-terms", action="store_true")
            add_resume(command)
        else:
            command.add_argument("--verify", action="store_true", help="Hash and validate installed files")

    viral = commands.add_parser(
        "viral", help="Download the current human-targeted viral collection and paired FASTAs")
    viral.add_argument("--species", default="human",
                       help="Target analysis species (currently human; future collections may differ)")
    viral.add_argument(
        "--collection", default="gdc-viruses",
        help="Reference collection (default gdc-viruses: 200 GDC references for human analysis)")
    viral.add_argument("--context", action="append",
                       help="Build an index for this context; repeat for a union")
    viral.add_argument("--exclude-context", action="append",
                       help="Exclude this context from the derived index; repeat as needed")
    viral.add_argument("--context-match", choices=("any", "all"), default="any",
                       help="Match any selected context (default) or require all")
    viral.add_argument("--subset", action="append",
                       help="Build an index for this named viral subset; repeat for a union")
    viral.add_argument("--inventory-molecules", dest="prefer_reference_molecule",
                       action="store_false",
                       help="Carry the molecule each inventory row names first, even where the "
                            "reference the collection follows was built from the other one")
    viral.add_argument("--no-annotation-fill", dest="annotate_from_alternate",
                       action="store_false",
                       help="Keep strictly to each entry's own record, leaving one whose source "
                            "annotates nothing without transcripts")
    viral.add_argument("--alternate-annotations", dest="prefer_alternate_annotation",
                       action="store_true",
                       help="Take the alternate source's coding features wherever it publishes "
                            "the same molecule and annotates at least as much of it, for "
                            "uniform gene models in place of each submission's own")
    add_force_timeout(viral)
    add_json(viral, help="Write machine-readable JSON")

    catalog = commands.add_parser("catalog", help="Download or build a registered ORF catalogue")
    catalog.add_argument("name")
    catalog.add_argument("--version")
    catalog.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH")
    add_force_timeout(catalog)
    add_json(catalog, help="Write machine-readable JSON")

    dataset = commands.add_parser("dataset", help="Download or build a registered quantitative dataset")
    dataset.add_argument("name")
    dataset.add_argument("--version")
    dataset.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH")
    dataset.add_argument("--quiet", action="store_true")
    add_force_timeout(dataset)
    add_json(dataset, help="Write machine-readable JSON")

    reference = commands.add_parser("reference", help="Download a dataset's native protein search FASTA")
    reference.add_argument("dataset")
    reference.add_argument("--version")
    reference.add_argument("--reference", help="Exact native reference context key")
    reference.add_argument("--asset", help="Declared source asset ID")
    reference.add_argument("--source-path", metavar="FILE")
    reference.add_argument("--quiet", action="store_true")
    add_force_timeout(reference)
    add_json(reference, help="Write machine-readable JSON")

    genome = commands.add_parser("genome", help="Download and index a pinned PyEnsembl annotation")
    genome.add_argument("name")
    genome.add_argument("--source-path", action="append", default=[], metavar="ROLE=PATH")
    add_force_timeout(genome)
    add_json(genome, help="Write machine-readable JSON")

    chains = commands.add_parser("chains", help="Download pinned assembly chain files")
    chains.add_argument("source_build", metavar="FROM", nargs="?")
    chains.add_argument("target_build", metavar="TO", nargs="?")
    chains.add_argument("--species")
    chains.add_argument("--accept-terms", action="store_true")
    chains.add_argument("--source-path", metavar="FILE")
    add_force_timeout(chains)
    add_json(chains, help="Write machine-readable JSON")

    share_cache_root(parser)
    return parser


def add_force_timeout(parser):
    """Add ``--force``, which replaces or repairs an installation, and ``--timeout``.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser of a command that installs something.
    """
    parser.add_argument("--force", action="store_true", help="Replace or repair an existing installation")
    add_timeout(parser)


def add_discovery(parser, *, status=False):
    """Add the filters ``data list`` and ``data status`` share, such as ``--kind`` and ``--species``.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The ``list`` or ``status`` parser.
    status : bool, optional
        True for ``status``, which shows anything with local files; False
        for ``list``, which also takes ``--installed`` or ``--missing``.
    """
    parser.add_argument(
        "--kind", choices=("viral", "catalog", "dataset", "reference", "genome", "chain", "records"))
    parser.add_argument("--species", help="Only resources for this species, e.g. human or mouse")
    parser.add_argument("--context", action="append",
                        help="Only resources in this context; repeat for a union")
    parser.add_argument("--context-match", choices=("any", "all"), default="any",
                        help="Match any selected context (default) or require all")
    add_json(parser, help="Write machine-readable JSON")
    add_limit(parser)
    parser.add_argument("--verify", action="store_true", help="Hash and validate installed files")
    if status:
        parser.set_defaults(installed=None, local=True)  # status: anything with local files
    else:
        selection = parser.add_mutually_exclusive_group()
        selection.add_argument("--installed", dest="installed", action="store_const", const=True,
                               help="Show only locally installed or bundled resources")
        selection.add_argument("--missing", dest="installed", action="store_const", const=False,
                               help="Show only missing, corrupt or inaccessible resources")
        parser.set_defaults(installed=None, local=False)


def query_parser(prog="orfdata query", *, species=True):
    """Return the parser of ``orfdata query``, which finds ORFs and exports catalogue records.

    Parameters
    ----------
    prog : str, optional
        The program name that help, errors and examples show.
    species : bool, optional
        Whether ``--species`` picks the catalogue, as it does with ORFdata's
        own catalogue loader. Another tool's loader reads one catalogue, so it
        passes False and the option isn't offered.

    Returns
    -------
    ArgumentParser
    """
    epilog = (f"Examples:\n  {prog} --gene TP53\n"
              f"  {prog} --orf-id c7riboseqorf56 --json\n")
    if species:
        epilog += (f"  {prog} --species mouse --format fasta\n\n"
                   "Find optional sources with 'orfdata catalogs list'.")
    parser = ArgumentParser(
        description="Query an installed ORF catalogue; preview records or export them explicitly.",
        prog=prog, epilog=epilog)
    add_query_arguments(parser, species=species)
    return parser


def add_query_arguments(parser, *, catalog_path=True, species=True):
    """Add the options of an ORF query: what to select, the evidence it needs, and the output format.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser of ``orfdata query`` or ``orfdata catalogs query``.
    catalog_path : bool, optional
        Whether ``--catalog`` (and ``--species``) choose the catalogue. False
        for ``catalogs query``, which names its catalogue; both are then None.
    species : bool, optional
        Whether ``--species`` is offered; when False it is None.
    """
    if catalog_path:
        parser.add_argument("--catalog", action="append",
                            help="Use a compatible local JSON/JSON.gz catalogue; repeat to merge several")
    else:
        parser.set_defaults(catalog=None)
    if catalog_path and species:
        parser.add_argument("--species", help="Use this species' catalogue, e.g. human or mouse (default human)")
    else:
        parser.set_defaults(species=None)
    parser.add_argument("--gene", help="Exact source gene symbol or identifier")
    parser.add_argument("--transcript", help="Exact source transcript ID, including its version when supplied")
    parser.add_argument("--orf-id", help="Exact catalogue ORF identifier")
    parser.add_argument("--evidence", type=evidence_name(), metavar="NAME",
                        help="Evidence an ORF needs: ms, ribosome (or riboseq), proteomics, immunopeptidomics, "
                             "ms_unclassified, or all (the default)")
    add_rule_arguments(parser)
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


def evidence_name(*extra):
    """Return an argparse type for evidence names: those :func:`~orfdata.evidence.evidence_rule` takes, then ``extra``.

    Names from before 0.28.0 are accepted but not listed.

    Parameters
    ----------
    *extra : str
        Further names to accept, after ``"all"``.

    Returns
    -------
    callable
        Returns the name it is given; raises :class:`argparse.ArgumentTypeError`,
        listing the choices, for any other.
    """
    names = EVIDENCE_METHODS + tuple(SYNONYMS) + ("all",) + extra

    def name(value):
        if value not in names and value not in _RENAMED:
            raise argparse.ArgumentTypeError(f"invalid choice: {value!r} (choose from {', '.join(names)})")
        return value

    return name


def add_rule_arguments(parser):
    """Add the Ribo-seq and MS requirements an ORF query selects ORFs by.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser of an ORF query.
    """
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


def setup_parser(prog="orfdata setup"):
    """Return the parser of ``orfdata setup``, which installs a species' catalogues, chains and collections.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show.

    Returns
    -------
    ArgumentParser
    """
    parser = ArgumentParser(
        prog=prog, description="Install a species' ORF catalogues, chains and default reference collections, so its "
                               "proteomes work offline on every assembly the chains reach. Human setup downloads "
                               "the GDC viral collection too; '--collections none' skips it.")
    parser.add_argument("--species", required=True, help="e.g. human or mouse")
    parser.add_argument("--accept-terms", action="store_true",
                        help="Accept usage terms of restricted chain files, such as UCSC's non-commercial ones")
    parser.add_argument("--force", action="store_true", help="Install again even what is already installed")
    parser.add_argument("--catalogs", metavar="NAMES", type=catalog_names,
                        help="'all' to add the species' optional catalogues, such as human's Ribo-seq ones, or "
                             "comma-separated catalogue names; the species' default set by default")
    parser.add_argument("--collections", metavar="NAMES", type=catalog_names,
                        help="'none' to skip the species' reference collections, such as human's GDC viral "
                             "collection, 'all', or comma-separated names; the default set by default")
    add_cache_root(parser)
    add_timeout(parser)
    add_resume(parser)
    add_json(parser)
    return parser


def catalog_names(value):
    """Read ``setup``'s ``--catalogs`` or ``--collections``: "all", or comma-separated names.

    Parameters
    ----------
    value : str
        The option's value.

    Returns
    -------
    str or list of str
        ``"all"`` as given; otherwise the nonempty names, in order. ``"none"``
        comes back as ``["none"]``, which ``setup`` reads as no reference
        collections when ``--collections`` gives it.
    """
    return value if value == "all" else [name for name in value.split(",") if name]


def catalogs_parser(prog=None):
    """Return the parser of ``orfdata catalogs``, which discovers, installs and queries ORF catalogues.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``list``, ``status``, ``info``, ``inspect``,
        ``download``, ``build`` and ``query``.
    """
    parser = ArgumentParser(description="Explicit optional catalogue discovery, acquisition and offline queries.",
                            prog=prog,
                            epilog="Start with 'orfdata catalogs list', then 'orfdata catalogs info NAME'.\n"
                                   "Use 'orfdata data fetch NAME' to select the published download or source build.")
    parser.add_argument("--registry", help="Local JSON catalogue registry definitions")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Discover release manifests without downloading")
    listing.add_argument("--species", help="Only this species' catalogues, e.g. human or mouse")
    listing.add_argument("--taxon-id", type=int)
    add_json(listing)
    add_limit(listing)
    status = commands.add_parser("status", help="Show bundled and installed catalogue status for every species")
    add_json(status)
    add_limit(status)
    for command, description in {
        "info": "Describe a catalogue's references, sources, and evidence coverage",
        "inspect": "Check the local artifact without downloading",
        "download": "Download and verify a published catalogue",
        "build": "Build a catalogue from pinned source files",
        "query": "Find ORFs in an installed catalogue and export records",
    }.items():
        subparser = commands.add_parser(command, help=description, description=description)
        subparser.add_argument("name", help="Catalogue name, as 'orfdata catalogs list' shows it")
        subparser.add_argument("--version", help="Release to use instead of the catalogue's default")
        if command in ("inspect", "query"):
            subparser.add_argument("--built", action="store_true", help="Select a local source build explicitly")
        if command in ("download", "build"):
            subparser.add_argument("--force", action="store_true", help="Replace an existing installation")
            subparser.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH",
                                   help="Use this local copy of a pinned source asset instead of downloading it")
            add_timeout(subparser)
        if command == "query":
            add_query_arguments(subparser, catalog_path=False)
        else:
            add_json(subparser)
    share_cache_root(parser)
    return parser


def datasets_parser(prog=None):
    """Return the parser of ``orfdata datasets``, which discovers, installs and queries quantitative datasets.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``list`` and the commands on one dataset, such
        as ``info``, ``download``, ``query`` and ``pairs``.
    """
    parser = ArgumentParser(description="Discover, download, inspect and query quantitative datasets.", prog=prog,
                            epilog="Start with 'orfdata datasets list', then 'orfdata datasets info NAME'.\n"
                                   "Queries preserve native units and reference identities.")
    parser.add_argument("--registry", help="Local registry JSON (otherwise use packaged releases)")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="List published releases and assay coverage")
    listing.add_argument("--species")
    listing.add_argument("--taxon-id", type=int)
    listing.add_argument("--assay")
    add_json(listing)
    add_limit(listing)
    for command, description in {
        "info": "Describe a dataset's sources, references, and assay coverage",
        "download": "Download and verify a published dataset",
        "build": "Build a dataset from pinned source files",
        "inspect": "Check the local dataset without downloading",
        "query": "Preview or export measurements in their native units",
        "pairs": "Find measurements with an explicit supported pairing basis",
        "references": "List the dataset's native protein search references",
        "download-reference": "Download a pinned native protein search FASTA",
        "inspect-reference": "Check an installed search FASTA without downloading",
        "query-reference": "Check an exact accession in an installed search reference",
    }.items():
        subparser = commands.add_parser(command, help=description, description=description)
        subparser.add_argument("name", help="Dataset name from 'orfdata datasets list'")
        subparser.add_argument("--version", help="Exact dataset release (default: registered version)")
        if command in ("query", "pairs"):
            output = subparser.add_mutually_exclusive_group()
            output.add_argument("--format", choices=("table", "json", "jsonl"), default="table",
                                help="Output format (default table; exports contain every matching record)")
            output.add_argument("--json", dest="format", action="store_const", const="json",
                                help="Write a complete machine-readable JSON array")
        else:
            add_json(subparser)
        if command in ("query", "pairs", "references"):
            add_limit(subparser)
        if command in ("inspect", "query", "pairs"):
            subparser.add_argument("--built", action="store_true",
                                   help="Use a local source build explicitly, instead of a published artifact")
        if command in ("download", "download-reference"):
            subparser.add_argument("--force", action="store_true")
            add_timeout(subparser)
            subparser.add_argument("--quiet", action="store_true")
            subparser.add_argument("--source-path", help="Use a local pinned file (compressed when declared)")
        if command in ("download-reference", "inspect-reference", "query-reference"):
            subparser.add_argument("--reference", help="Exact native reference context key")
            subparser.add_argument("--asset", help="Declared source asset ID")
        if command == "query-reference":
            subparser.add_argument("--target", required=True, help="Exact native accession")
            subparser.add_argument("--kind", default="protein")
            subparser.add_argument("--namespace", help="Default: the declared source namespace")
        if command == "build":
            subparser.add_argument("--force", action="store_true")
            add_timeout(subparser)
            subparser.add_argument("--quiet", action="store_true")
            subparser.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH",
                                   help="Use a local pinned source; repeat for each input asset")
        if command == "query":
            subparser.add_argument("--reference", help="Exact ref:SHA256 context key")
            subparser.add_argument("--sample")
            subparser.add_argument("--assay-id")
            subparser.add_argument("--quantity")
            subparser.add_argument("--status")
            subparser.add_argument("--unit", help="Exact native unit, e.g. iBAQ, rIBAQ, TPM or counts")
            subparser.add_argument("--scale", choices=("linear", "log"))
            subparser.add_argument("--log-base", type=float, help="Exact native logarithm base; no transformation")
            subparser.add_argument("--target")
            subparser.add_argument("--kind")
            subparser.add_argument("--namespace")
            subparser.add_argument("--mapped", action="store_true",
                                   help="Return source observations with explicit reference mappings")
            subparser.add_argument("--mapping-status", action="append",
                                   help="Mapping status to include; default unique")
            subparser.add_argument("--require-sequence", action="store_true",
                                   help="Require recorded exact sequence validation for mapped queries")
        if command == "pairs":
            subparser.add_argument("--basis", action="append", required=True,
                                   help="Allowed pairing basis; repeat for multiple bases")
            subparser.add_argument("--temporal", action="append",
                                   help="Allowed temporal relation; default concurrent")
            subparser.add_argument("--reference")
            subparser.add_argument("--pair-id", action="append")
            subparser.add_argument("--relation-id", action="append")
            for side in ("left", "right"):
                for field in ("target", "kind", "namespace", "reference", "quantity", "status", "unit"):
                    subparser.add_argument(f"--{side}-{field}")
                subparser.add_argument(f"--{side}-scale", choices=("linear", "log"))
                subparser.add_argument(f"--{side}-log-base", type=float)
    share_cache_root(parser)
    return parser


def proteome_parser(prog=None):
    """Return the parser of ``orfdata proteome``, which builds a reference proteome and writes FASTA.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
    """
    parser = ArgumentParser(
        description="Write a proteome from selected host protein sources and optional viruses.", prog=prog,
        epilog="The PyEnsembl annotation must already be installed, e.g. "
               "'pyensembl install --release 115 --species human', or, for --genome, with "
               "'orfdata genomes install NAME'.")
    parser.add_argument("--species", help="PyEnsembl species, e.g. human, mouse or dog")
    parser.add_argument("--release", type=int, help="Ensembl release")
    parser.add_argument("--genome", metavar="NAME",
                        help="A pinned annotation installed with 'orfdata genomes install', e.g. t2t-chm13, "
                             "instead of --species and --release")
    parser.add_argument("--isoforms", choices=("all", "longest"), default="all",
                        help="Reference CDS isoforms: all (default), or the longest protein per gene")
    parser.add_argument("--nc-orf-sources", action="append", metavar="SOURCE[,SOURCE...]",
                        help="Add ncORFs from 'all', evidence presets such as 'ms' or 'ribosome', "
                             "dataset names, or local JSON/JSON.gz files. Repeat or separate by commas "
                             "to combine sources. No ncORFs by default")
    parser.add_argument("--no-reference-cds", action="store_true",
                        help="Omit proteins translated from the reference annotation's CDS features")
    parser.add_argument("--min-protein-length", type=int, metavar="AA",
                        help="Add only ncORFs with at least this many amino acids; source catalogues stay unfiltered")
    parser.add_argument("--viral", action="append", metavar="REFERENCE",
                        help="Add a bundled virus by type (e.g. HPV16) or exact reference ID; "
                             "use 'all' for every distinct bundled virus recorded for this host; repeat for several")
    parser.add_argument("--snapshot", help="Also save a checksummed snapshot that Proteome.load reads offline")
    parser.add_argument("--cache", action="store_true",
                        help="Reuse an identical earlier build from the OpenVax DataCache, or keep this one there")
    add_cache_root(parser, help="Cache directory holding catalogues, chains, pinned genomes and cached builds, "
                                "instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", help="FASTA path; standard output by default")
    parser.add_argument("--annotation", metavar="DIRECTORY",
                        help="Also write the ncORF models there as annotation.gtf with matching "
                             "transcripts.fasta and proteins.fasta, for PyEnsembl-based tools")
    parser.add_argument("--json", action="store_true", help="Write the complete build summary as JSON to stderr")
    return parser


def viral_parser(prog="orfdata viral"):
    """Return the parser of ``orfdata viral``, which describes the viral collection published and installed.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``coverage``, ``contexts``, ``subsets``,
        ``summary``, ``completeness`` and ``inventory download``.
    """
    parser = ArgumentParser(prog=prog, description=(
        "Inspect GDC viral coverage, contexts, subsets and installed reference layers.\n\n"
        "Acquisition lives in 'orfdata download viral'; this command describes what is\n"
        "published and what is installed.\n"))
    commands = parser.add_subparsers(dest="command", required=True)
    coverage = commands.add_parser("coverage", help="List all GDC human-screen rows and coding-reference status")
    selection = coverage.add_mutually_exclusive_group()
    selection.add_argument("--available", action="store_true", help="Show rows with bundled source CDSs")
    selection.add_argument("--missing", action="store_true", help="Show rows not yet bundled in ORFdata")
    coverage.add_argument("--virus", help="Exact GDC inventory abbreviation, e.g. HPV18")
    coverage.add_argument("--search", help="Case-insensitive name, label or source-accession search")
    add_viral_selection(coverage)
    output = coverage.add_mutually_exclusive_group()
    output.add_argument("--format", choices=("table", "json", "tsv"), default="table")
    output.add_argument("--json", dest="format", action="store_const", const="json",
                        help="Write machine-readable JSON")
    add_limit(coverage)

    contexts = commands.add_parser("contexts", help="List selectable viral-reference contexts")
    contexts.add_argument("--json", action="store_true", help="Write machine-readable JSON")
    add_limit(contexts)

    subsets = commands.add_parser("subsets", help="List named viral-reference subsets")
    subsets.add_argument("--json", action="store_true", help="Write machine-readable JSON")
    add_limit(subsets)

    summary = commands.add_parser(
        "summary", help="Summarize the installed collection's genome, transcriptome and proteome")
    add_viral_selection(summary)
    summary.add_argument("--subset", action="append",
                         help="Count only this named viral subset; repeat for a union")
    summary.add_argument("--species", default="human", help="Target analysis species (currently human)")
    add_cache_root(summary)  # share_cache_root also accepts it before the subcommand.
    summary.add_argument("--json", action="store_true", help="Write machine-readable JSON")

    completeness = commands.add_parser(
        "completeness", help="Report how completely each installed genome is annotated")
    add_viral_selection(completeness)
    completeness.add_argument("--subset", action="append",
                              help="Survey only this named viral subset; repeat for a union")
    completeness.add_argument("--species", default="human",
                              help="Target analysis species (currently human)")
    add_cache_root(completeness)
    completeness.add_argument("--incomplete", action="store_true",
                              help="Show only genomes with a gap worth looking at")
    completeness.add_argument("--json", action="store_true", help="Write machine-readable JSON")
    add_limit(completeness)

    inventory = commands.add_parser("inventory", help="Manage the published GDC inventory metadata")
    inventory_commands = inventory.add_subparsers(dest="inventory_command", required=True)
    download = inventory_commands.add_parser("download", help="Download the pinned inventory via DataCache")
    download.add_argument("--force", action="store_true", help="Replace the cached copy from the pinned source")
    add_timeout(download)

    add_json(download)
    share_cache_root(parser)
    return parser


def add_viral_selection(parser):
    """Add ``--context``, ``--exclude-context`` and ``--context-match``, which select viral references.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser of a ``viral`` subcommand.
    """
    parser.add_argument("--context", action="append",
                        help="Select references in this context; repeat for a union")
    parser.add_argument("--exclude-context", action="append",
                        help="Exclude references in this context; repeat as needed")
    parser.add_argument("--context-match", choices=("any", "all"), default="any",
                        help="Match any selected context (default) or require all")


def genomes_parser(prog=None):
    """Return the parser of ``orfdata genomes``, which lists and installs pinned genome annotations.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``list`` and ``install``.
    """
    parser = ArgumentParser(description="List and install the pinned annotations PyEnsembl's EnsemblRelease "
                                        "doesn't cover, such as T2T-CHM13.", prog=prog)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Show pinned annotations, without downloading")
    listing.add_argument("--species", help="Only this species' annotations, e.g. human")
    add_json(listing)
    add_limit(listing)
    install = commands.add_parser("install", help="Download, check and index a pinned annotation")
    install.add_argument("name", help="Annotation name, as 'orfdata genomes list' shows it, e.g. t2t-chm13")
    install.add_argument("--force", action="store_true", help="Download again even when present")
    install.add_argument("--genome-fasta", help="Explicit chromosome FASTA path or URL to install through PyEnsembl")
    install.add_argument("--source-path", action="append", default=[], metavar="ROLE=PATH",
                         help="Install a local copy of the gtf, transcripts or proteins file instead of downloading")
    add_timeout(install, note="; PyEnsembl downloads --genome-fasta without one")
    add_json(install)
    share_cache_root(parser)
    return parser


def lift_parser(prog=None):
    """Return the parser of ``orfdata lift``, which maps an ORF catalogue to another genome assembly.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
    """
    parser = ArgumentParser(
        description="Lift an ORF catalogue to another assembly of its species with a chain file.", prog=prog,
        epilog="Without --chain, the pinned chain files, in as many hops as needed, must already be downloaded "
               "with 'orfdata chains download'. Placement on the new assembly's annotation re-checks each ORF's "
               "coding path and protein.")
    parser.add_argument("--to", required=True, metavar="BUILD", help="Target assembly, e.g. GRCh37")
    parser.add_argument("--catalog", help="Local JSON/JSON.gz catalogue to lift instead of a species' catalogue")
    parser.add_argument("--species", help="Lift this species' catalogue, e.g. human or mouse (default human)")
    parser.add_argument("--genome-build", metavar="BUILD",
                        help="Assembly of the species' catalogue to lift, when it has several")
    parser.add_argument("--chain", metavar="FILE",
                        help="UCSC or Ensembl chain file from the catalogue's assembly to --to, instead of the "
                             "pinned one; its direction is checked only when it is a pinned file")
    add_cache_root(parser, help="Cache directory holding catalogues and downloaded chain files, instead of "
                                 "ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", required=True, metavar="FILE",
                        help="Write the lifted catalogue here; gzipped when it ends in .gz")
    parser.add_argument("--json", action="store_true", help="Write the complete build summary as JSON to stderr")
    return parser


def chains_parser(prog=None):
    """Return the parser of ``orfdata chains``, which lists and downloads the pinned assembly chain files.

    Parameters
    ----------
    prog : str, optional
        The program name that help and errors show; argparse's default when None.

    Returns
    -------
    ArgumentParser
        Its subcommands are ``list`` and ``download``.
    """
    parser = ArgumentParser(
        description="List and download the pinned chain files 'orfdata lift' uses.", prog=prog)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Show pinned chain files and their terms, without downloading")
    listing.add_argument("--species", help="Only this species' chains, e.g. human or mouse")
    add_json(listing)
    add_limit(listing)
    download = commands.add_parser(
        "download", help="Download every pinned chain file (or a species', or one FROM TO) and verify checksums")
    download.add_argument("source_build", metavar="FROM", nargs="?", help="Assembly the chain maps from, e.g. GRCh38")
    download.add_argument("target_build", metavar="TO", nargs="?", help="Assembly the chain maps to, e.g. GRCh37")
    download.add_argument("--species", help="Only this species' chains, or the species when names are shared")
    download.add_argument("--accept-terms", action="store_true",
                          help="Accept the provider's usage terms, required for restricted files such as UCSC's")
    download.add_argument("--force", action="store_true", help="Download again even when present")
    download.add_argument("--source-path", metavar="FILE",
                          help="Install this local copy instead of downloading; it must match the pin")
    add_timeout(download)
    add_json(download)
    share_cache_root(parser)
    return parser
