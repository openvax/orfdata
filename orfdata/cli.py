"""Query and export the offline ORF evidence reference."""

import argparse
import csv
import json
import re
import sys

from .catalog import EVIDENCE_METHODS, UORFCatalog, load_catalog, merge_catalogs
from .evidence import MS, RiboSeq, evidence_rule


def run(argv=None, prog=None, *, catalog_loader=load_catalog):
    """Export reference records as JSON, TSV or FASTA."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "catalogs":
        from .catalog_cli import run as catalog_run

        return catalog_run(argv[1:], prog=(prog or "orfdata") + " catalogs")
    if argv and argv[0] == "datasets":
        from .dataset_cli import run as dataset_run

        return dataset_run(argv[1:], prog=(prog or "orfdata") + " datasets")
    if argv and argv[0] == "proteome":
        from .proteome_cli import run as proteome_run

        return proteome_run(argv[1:], prog=(prog or "orfdata") + " proteome")
    if argv and argv[0] == "setup":
        return _setup(argv[1:], prog=(prog or "orfdata") + " setup")
    if argv and argv[0] == "genomes":
        from .genomes_cli import run as genomes_run

        return genomes_run(argv[1:], prog=(prog or "orfdata") + " genomes")
    if argv and argv[0] in ("lift", "chains"):
        from . import liftover_cli

        command = liftover_cli.lift if argv[0] == "lift" else liftover_cli.chains
        return command(argv[1:], prog=f"{prog or 'orfdata'} {argv[0]}")
    parser = argparse.ArgumentParser(
        description=__doc__, prog=prog,
        epilog="Install a species' catalogues and chains with 'orfdata setup --species S'. "
               "Write a reference proteome with 'orfdata proteome --help'; manage optional ORF "
               "catalogues with 'orfdata catalogs --help'; lift a catalogue to another assembly with "
               "'orfdata lift --help' and 'orfdata chains --help'; annotations PyEnsembl lacks, such as "
               "T2T-CHM13, with 'orfdata genomes --help'; quantitative data with "
               "'orfdata datasets --help'.")
    _add_query_arguments(parser)
    args = parser.parse_args(argv)
    return _query(args, parser, catalog_loader)


def _setup(argv, prog):
    from .install import InstallError, install_species

    parser = argparse.ArgumentParser(
        prog=prog, description="Install a species' ORF catalogues and chain files, so its proteomes work offline "
                               "on every assembly the chains reach.")
    parser.add_argument("--species", required=True, help="e.g. human or mouse")
    parser.add_argument("--accept-terms", action="store_true",
                        help="Accept usage terms of restricted chain files, such as UCSC's non-commercial ones")
    parser.add_argument("--force", action="store_true", help="Install again even what is already installed")
    parser.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--timeout", type=float, default=60,
                        help="Seconds a download may stall before it is retried, then fails (default 60)")
    args = parser.parse_args(argv)
    try:
        result = install_species(args.species, accept_terms=args.accept_terms, force=args.force,
                                 cache_root=args.cache_root, timeout=args.timeout)
    except InstallError as error:  # What did install is still reported.
        json.dump(error.result, sys.stdout, indent=2)
        sys.stdout.write("\n")
        parser.exit(1, f"{parser.prog}: error: {error}\n")
    except (ValueError, OSError, ImportError) as error:
        parser.error(str(error))
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")


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
    parser.add_argument("--evidence", choices=EVIDENCE_METHODS)
    _add_rule_arguments(parser)
    parser.add_argument("--biotype", help="Source ORF class, as each source labels it: uORF, lncRNA, intORF, dORF "
                                          "(GENCODE, human); uORF, uoORF, iORF, dORF, doORF, lncORF (mouse Ribo-seq; "
                                          "RibORF also truncation, extension); AltProt, Isoform (OpenProt)")
    parser.add_argument("--region", help="1-based inclusive interval: chr7:66205483-66205522")
    parser.add_argument(
        "--position", help="Map a 1-based variant position: chr7:66205483 (JSON only)"
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
    parser.add_argument("--format", choices=("json", "tsv", "fasta"), default="json")


def _add_rule_arguments(parser):
    """Ribo-seq and MS requirements, shared by ORF queries and proteome builds."""
    parser.add_argument("--min-ribo-seq-studies", type=int, metavar="N",
                        help="Require Ribo-seq in at least N studies the source names")
    parser.add_argument("--min-ms-peptides", type=int, metavar="N",
                        help="Require at least N distinct supporting MS peptides")
    parser.add_argument("--ms-methods", metavar="KINDS",
                        help="Count only these MS peptide kinds, comma-separated: shotgun, mhc, hla, ms_unclassified")
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
        raise ValueError("--either needs both --min-ribo-seq-studies and an MS requirement")
    if not (ribo or ms):
        method = default if method is None else method
        return evidence_rule(None if method == "all" else method, require_unique=args.require_unique,
                             require_reviewed=args.require_reviewed, require_resolved=args.require_resolved)
    if method is not None:
        raise ValueError(f"Use {flag} or the Ribo-seq/MS requirements, not both")
    if not ms and (args.require_unique or args.require_reviewed):
        raise ValueError("--require-unique and --require-reviewed need an MS requirement")
    rules = [RiboSeq(args.min_ribo_seq_studies, resolved=args.require_resolved)] if ribo else []
    if ms:
        rules.append(MS(1 if args.min_ms_peptides is None else args.min_ms_peptides,
                        unique=args.require_unique, reviewed=args.require_reviewed,
                        resolved=args.require_resolved,
                        methods=args.ms_methods.split(",") if args.ms_methods is not None else None))
    if len(rules) == 1:
        return rules[0]
    return rules[0] | rules[1] if args.either else rules[0] & rules[1]


def _query(args, parser, catalog_loader):
    if args.region and args.position:
        parser.error("Choose --region or --position")
    if args.position and args.format != "json":
        parser.error("--position requires --format json")
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
            "genome_build": args.genome_build,
            "include_initiation_context": args.include_initiation_context or bool(args.ref),
        }
        locus = args.region or args.position
        if locus:
            pattern = r"([^:]+):([0-9]+)-([0-9]+)" if args.region else r"([^:]+):([0-9]+)"
            match = re.fullmatch(pattern, locus)
            if match is None:
                raise ValueError(f"Invalid locus: {locus}")
            start = int(match.group(2))
            options.update(
                contig=match.group(1),
                start=start,
                end=int(match.group(3)) if args.region else start,
            )
        records = catalog.query(**options)
        if args.orf_id:
            orf_id = catalog.get(args.orf_id).orf_id  # Any source's ID; a misspelled one is an error.
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
    except (ValueError, KeyError, OSError, ImportError) as exc:
        parser.error(str(exc))
    if args.format == "json":
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
    elif args.format == "fasta":
        for r in records:
            reference = "" if r.reference_context is None else " reference=" + r.reference_context.key
            scope = "" if isinstance(catalog, UORFCatalog) else f" taxon={r.taxon_id}{reference}"
            sys.stdout.write(
                f">{r.orf_id} gene={r.gene_name} transcript={r.transcript_id} assembly={r.genome_build} biotype={r.biotype}{scope}\n{r.protein_sequence}\n"
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
        provenance = ("members", "liftover", "lifted_from")
        if scoped:
            fields += ("start_codon", "reference_context_json", "summaries_json") + tuple(
                f + "_json" for f in provenance)
        elif any(r.summaries for r in records):
            # The released upstream-only columns have no place for aggregate counts.
            parser.error("The legacy upstream-ORF TSV cannot export aggregate evidence; use --format json")
        elif any(r.members or r.liftover for r in records):
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
