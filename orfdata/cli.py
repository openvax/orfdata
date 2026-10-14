"""Query and export the offline ORF evidence reference."""

import argparse
import csv
import json
import re
import sys

from .catalog import EVIDENCE_METHODS, UORFCatalog, load_catalog


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
    parser = argparse.ArgumentParser(
        description=__doc__, prog=prog,
        epilog="Write a reference proteome with 'orfdata proteome --help'; manage optional ORF "
               "catalogues with 'orfdata catalogs --help'; quantitative data with 'orfdata datasets --help'.")
    _add_query_arguments(parser)
    args = parser.parse_args(argv)
    return _query(args, parser, catalog_loader)


def _add_query_arguments(parser, *, catalog_path=True):
    if catalog_path:
        parser.add_argument("--catalog", help="Use a compatible local JSON/JSON.gz catalogue")
        parser.add_argument("--species", help="Use this species' catalogue, e.g. human or mouse (default human)")
    else:
        parser.set_defaults(catalog=None, species=None)
    parser.add_argument("--gene", help="Exact source gene symbol or identifier")
    parser.add_argument("--transcript", help="Exact source transcript ID, including its version when supplied")
    parser.add_argument("--orf-id", help="Exact catalogue ORF identifier")
    parser.add_argument("--evidence", choices=EVIDENCE_METHODS)
    parser.add_argument("--require-resolved", action="store_true",
                        help="Require resolved peptides or named ribosome studies; exclude aggregate-only support")
    parser.add_argument(
        "--require-unique", action="store_true", help="Require one protein mapping in the MS source"
    )
    parser.add_argument(
        "--require-reviewed",
        action="store_true",
        help="Require convincing reviewed spectra (benchmark scores 4-5)",
    )
    parser.add_argument("--biotype", help="Source ORF class, e.g. uORF, lncRNA, intORF, dORF")
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
        catalog = (catalog_loader(args.catalog, species=args.species, genome_build=args.genome_build)
                   if args.species else catalog_loader(args.catalog))
        options = {
            "gene": args.gene,
            "transcript_id": args.transcript,
            "evidence": args.evidence,
            "require_resolved": args.require_resolved,
            "require_unique": args.require_unique,
            "require_reviewed": args.require_reviewed,
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
            catalog.get(args.orf_id)  # A misspelled ID is an error, not empty evidence.
            records = tuple(r for r in records if r.orf_id == args.orf_id)
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
    except (ValueError, KeyError, OSError) as exc:
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
        if scoped:
            fields += ("reference_context_json", "summaries_json")
        elif any(r.summaries for r in records):
            # The released upstream-only columns have no place for aggregate counts.
            parser.error("The legacy upstream-ORF TSV cannot export aggregate evidence; use --format json")
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
            writer.writerow(row)


if __name__ == "__main__":
    run()
