"""Discover, download, inspect and query quantitative datasets."""

import json
import sqlite3
import sys
from contextlib import redirect_stdout

from .cli_output import manifest_summary, print_fields, print_record_stream, print_result
from .cli_parser import ArgumentParser, add_json, add_limit, add_timeout, share_cache_root
from .datasets import DatasetRegistry
from .reference import ScopedTarget


def run(argv=None, prog=None):
    parser = ArgumentParser(description=__doc__, prog=prog,
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
    args = parser.parse_args(argv)
    try:
        definitions = None
        if args.registry:
            with open(args.registry) as handle:
                definitions = json.load(handle)
        registry = DatasetRegistry(definitions, cache_root=args.cache_root)
        if args.command == "list":
            result = [m.to_dict() for m in registry.list(
                species=args.species, taxon_id=args.taxon_id, assay=args.assay)]
        elif args.command == "info":
            result = registry.info(args.name, args.version).to_dict()
            result["download"] = registry.download_metadata(args.name, args.version)
        elif args.command == "references":
            result = [spec.to_dict() for spec in registry.search_references(args.name, args.version)]
        elif args.command == "inspect-reference":
            inspected = registry.inspect_reference(args.name, args.version, reference=args.reference, asset_id=args.asset)
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
        elif args.command == "download-reference":
            options = dict(show_progress=not args.quiet, timeout=args.timeout)
            with redirect_stdout(sys.stderr):
                result = registry.download_reference(args.name, args.version, reference=args.reference, asset_id=args.asset,
                                                     force=args.force, source_path=args.source_path, **options)
        elif args.command == "query-reference":
            reference = registry.open_reference(args.name, args.version, reference=args.reference, asset_id=args.asset)
            target = ScopedTarget(reference.reference, args.kind, args.target, args.namespace or reference.spec.namespace)
            result = dict(source=reference.spec.to_dict(), membership=reference.membership(target).to_dict())
            if not args.json:
                membership = result["membership"]
                result = dict(target=args.target, namespace=target.namespace, status=membership["status"],
                              reason=membership["reason"], reference=target.reference.key,
                              source_asset=reference.spec.source.asset_id)
        elif args.command == "inspect":
            inspect = registry.inspect_build if args.built else registry.inspect
            inspected = inspect(args.name, args.version)
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
            if args.built or not registry.published(args.name, args.version):
                result.update(origin="source_build", verified=inspected.verified)
        elif args.command == "download":
            options = dict(show_progress=not args.quiet, timeout=args.timeout)
            sources = {"dataset.sqlite": args.source_path} if args.source_path else None
            with redirect_stdout(sys.stderr):
                result = registry.download(args.name, args.version, force=args.force,
                                           source_paths=sources, **options)
        elif args.command == "build":
            sources = {}
            for value in args.source_path:
                asset, separator, path = value.partition("=")
                if not separator or not asset or not path or asset in sources:
                    parser.error("Source overrides require unique ASSET=PATH entries")
                sources[asset] = path
            options = dict(show_progress=not args.quiet, timeout=args.timeout)
            with redirect_stdout(sys.stderr):
                dataset = registry.build(args.name, args.version, force=args.force,
                                         source_paths=sources, download_options=options)
            result = dict(path=str(dataset.path), counts=dataset.counts,
                          dataset_id=dataset.manifest.dataset_id, version=dataset.manifest.version)
        elif args.command == "pairs":
            open_dataset = registry.open_build if args.built else registry.open
            dataset = open_dataset(args.name, args.version)
            references = {r.key: r for r in dataset.manifest.references}
            selections = {}
            for side in ("left", "right"):
                identifier, kind, namespace, reference = (
                    getattr(args, f"{side}_{field}") for field in ("target", "kind", "namespace", "reference"))
                if any((identifier, kind, namespace, reference)):
                    if not all((identifier, kind, namespace, reference)):
                        parser.error(f"{side.title()} target queries require --{side}-target, --{side}-kind, "
                                     f"--{side}-namespace and --{side}-reference")
                    selections[f"{side}_target"] = ScopedTarget(references[reference], kind, identifier, namespace)
                selections[f"{side}_filters"] = {
                    field: getattr(args, f"{side}_{field}")
                    for field in ("quantity", "status", "unit", "scale", "log_base")
                    if getattr(args, f"{side}_{field}") is not None}
            matches = dataset.matched_measurements(
                allowed_bases=args.basis, allowed_temporal=args.temporal or ["concurrent"],
                reference=args.reference, pair_ids=args.pair_id, relation_ids=args.relation_id,
                **selections,
            )
            print_record_stream((match.to_table_row() for match in matches), format=args.format, limit=args.limit,
                                columns=(("PAIR", "pairing.pair_id"), ("BASIS", "pairing.basis"),
                                         ("TEMPORAL", "pairing.temporal_relation"),
                                         ("RELATION", "target_relation.relation"),
                                         ("LEFT REF", "target_relation.left.reference"),
                                         ("LEFT", "target_relation.left.identifier"),
                                         ("RIGHT REF", "target_relation.right.reference"),
                                         ("RIGHT", "target_relation.right.identifier")))
            return
        else:
            open_dataset = registry.open_build if args.built else registry.open
            dataset = open_dataset(args.name, args.version)
            target = None
            if any((args.target, args.kind, args.namespace)):
                if not all((args.target, args.kind, args.namespace, args.reference)):
                    parser.error("Target queries require --target, --kind, --namespace and --reference")
                references = {r.key: r for r in dataset.manifest.references}
                target = ScopedTarget(references[args.reference], args.kind, args.target, args.namespace)
            filters = dict(sample_id=args.sample, assay_id=args.assay_id,
                           quantity=args.quantity, status=args.status,
                           unit=args.unit, scale=args.scale, log_base=args.log_base)
            if args.mapped:
                if target is None:
                    parser.error("--mapped requires an exact target and reference")
                measurements = dataset.mapped_measurements(
                    target, mapping_statuses=args.mapping_status or ["unique"],
                    require_sequence=args.require_sequence, **filters)
            else:
                if args.mapping_status or args.require_sequence:
                    parser.error("Mapping filters require --mapped")
                measurements = dataset.measurements(target=target, reference=args.reference, **filters)
            columns = (("REF", "target.reference"), ("TARGET", "target.identifier"), ("SAMPLE", "sample_id"),
                       ("QUANTITY", "quantity"), ("VALUE", "value"),
                       ("UNIT", "unit"), ("SCALE", "scale"), ("STATUS", "status"))
            if args.mapped:
                columns = tuple((label, "measurement." + path) for label, path in columns) + (
                    ("MAPPING", "attribution"),)
                if args.format == "table":
                    print_fields(dict(requested_target=target.identifier, namespace=target.namespace,
                                      requested_reference=target.reference.key))
                    print("Values retain their original targets; group membership supplies context only.\n")
            print_record_stream((measurement.to_dict() for measurement in measurements),
                                format=args.format, limit=args.limit, columns=columns)
            return
        columns = (("dataset_id", "version", "assay_kinds") if args.command == "list"
                   else ("source.asset_id", "namespace", "reference.species", "reference.genome_build"))
        if args.command == "info" and not args.json:
            result = dict(manifest_summary(result), download=result["download"])
        print_result(result, json_output=args.json, columns=columns, limit=getattr(args, "limit", 25))
        if args.command in ("info", "query-reference") and not args.json:
            print("\nUse --json for full source, reference, and provenance metadata.")
        if not args.json and args.command == "list":
            print("\nDetails: orfdata datasets info NAME\nInstall: orfdata data fetch NAME")
    except (ValueError, TypeError, KeyError, OSError, ImportError, RuntimeError, sqlite3.Error) as exc:
        parser.error(str(exc))
