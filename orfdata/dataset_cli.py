"""Discover, download, inspect and query quantitative datasets."""

import argparse
import json
import sqlite3
import sys

from .datasets import DatasetRegistry
from .reference import ScopedTarget


def run(argv=None, prog=None):
    parser = argparse.ArgumentParser(description=__doc__, prog=prog)
    parser.add_argument("--registry", help="Local registry JSON (otherwise use packaged releases)")
    parser.add_argument("--cache-root", help="Override the OpenVax data cache root")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="List published releases and assay coverage")
    listing.add_argument("--species")
    listing.add_argument("--taxon-id", type=int)
    listing.add_argument("--assay")
    for command in ("info", "download", "build", "inspect", "query", "pairs", "references",
                    "download-reference", "inspect-reference", "query-reference"):
        subparser = commands.add_parser(command)
        subparser.add_argument("name")
        subparser.add_argument("--version")
        if command in ("inspect", "query", "pairs"):
            subparser.add_argument("--built", action="store_true",
                                   help="Use a local source build explicitly, instead of a published artifact")
        if command in ("download", "download-reference"):
            subparser.add_argument("--force", action="store_true")
            subparser.add_argument("--timeout", type=float)
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
            subparser.add_argument("--timeout", type=float)
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
            options = dict(show_progress=not args.quiet)
            if args.timeout is not None:
                options["timeout"] = args.timeout
            result = registry.download_reference(args.name, args.version, reference=args.reference, asset_id=args.asset,
                                                 force=args.force, source_path=args.source_path, **options)
        elif args.command == "query-reference":
            reference = registry.open_reference(args.name, args.version, reference=args.reference, asset_id=args.asset)
            target = ScopedTarget(reference.reference, args.kind, args.target, args.namespace or reference.spec.namespace)
            result = dict(source=reference.spec.to_dict(), membership=reference.membership(target).to_dict())
        elif args.command == "inspect":
            inspect = registry.inspect_build if args.built else registry.inspect
            inspected = inspect(args.name, args.version)
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
            if args.built:
                result.update(origin="source_build", verified=inspected.verified)
        elif args.command == "download":
            options = dict(show_progress=not args.quiet)
            if args.timeout is not None:
                options["timeout"] = args.timeout
            sources = {"dataset.sqlite": args.source_path} if args.source_path else None
            result = registry.download(args.name, args.version, force=args.force,
                                       source_paths=sources, **options)
        elif args.command == "build":
            sources = {}
            for value in args.source_path:
                asset, separator, path = value.partition("=")
                if not separator or not asset or not path or asset in sources:
                    parser.error("Source overrides require unique ASSET=PATH entries")
                sources[asset] = path
            options = dict(show_progress=not args.quiet)
            if args.timeout is not None:
                options["timeout"] = args.timeout
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
            for match in dataset.matched_measurements(
                allowed_bases=args.basis, allowed_temporal=args.temporal or ["concurrent"],
                reference=args.reference, pair_ids=args.pair_id, relation_ids=args.relation_id,
                **selections,
            ):
                sys.stdout.write(json.dumps(match.to_table_row()) + "\n")
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
            for measurement in measurements:
                sys.stdout.write(json.dumps(measurement.to_dict()) + "\n")
            return
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
    except (ValueError, TypeError, KeyError, OSError, ImportError, sqlite3.Error) as exc:
        parser.error(str(exc))
