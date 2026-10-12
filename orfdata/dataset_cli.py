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
    for command in ("info", "download", "inspect", "query"):
        subparser = commands.add_parser(command)
        subparser.add_argument("name")
        subparser.add_argument("--version")
        if command == "download":
            subparser.add_argument("--force", action="store_true")
            subparser.add_argument("--timeout", type=float)
            subparser.add_argument("--quiet", action="store_true")
            subparser.add_argument("--source-path", help="Use a local dataset.sqlite with the registered hash")
        if command == "query":
            subparser.add_argument("--reference", help="Exact ref:SHA256 context key")
            subparser.add_argument("--sample")
            subparser.add_argument("--assay-id")
            subparser.add_argument("--quantity")
            subparser.add_argument("--status")
            subparser.add_argument("--target")
            subparser.add_argument("--kind")
            subparser.add_argument("--namespace")
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
        elif args.command == "inspect":
            inspected = registry.inspect(args.name, args.version)
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
        elif args.command == "download":
            options = dict(show_progress=not args.quiet)
            if args.timeout is not None:
                options["timeout"] = args.timeout
            sources = {"dataset.sqlite": args.source_path} if args.source_path else None
            result = registry.download(args.name, args.version, force=args.force,
                                       source_paths=sources, **options)
        else:
            dataset = registry.open(args.name, args.version)
            target = None
            if any((args.target, args.kind, args.namespace)):
                if not all((args.target, args.kind, args.namespace, args.reference)):
                    parser.error("Target queries require --target, --kind, --namespace and --reference")
                references = {r.key: r for r in dataset.manifest.references}
                target = ScopedTarget(references[args.reference], args.kind, args.target, args.namespace)
            for measurement in dataset.measurements(
                target=target, reference=args.reference, sample_id=args.sample,
                assay_id=args.assay_id, quantity=args.quantity, status=args.status,
            ):
                sys.stdout.write(json.dumps(measurement.to_dict()) + "\n")
            return
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
    except (ValueError, TypeError, KeyError, OSError, ImportError, sqlite3.Error) as exc:
        parser.error(str(exc))
