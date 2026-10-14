"""Explicit optional catalogue discovery, acquisition and offline queries."""

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .catalog import catalog_status
from .catalogs import CatalogRegistry
from .cli import _add_query_arguments, _query


def run(argv=None, prog=None):
    parser = argparse.ArgumentParser(description=__doc__, prog=prog)
    parser.add_argument("--registry", help="Local JSON catalogue registry definitions")
    parser.add_argument("--cache-root", help="Override the OpenVax cache root")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Discover release manifests without downloading")
    listing.add_argument("--species", help="Only this species' catalogues, e.g. human or mouse")
    listing.add_argument("--taxon-id", type=int)
    commands.add_parser("status", help="Every species' catalogues, including the bundled one, and what is installed")
    query_parser = None
    for command in ("info", "inspect", "download", "build", "query"):
        subparser = commands.add_parser(command)
        subparser.add_argument("name", help="Catalogue name, as 'orfdata catalogs list' shows it")
        subparser.add_argument("--version", help="Release to use instead of the catalogue's default")
        if command in ("inspect", "query"):
            subparser.add_argument("--built", action="store_true", help="Select a local source build explicitly")
        if command in ("download", "build"):
            subparser.add_argument("--force", action="store_true", help="Replace an existing installation")
            subparser.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH",
                                   help="Use this local copy of a pinned source asset instead of downloading it")
        if command == "query":
            query_parser = subparser
            _add_query_arguments(subparser, catalog_path=False)
    args = parser.parse_args(argv)
    try:
        definitions = json.loads(Path(args.registry).read_text()) if args.registry else None
        registry = CatalogRegistry(definitions, cache_root=args.cache_root)
        if args.command == "list":
            output = [m.to_dict() for m in registry.list(species=args.species, taxon_id=args.taxon_id)]
        elif args.command == "status":
            output = [asdict(row) for row in catalog_status(registry=registry)]
        elif args.command == "info":
            output = registry.info(args.name, args.version).to_dict()
        elif args.command == "inspect":
            method = registry.inspect_build if args.built else registry.inspect
            inspection = method(args.name, args.version)
            output = dict(status=inspection.status, path=str(inspection.path), verified=inspection.verified)
            if inspection.error is not None:
                output["error"] = str(inspection.error)
        elif args.command in ("download", "build"):
            paths = {}
            for item in args.source_path:
                asset, separator, path = item.partition("=")
                if not separator or not asset or not path or asset in paths:
                    raise ValueError("Source overrides require unique ASSET=PATH entries")
                paths[asset] = path
            method = registry.build if args.command == "build" else registry.download
            method(args.name, args.version, force=args.force, source_paths=paths or None)
            inspect = registry.inspect_build if args.command == "build" else registry.inspect
            result = inspect(args.name, args.version)
            output = dict(status=result.status, path=str(result.path), verified=result.verified)
        else:
            method = registry.open_build if args.built else registry.open
            return _query(args, query_parser, lambda _: method(args.name, args.version))
    except (ValueError, TypeError, KeyError, OSError, ImportError) as error:
        parser.error(str(error))
    print(json.dumps(output, indent=2, sort_keys=True))
