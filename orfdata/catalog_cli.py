"""Explicit optional catalogue discovery, acquisition and offline queries."""

import sys
from contextlib import redirect_stdout
from dataclasses import asdict

from .catalog import catalog_status
from .catalogs import CatalogRegistry
from .cli import _query
from .cli_output import manifest_summary, print_result
from .cli_parser import parse_source_paths
from .command_parsers import catalogs_parser
from .downloads import normalize_status
from .registry import load_registry


def run(argv=None, prog=None):
    parser = catalogs_parser(prog)
    args = parser.parse_args(argv)
    try:
        registry = load_registry(CatalogRegistry, args.registry, cache_root=args.cache_root)
        if args.command == "list":
            output = [m.to_dict() for m in registry.list(species=args.species, taxon_id=args.taxon_id)]
        elif args.command == "status":
            output = [asdict(row) for row in catalog_status(registry=registry)]
        elif args.command == "info":
            output = registry.info(args.name, args.version).to_dict()
        elif args.command == "inspect":
            method = registry.inspect_build if args.built else registry.inspect
            inspection = method(args.name, args.version)
            output = dict(
                status=normalize_status(inspection.status),
                datacache_status=inspection.status,
                path=str(inspection.path),
                verified=inspection.verified,
            )
            if inspection.error is not None:
                output["error"] = str(inspection.error)
        elif args.command in ("download", "build"):
            paths = parse_source_paths(args.source_path, metavar="ASSET=PATH")
            options = dict(timeout=args.timeout)
            with redirect_stdout(sys.stderr):
                if args.command == "build":
                    registry.build(args.name, args.version, force=args.force, source_paths=paths,
                                   download_options=options)
                else:
                    registry.download(args.name, args.version, force=args.force, source_paths=paths, **options)
            inspect = registry.inspect_build if args.command == "build" else registry.inspect
            result = inspect(args.name, args.version)
            output = dict(
                status=normalize_status(result.status),
                datacache_status=result.status,
                path=str(result.path),
            )
            if args.command == "download":  # A local build is checked by its receipt, not a published pin.
                output["verified"] = result.verified
        else:
            method = registry.open_build if args.built else registry.open
            # The parser reports errors under the subcommand it parsed, 'orfdata catalogs query'.
            return _query(args, parser, lambda _: method(args.name, args.version))
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    columns = ("catalog_id", "version", "title") if args.command == "list" else (
        "species", "name", "version", "genome_build", "status", "acquire")
    if args.command == "info" and not args.json:
        output = manifest_summary(output)
    print_result(output, json_output=args.json, columns=columns, limit=getattr(args, "limit", 25))
    if args.command == "info" and not args.json:
        print("\nUse --json for full source pins, citations, and processing metadata.")
    if not args.json and args.command in ("list", "status"):
        print("\nDetails: orfdata catalogs info NAME\nInstall: orfdata data fetch NAME")
