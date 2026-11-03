"""Explicit optional catalogue discovery, acquisition and offline queries."""

import json
import sys
from contextlib import redirect_stdout
from dataclasses import asdict
from pathlib import Path

from .catalog import catalog_status
from .catalogs import CatalogRegistry
from .cli import _query
from .cli_output import manifest_summary, print_result
from .cli_parser import UsageError
from .command_parsers import catalogs_parser


def run(argv=None, prog=None):
    parser = catalogs_parser(prog)
    args = parser.parse_args(argv)
    try:
        definitions = None
        if args.registry:
            try:
                definitions = json.loads(Path(args.registry).read_text())
            except ValueError as error:
                raise ValueError(f"Invalid catalogue registry {args.registry}: {error}") from error
        try:
            registry = CatalogRegistry(definitions, cache_root=args.cache_root)
        except (KeyError, TypeError, AttributeError) as error:
            if not args.registry:  # The packaged registry: a bug, which keeps its traceback.
                raise
            raise ValueError(f"Invalid catalogue registry {args.registry}: {type(error).__name__}: {error}") from error
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
                    raise UsageError("Source overrides require unique ASSET=PATH entries")
                paths[asset] = path
            options = dict(timeout=args.timeout)
            with redirect_stdout(sys.stderr):
                if args.command == "build":
                    registry.build(args.name, args.version, force=args.force, source_paths=paths or None,
                                   download_options=options)
                else:
                    registry.download(args.name, args.version, force=args.force, source_paths=paths or None, **options)
            inspect = registry.inspect_build if args.command == "build" else registry.inspect
            result = inspect(args.name, args.version)
            output = dict(status=result.status, path=str(result.path))
            if args.command == "download":  # A local build is checked by its receipt, not a published pin.
                output["verified"] = result.verified
        else:
            method = registry.open_build if args.built else registry.open
            # The parser reports errors under the subcommand it parsed, 'orfdata catalogs query'.
            return _query(args, parser, lambda _: method(args.name, args.version))
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    columns = ("catalog_id", "version", "title") if args.command == "list" else (
        "species", "name", "version", "genome_build", "installed", "acquire")
    if args.command == "info" and not args.json:
        output = manifest_summary(output)
    print_result(output, json_output=args.json, columns=columns, limit=getattr(args, "limit", 25))
    if args.command == "info" and not args.json:
        print("\nUse --json for full source pins, citations, and processing metadata.")
    if not args.json and args.command in ("list", "status"):
        print("\nDetails: orfdata catalogs info NAME\nInstall: orfdata data fetch NAME")
