"""Explicit optional catalogue discovery, acquisition and offline queries."""

import json
import sys
from contextlib import redirect_stdout
from dataclasses import asdict
from pathlib import Path

from .catalog import catalog_status
from .catalogs import CatalogRegistry
from .cli import _add_query_arguments, _query
from .cli_output import manifest_summary, print_result
from .cli_parser import ArgumentParser, UsageError, add_json, add_limit


def run(argv=None, prog=None):
    parser = ArgumentParser(description=__doc__, prog=prog,
                            epilog="Start with 'orfdata catalogs list', then 'orfdata catalogs info NAME'.\n"
                                   "Use 'orfdata data fetch NAME' to select the published download or source build.")
    parser.add_argument("--registry", help="Local JSON catalogue registry definitions")
    parser.add_argument("--cache-root", help="Override the OpenVax cache root")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Discover release manifests without downloading")
    listing.add_argument("--species", help="Only this species' catalogues, e.g. human or mouse")
    listing.add_argument("--taxon-id", type=int)
    add_json(listing)
    add_limit(listing)
    status = commands.add_parser("status", help="Show bundled and installed catalogue status for every species")
    add_json(status)
    add_limit(status)
    query_parser = None
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
            subparser.add_argument("--timeout", type=float, default=60,
                                   help="Seconds a download may stall before it is retried, then fails (default 60)")
        if command == "query":
            query_parser = subparser
            _add_query_arguments(subparser, catalog_path=False)
        else:
            add_json(subparser)
    args = parser.parse_args(argv)
    try:
        definitions = json.loads(Path(args.registry).read_text()) if args.registry else None
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
            return _query(args, query_parser, lambda _: method(args.name, args.version))
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
