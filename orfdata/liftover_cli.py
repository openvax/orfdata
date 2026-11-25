"""Lift an ORF catalogue to another assembly, and manage the pinned chain files that do it."""

import sys
from collections import Counter
from contextlib import redirect_stdout
from dataclasses import asdict
from pathlib import Path

from .catalog import load_catalog
from .catalogs import CatalogRegistry
from .cli_output import print_result, print_run_summary
from .cli_parser import UsageError
from .command_parsers import chains_parser, lift_parser
from .liftover import (
    ChainFile,
    chain_downloaded,
    download_chain,
    download_chains,
    lift_catalog,
    lift_to,
    pinned_chains,
)


def lift(argv=None, prog=None):
    parser = lift_parser(prog)
    args = parser.parse_args(argv)
    try:
        if args.catalog and (args.species or args.genome_build):
            raise UsageError("Pass --catalog or --species/--genome-build, not both")
        if args.summary is not None and Path(args.summary).resolve() == Path(args.output).resolve():
            raise UsageError("--summary and --output name the same file")
        registry = CatalogRegistry(cache_root=args.cache_root) if args.cache_root else None
        catalog = (load_catalog(args.catalog) if args.catalog
                   else load_catalog(species=args.species, genome_build=args.genome_build, registry=registry))
        if args.to == catalog.genome_build:
            raise ValueError(f"The catalogue is already on {args.to}")
        lifted = (lift_catalog(catalog, ChainFile(args.chain, source_build=catalog.genome_build, target_build=args.to))
                  if args.chain else lift_to(catalog, args.to, cache_root=args.cache_root))
        lifted.save(args.output)
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    hops, metadata = [], lifted.metadata
    while metadata.get("liftover") and metadata is not catalog.metadata:  # This command's hops, last first.
        hops.insert(0, metadata["liftover"])
        if metadata["liftover"]["source"] == catalog.source:
            break
        metadata = metadata["liftover"]["input"]["metadata"]
    unmapped = [u for hop in hops for u in hop["unmapped"]]
    summary = dict(source=catalog.source, lifted=lifted.source,
                   chains=[dict(chain=h["chain"], sha256=h["sha256"], lifted=h["lifted"]) for h in hops],
                   records=len(catalog.records), lifted_records=len(lifted.records), unmapped=len(unmapped),
                   unmapped_reasons=dict(Counter(u["reason"] for u in unmapped)), output=args.output)
    try:
        print_run_summary(summary, json_output=args.json, path=args.summary)
    except OSError as error:
        parser.fail(error)


def require_terms(source_build, target_build, species, accept_terms):
    """Raise :class:`~orfdata.errors.TermsNotAccepted` when a restricted chain lacks ``accept_terms``."""
    from .errors import TermsNotAccepted
    from .liftover import _pinned

    chain = _pinned(source_build, target_build, species)
    if chain.restricted and not accept_terms:
        raise TermsNotAccepted(
            f"The {chain.name} chain file requires --accept-terms: {chain.terms} ({chain.terms_url})"
        )


def chains(argv=None, prog=None):
    parser = chains_parser(prog)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            output = []
            for chain in pinned_chains(species=args.species):
                downloaded = chain_downloaded(chain, cache_root=args.cache_root)
                output.append(dict(asdict(chain), downloaded=downloaded))
        elif args.source_build is None:
            if args.source_path:
                raise UsageError("--source-path installs one chain; name it with FROM TO")
            with redirect_stdout(sys.stderr):
                result = download_chains(species=args.species, accept_terms=args.accept_terms, force=args.force,
                                         cache_root=args.cache_root, download_options=dict(timeout=args.timeout))
            output = dict(downloaded={k: str(v) for k, v in result["downloaded"].items()}, skipped=result["skipped"])
        elif args.target_build is None:
            raise UsageError("Name both assemblies, FROM and TO, or neither to download every pinned chain")
        else:
            with redirect_stdout(sys.stderr):
                output = dict(path=str(download_chain(args.source_build, args.target_build, species=args.species,
                                                      accept_terms=args.accept_terms, force=args.force,
                                                      source_path=args.source_path, cache_root=args.cache_root,
                                                      download_options=dict(timeout=args.timeout))))
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    print_result(output, json_output=args.json,
                 columns=("source_build", "target_build", "species", "provider", "downloaded", "restricted"),
                 limit=getattr(args, "limit", 25))
    if not args.json and args.command == "list":
        print("\nInstall: orfdata chains download FROM TO\nRestricted chains require --accept-terms.")
