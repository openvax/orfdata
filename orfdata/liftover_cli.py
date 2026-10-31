"""Lift an ORF catalogue to another assembly, and manage the pinned chain files that do it."""

import sys
from collections import Counter
from contextlib import redirect_stdout
from dataclasses import asdict

from .catalog import load_catalog
from .catalogs import CatalogRegistry
from .cli_output import print_result
from .cli_parser import (
    ArgumentParser,
    add_cache_root,
    add_json,
    add_limit,
    add_timeout,
    share_cache_root,
)
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
    parser = ArgumentParser(
        description="Lift an ORF catalogue to another assembly of its species with a chain file.", prog=prog,
        epilog="Without --chain, the pinned chain files, in as many hops as needed, must already be downloaded "
               "with 'orfdata chains download'. Placement on the new assembly's annotation re-checks each ORF's "
               "coding path and protein.")
    parser.add_argument("--to", required=True, metavar="BUILD", help="Target assembly, e.g. GRCh37")
    parser.add_argument("--catalog", help="Local JSON/JSON.gz catalogue to lift instead of a species' catalogue")
    parser.add_argument("--species", help="Lift this species' catalogue, e.g. human or mouse (default human)")
    parser.add_argument("--genome-build", metavar="BUILD",
                        help="Assembly of the species' catalogue to lift, when it has several")
    parser.add_argument("--chain", metavar="FILE",
                        help="UCSC or Ensembl chain file from the catalogue's assembly to --to, instead of the "
                             "pinned one; its direction is checked only when it is a pinned file")
    add_cache_root(parser, help="Cache directory holding catalogues and downloaded chain files, instead of "
                                 "ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", required=True, metavar="FILE",
                        help="Write the lifted catalogue here; gzipped when it ends in .gz")
    parser.add_argument("--json", action="store_true", help="Write the complete build summary as JSON to stderr")
    args = parser.parse_args(argv)
    try:
        if args.catalog and (args.species or args.genome_build):
            raise ValueError("Pass --catalog or --species/--genome-build, not both")
        registry = CatalogRegistry(cache_root=args.cache_root) if args.cache_root else None
        catalog = (load_catalog(args.catalog) if args.catalog
                   else load_catalog(species=args.species, genome_build=args.genome_build, registry=registry))
        if args.to == catalog.genome_build:
            raise ValueError(f"The catalogue is already on {args.to}")
        lifted = (lift_catalog(catalog, ChainFile(args.chain, source_build=catalog.genome_build, target_build=args.to))
                  if args.chain else lift_to(catalog, args.to, cache_root=args.cache_root))
        lifted.save(args.output)
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.error(str(error))
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
    print_result(summary, json_output=args.json, stream=sys.stderr)


def chains(argv=None, prog=None):
    parser = ArgumentParser(
        description="List and download the pinned chain files 'orfdata lift' uses.", prog=prog)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Show pinned chain files and their terms, without downloading")
    listing.add_argument("--species", help="Only this species' chains, e.g. human or mouse")
    add_json(listing)
    add_limit(listing)
    download = commands.add_parser(
        "download", help="Download every pinned chain file (or a species', or one FROM TO) and verify checksums")
    download.add_argument("source_build", metavar="FROM", nargs="?", help="Assembly the chain maps from, e.g. GRCh38")
    download.add_argument("target_build", metavar="TO", nargs="?", help="Assembly the chain maps to, e.g. GRCh37")
    download.add_argument("--species", help="Only this species' chains, or the species when names are shared")
    download.add_argument("--accept-terms", action="store_true",
                          help="Accept the provider's usage terms, required for restricted files such as UCSC's")
    download.add_argument("--force", action="store_true", help="Download again even when present")
    download.add_argument("--source-path", metavar="FILE",
                          help="Install this local copy instead of downloading; it must match the pin")
    add_timeout(download)
    add_json(download)
    share_cache_root(parser)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            output = []
            for chain in pinned_chains(species=args.species):
                downloaded = chain_downloaded(chain, cache_root=args.cache_root)
                output.append(dict(asdict(chain), downloaded=downloaded))
        elif args.source_build is None:
            if args.source_path:
                raise ValueError("--source-path installs one chain; name it with FROM TO")
            with redirect_stdout(sys.stderr):
                result = download_chains(species=args.species, accept_terms=args.accept_terms, force=args.force,
                                         cache_root=args.cache_root, download_options=dict(timeout=args.timeout))
            output = dict(downloaded={k: str(v) for k, v in result["downloaded"].items()}, skipped=result["skipped"])
        elif args.target_build is None:
            raise ValueError("Name both assemblies, FROM and TO, or neither to download every pinned chain")
        else:
            with redirect_stdout(sys.stderr):
                output = dict(path=str(download_chain(args.source_build, args.target_build, species=args.species,
                                                      accept_terms=args.accept_terms, force=args.force,
                                                      source_path=args.source_path, cache_root=args.cache_root,
                                                      download_options=dict(timeout=args.timeout))))
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.error(str(error))
    print_result(output, json_output=args.json,
                 columns=("source_build", "target_build", "species", "provider", "downloaded", "restricted"),
                 limit=getattr(args, "limit", 25))
    if not args.json and args.command == "list":
        print("\nInstall: orfdata chains download FROM TO\nRestricted chains require --accept-terms.")
