"""Discover and acquire ORFdata resources through one command."""

import argparse

from .cli_output import print_fields, print_json, print_table
from .downloads import list_downloads


def _source_paths(items, label):
    paths = {}
    for item in items:
        key, separator, path = item.partition("=")
        if not separator or not key or not path or key in paths:
            raise ValueError(f"{label} require unique NAME=PATH entries")
        paths[key] = path
    return paths or None


def _add_cache(parser):
    parser.add_argument(
        "--cache-root", help="Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")


def _add_force_timeout(parser, timeout=60):
    parser.add_argument("--force", action="store_true", help="Replace or repair an existing installation")
    parser.add_argument("--timeout", type=float, default=timeout, help="Per-attempt timeout in seconds")


def _add_json(parser):
    parser.add_argument("--json", action="store_true", help="Write machine-readable JSON")


def _show(rows, json_output):
    if json_output:
        print_json([row.to_dict() for row in rows])
        return
    headings = ("KIND", "NAME", "VERSION", "SPECIES", "CONTEXTS", "STATUS", "ACQUIRE", "COMMAND")
    values = []
    for row in rows:
        version = row.version or "-"
        if len(version) > 16:
            version = version[:12] + "\u2026"
        values.append((row.kind, row.name, version, ",".join(row.species),
                       ",".join(row.contexts), row.status, row.acquire, row.command))
    print_table(headings, values, empty="No matching resources")
    installed = [row for row in rows if row.installed and row.path]
    if installed:
        print("\nInstalled paths")
        for row in installed:
            print(f"{row.kind}/{row.name}: {row.path}")


def _show_result(result, json_output):
    if json_output:
        print_json(result)
    else:
        print_fields(result)


def _add_discovery(parser, *, status=False):
    _add_cache(parser)
    parser.add_argument(
        "--kind", choices=("viral", "catalog", "dataset", "reference", "genome", "chain"))
    parser.add_argument("--species", help="Only resources for this species, e.g. human or mouse")
    parser.add_argument("--context", action="append",
                        help="Only resources in this context; repeat for a union")
    parser.add_argument("--context-match", choices=("any", "all"), default="any",
                        help="Match any selected context (default) or require all")
    _add_json(parser)
    parser.add_argument("--verify", action="store_true", help="Hash and validate installed files")
    if status:
        parser.set_defaults(installed=True)
    else:
        selection = parser.add_mutually_exclusive_group()
        selection.add_argument("--installed", dest="installed", action="store_const", const=True,
                               help="Show only locally installed or bundled resources")
        selection.add_argument("--missing", dest="installed", action="store_const", const=False,
                               help="Show only missing, corrupt or inaccessible resources")
        parser.set_defaults(installed=None)


def run(argv=None, prog="orfdata download"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    listing = commands.add_parser("list", help="List every bundled, downloadable or buildable resource")
    _add_discovery(listing)
    status = commands.add_parser("status", help="Show resources already available locally")
    _add_discovery(status, status=True)

    viral = commands.add_parser(
        "viral", help="Download the current human-targeted viral collection and paired FASTAs")
    viral.add_argument("--species", default="human",
                       help="Target analysis species (currently human; future collections may differ)")
    viral.add_argument(
        "--collection", default="gdc-viruses",
        help="Reference collection (default gdc-viruses: 200 GDC references for human analysis)")
    viral.add_argument("--context", action="append",
                       help="Build an index for this context; repeat for a union")
    viral.add_argument("--exclude-context", action="append",
                       help="Exclude this context from the derived index; repeat as needed")
    viral.add_argument("--context-match", choices=("any", "all"), default="any",
                       help="Match any selected context (default) or require all")
    viral.add_argument("--subset", action="append",
                       help="Build an index for this named viral subset; repeat for a union")
    viral.add_argument("--inventory-molecules", dest="prefer_reference_molecule",
                       action="store_false",
                       help="Carry the molecule each inventory row names first, even where the "
                            "reference the collection follows was built from the other one")
    viral.add_argument("--no-annotation-fill", dest="annotate_from_alternate",
                       action="store_false",
                       help="Keep strictly to each entry's own record, leaving one whose source "
                            "annotates nothing without transcripts")
    viral.add_argument("--alternate-annotations", dest="prefer_alternate_annotation",
                       action="store_true",
                       help="Take the alternate source's coding features wherever it publishes "
                            "the same molecule and annotates at least as much of it, for "
                            "uniform gene models in place of each submission's own")
    _add_cache(viral)
    _add_force_timeout(viral, timeout=120)
    _add_json(viral)

    catalog = commands.add_parser("catalog", help="Download or build a registered ORF catalogue")
    catalog.add_argument("name")
    catalog.add_argument("--version")
    catalog.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH")
    _add_cache(catalog)
    _add_force_timeout(catalog)
    _add_json(catalog)

    dataset = commands.add_parser("dataset", help="Download or build a registered quantitative dataset")
    dataset.add_argument("name")
    dataset.add_argument("--version")
    dataset.add_argument("--source-path", action="append", default=[], metavar="ASSET=PATH")
    dataset.add_argument("--quiet", action="store_true")
    _add_cache(dataset)
    _add_force_timeout(dataset)
    _add_json(dataset)

    reference = commands.add_parser("reference", help="Download a dataset's native protein search FASTA")
    reference.add_argument("dataset")
    reference.add_argument("--version")
    reference.add_argument("--reference", help="Exact native reference context key")
    reference.add_argument("--asset", help="Declared source asset ID")
    reference.add_argument("--source-path", metavar="FILE")
    reference.add_argument("--quiet", action="store_true")
    _add_cache(reference)
    _add_force_timeout(reference)
    _add_json(reference)

    genome = commands.add_parser("genome", help="Download and index a pinned PyEnsembl annotation")
    genome.add_argument("name")
    genome.add_argument("--source-path", action="append", default=[], metavar="ROLE=PATH")
    _add_cache(genome)
    genome.add_argument("--force", action="store_true")
    _add_json(genome)

    chains = commands.add_parser("chains", help="Download pinned assembly chain files")
    chains.add_argument("source_build", metavar="FROM", nargs="?")
    chains.add_argument("target_build", metavar="TO", nargs="?")
    chains.add_argument("--species")
    chains.add_argument("--accept-terms", action="store_true")
    chains.add_argument("--source-path", metavar="FILE")
    _add_cache(chains)
    chains.add_argument("--force", action="store_true")
    _add_json(chains)

    args = parser.parse_args(argv)
    try:
        if args.command in ("list", "status"):
            rows = list_downloads(
                kind=args.kind, species=args.species, contexts=args.context,
                context_match=args.context_match, installed=args.installed,
                cache_root=args.cache_root, verify_files=args.verify)
            _show(rows, args.json)
            return
        if args.command == "viral":
            from .viral_panel import download_viral_references, index_viral_references

            bundle = download_viral_references(
                species=args.species, collection=args.collection, cache_root=args.cache_root,
                force=args.force, timeout=args.timeout,
                prefer_reference_molecule=args.prefer_reference_molecule,
                annotate_from_alternate=args.annotate_from_alternate,
                prefer_alternate_annotation=args.prefer_alternate_annotation)
            if args.context or args.exclude_context or args.subset:
                indexed = index_viral_references(
                    contexts=args.context, exclude_contexts=args.exclude_context,
                    context_match=args.context_match, subsets=args.subset,
                    species=args.species, collection=args.collection,
                    cache_root=args.cache_root, force=args.force, timeout=args.timeout)
                output = dict(indexed.to_dict(), kind="viral-index", status="available")
            else:
                output = dict(bundle.to_dict(), kind="viral", status="available")
        elif args.command == "catalog":
            from .catalogs import CatalogRegistry

            registry = CatalogRegistry(cache_root=args.cache_root)
            paths = _source_paths(args.source_path, "Source overrides")
            if registry.published(args.name, args.version):
                registry.download(args.name, args.version, force=args.force, source_paths=paths,
                                  timeout=args.timeout)
                inspection = registry.inspect(args.name, args.version)
                acquired = "downloaded"
            else:
                registry.build(args.name, args.version, force=args.force, source_paths=paths,
                               download_options={"timeout": args.timeout})
                inspection = registry.inspect_build(args.name, args.version)
                acquired = "built"
            output = dict(kind="catalog", name=args.name, version=registry.info(args.name, args.version).version,
                          acquired=acquired, status=inspection.status, path=str(inspection.path))
        elif args.command == "dataset":
            from .datasets import DatasetRegistry

            registry = DatasetRegistry(cache_root=args.cache_root)
            paths = _source_paths(args.source_path, "Source overrides")
            options = {"timeout": args.timeout, "show_progress": not args.quiet}
            if registry.published(args.name, args.version):
                registry.download(args.name, args.version, force=args.force, source_paths=paths, **options)
                inspection = registry.inspect(args.name, args.version)
                acquired = "downloaded"
            else:
                registry.build(args.name, args.version, force=args.force, source_paths=paths,
                               download_options=options)
                inspection = registry.inspect_build(args.name, args.version)
                acquired = "built"
            output = dict(kind="dataset", name=args.name, version=registry.info(args.name, args.version).version,
                          acquired=acquired, status=inspection.status, path=str(inspection.path))
        elif args.command == "reference":
            from .datasets import DatasetRegistry

            registry = DatasetRegistry(cache_root=args.cache_root)
            options = {"timeout": args.timeout, "show_progress": not args.quiet}
            paths = registry.download_reference(
                args.dataset, args.version, reference=args.reference, asset_id=args.asset,
                force=args.force, source_path=args.source_path, **options)
            inspection = registry.inspect_reference(
                args.dataset, args.version, reference=args.reference, asset_id=args.asset)
            output = dict(kind="reference", dataset=args.dataset,
                          version=registry.info(args.dataset, args.version).version,
                          status=inspection.status, path=str(inspection.path), files=paths)
        elif args.command == "genome":
            from .genomes import install_genome

            paths = _source_paths(args.source_path, "Source overrides")
            installed = install_genome(args.name, force=args.force, source_paths=paths,
                                       cache_root=args.cache_root)
            output = dict(kind="genome", name=args.name, reference_name=installed.reference_name,
                          annotation=f"{installed.annotation_name} {installed.annotation_version}",
                          transcripts=len(installed.transcript_ids()))
        else:
            from .liftover import download_chain, download_chains

            if args.source_build is None:
                if args.source_path:
                    raise ValueError("--source-path installs one chain; name it with FROM TO")
                result = download_chains(species=args.species, accept_terms=args.accept_terms,
                                         force=args.force, cache_root=args.cache_root)
                output = dict(downloaded={key: str(path) for key, path in result["downloaded"].items()},
                              skipped=result["skipped"])
            elif args.target_build is None:
                raise ValueError("Name both assemblies, FROM and TO, or neither")
            else:
                path = download_chain(
                    args.source_build, args.target_build, species=args.species,
                    accept_terms=args.accept_terms, force=args.force, source_path=args.source_path,
                    cache_root=args.cache_root)
                output = dict(path=str(path))
    except (ImportError, KeyError, RuntimeError, TypeError, ValueError, OSError) as error:
        parser.error(str(error))
    _show_result(output, args.json)
