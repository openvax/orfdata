"""Discover and acquire ORFdata resources through one command."""

import sys
from contextlib import redirect_stdout

from .cli_output import manifest_summary, print_fields, print_json, print_table
from .cli_parser import UsageError
from .command_parsers import data_parser
from .downloads import list_data


def _source_paths(items, label):
    paths = {}
    for item in items:
        key, separator, path = item.partition("=")
        if not separator or not key or not path or key in paths:
            raise UsageError(f"{label} require unique NAME=PATH entries")
        paths[key] = path
    return paths or None


def _preview(rows, limit):
    """At most ``limit`` rows, taken from each kind in turn, in their listed order; also what each kind hides.

    Every kind appears whenever ``limit`` leaves room for one row of each.
    """
    if limit < 0:
        raise ValueError("limit must be nonnegative")
    if limit == 0 or len(rows) <= limit:
        return list(rows), {}
    counts, rank = {}, []
    for row in rows:  # A row's place within its kind: the first of each kind comes before any second.
        rank.append(counts.get(row.kind, 0))
        counts[row.kind] = rank[-1] + 1
    kinds = list(counts)
    keep = set(sorted(range(len(rows)), key=lambda i: (rank[i], kinds.index(rows[i].kind)))[:limit])
    shown = [row for index, row in enumerate(rows) if index in keep]
    hidden = {kind: counts[kind] - sum(row.kind == kind for row in shown) for kind in kinds}
    return shown, {kind: count for kind, count in hidden.items() if count}


def _show(rows, json_output, limit=25):
    if json_output:
        print_json([row.to_dict() for row in rows])
        return
    headings = ("KIND", "NAME", "VERSION", "STATUS", "ACQUIRE")
    preview, hidden = _preview(rows, limit)
    values = [(row.kind, row.name, row.version, row.status, row.acquire) for row in preview]
    print_table(headings, values, empty="No matching resources", limit=0)
    if hidden:
        more = ", ".join(f"{count} {kind}" for kind, count in hidden.items())
        print(f"\nShowing {len(preview)} of {len(rows)} entries; not shown: {more}. Use --limit 0 for all rows.")
    _repair_hints(rows)
    print("\nDetails: orfdata data info NAME\nInstall: orfdata data fetch NAME\n"
          "Local path: orfdata data path NAME")
    print("Filter with --species or --kind; use --json for paths and source metadata.")


def _repair_hints(rows):
    """The exact command for each damaged or unreadable row, which a bare NAME can't give."""
    damaged = [row for row in rows if row.command and row.status not in ("available", "missing", "inaccessible")]
    unreadable = [row for row in rows if row.command and row.status == "inaccessible"]
    for found, heading, suffix in ((damaged, "Damaged; 'data fetch' replaces them:", ""),
                                   (unreadable, "Unreadable; fix their permissions, or replace them:", " --force")):
        if found:
            print(f"\n{heading}")
            for row in found[:5]:
                print(f"  {row.command}{suffix}")
            if len(found) > 5:
                print(f"  ... {len(found) - 5} more; see --json")


def _show_result(result, json_output):
    if json_output:
        print_json(result)
    else:
        if "manifest" in result:
            summary = {key: result[key] for key in
                       ("kind", "name", "version", "title", "status", "acquired", "acquire", "installed",
                        "species", "contexts", "path", "restricted", "command") if key in result}
            summary["distribution"] = {key: value for key, value in result.get("distribution", {}).items()
                                       if key in ("access_terms", "redistributable", "download_size", "installed_size")}
            manifest = manifest_summary(result["manifest"])
            summary.update((key, manifest[key]) for key in ("references", "source_assets", "coverage"))
            print_fields(summary)
            print("\nUse --json for complete source pins and metadata.")
            return
        if "acquired" in result:  # What a fetch did, whatever the kind.
            result = {key: result[key] for key in _FETCH_FIELDS if result.get(key) not in (None, (), [], False)}
        print_fields(result)


# What every fetch summary shows; empty contexts and an unrestricted resource are left out.
_FETCH_FIELDS = ("kind", "name", "version", "title", "species", "contexts", "acquired", "status", "path",
                 "restricted")
# What `data viral` shows of the collection or selection it installed.
_COLLECTION_FIELDS = ("kind", "status", "collection", "path", "bundle_version", "genome_count",
                      "transcript_count", "protein_count", "target_species", "contexts", "subsets")


def run(argv=None, prog="orfdata download"):
    parser = data_parser(prog)
    args = parser.parse_args(argv)
    try:
        if args.command in ("info", "path", "fetch"):
            from .data_resources import check_fetch, data_info, data_path, fetch_data

            options = dict(kind=args.kind, version=args.version, cache_root=args.cache_root)
            if args.command == "info":
                _show_result(data_info(args.names, verify_files=args.verify, **options), args.json)
            elif args.command == "path":
                path = str(data_path(args.names, verify_files=args.verify, **options))
                print_json(dict(path=path)) if args.json else print(path)
            else:
                names = list(dict.fromkeys(args.names))
                if len(names) > 1 and (args.source_path or args.version):
                    raise UsageError("--source-path and --version require exactly one resource")
                # Terms are a command-line omission, found without verifying anything.
                for name in names:
                    info = data_info(name, **options)
                    if info["restricted"] and not args.accept_terms:
                        raise UsageError(f"{name!r} requires --accept-terms; see 'orfdata data info {name}'")
                paths = _source_paths(args.source_path, "Source overrides")
                if len(names) > 1:  # Find every name, terms and permissions problem before downloading anything.
                    for name in names:
                        check_fetch(name, **options, force=args.force, accept_terms=args.accept_terms)
                # Source processors may print progress; keep stdout machine-readable.
                with redirect_stdout(sys.stderr):
                    output = [fetch_data(name, **options, force=args.force, source_paths=paths,
                                         timeout=args.timeout, show_progress=not args.quiet and not args.json,
                                         accept_terms=args.accept_terms, resume=args.resume) for name in names]
                if args.json:
                    print_json(output)
                else:
                    for index, result in enumerate(output):
                        if index:
                            print()
                        _show_result(result, False)
            return
        if args.command in ("list", "status"):
            rows = list_data(
                kind=args.kind, species=args.species, contexts=args.context,
                context_match=args.context_match, installed=args.installed, local=args.local,
                cache_root=args.cache_root, verify_files=args.verify)
            _show(rows, args.json, args.limit)
            return
        with redirect_stdout(sys.stderr):
            if args.command == "viral":
                from .viral_panel import download_viral_references, select_viral_references

                bundle = download_viral_references(
                    species=args.species, collection=args.collection, cache_root=args.cache_root,
                    force=args.force, timeout=args.timeout,
                    prefer_reference_molecule=args.prefer_reference_molecule,
                    annotate_from_alternate=args.annotate_from_alternate,
                    prefer_alternate_annotation=args.prefer_alternate_annotation)
                if args.context or args.exclude_context or args.subset:
                    indexed = select_viral_references(
                        contexts=args.context, exclude_contexts=args.exclude_context,
                        context_match=args.context_match, subsets=args.subset,
                        species=args.species, collection=args.collection,
                        cache_root=args.cache_root, force=args.force, timeout=args.timeout)
                    output = dict(indexed.to_dict(), kind="viral-selection", status="available")
                else:
                    output = dict(bundle.to_dict(), kind="viral", status="available")
                if not args.json:
                    output = {key: value for key, value in output.items() if key in _COLLECTION_FIELDS}
            elif args.command == "catalog":
                from .catalogs import CatalogRegistry

                registry = CatalogRegistry(cache_root=args.cache_root)
                paths = _source_paths(args.source_path, "Source overrides")
                reused = registry.inspect(args.name, args.version, verify_files=False).status == "available"
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
                manifest = registry.info(args.name, args.version)
                output = dict(kind="catalog", name=args.name, version=manifest.version, title=manifest.title,
                              acquired="already installed" if reused and not args.force else acquired,
                              status=inspection.status, path=str(inspection.path))
            elif args.command == "dataset":
                from .datasets import DatasetRegistry

                registry = DatasetRegistry(cache_root=args.cache_root)
                paths = _source_paths(args.source_path, "Source overrides")
                options = {"timeout": args.timeout, "show_progress": not args.quiet}
                reused = registry.inspect(args.name, args.version, verify_files=False).status == "available"
                if registry.published(args.name, args.version):
                    registry.download(args.name, args.version, force=args.force, source_paths=paths, **options)
                    inspection = registry.inspect(args.name, args.version)
                    acquired = "downloaded"
                else:
                    registry.build(args.name, args.version, force=args.force, source_paths=paths,
                                   download_options=options)
                    inspection = registry.inspect_build(args.name, args.version)
                    acquired = "built"
                manifest = registry.info(args.name, args.version)
                output = dict(kind="dataset", name=args.name, version=manifest.version, title=manifest.title,
                              acquired="already installed" if reused and not args.force else acquired,
                              status=inspection.status, path=str(inspection.path))
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
                                           cache_root=args.cache_root,
                                           download_options=dict(timeout=args.timeout))
                output = dict(kind="genome", name=args.name, reference_name=installed.reference_name,
                              annotation=f"{installed.annotation_name} {installed.annotation_version}",
                              transcripts=len(installed.transcript_ids()))
            else:
                from .liftover import download_chain, download_chains

                if args.source_build is None:
                    if args.source_path:
                        raise UsageError("--source-path installs one chain; name it with FROM TO")
                    result = download_chains(species=args.species, accept_terms=args.accept_terms,
                                             force=args.force, cache_root=args.cache_root,
                                             download_options=dict(timeout=args.timeout))
                    output = dict(downloaded={key: str(path) for key, path in result["downloaded"].items()},
                                  skipped=result["skipped"])
                elif args.target_build is None:
                    raise UsageError("Name both assemblies, FROM and TO, or neither")
                else:
                    from .liftover_cli import require_terms

                    require_terms(args.source_build, args.target_build, args.species, args.accept_terms)
                    path = download_chain(
                        args.source_build, args.target_build, species=args.species,
                        accept_terms=args.accept_terms, force=args.force, source_path=args.source_path,
                        cache_root=args.cache_root, download_options=dict(timeout=args.timeout))
                    output = dict(path=str(path))
    except (ImportError, RuntimeError, ValueError, OSError) as error:
        parser.fail(error)
    _show_result(output, args.json)
