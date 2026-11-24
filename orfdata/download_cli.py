"""Discover and acquire ORFdata resources through one command."""

import os
import sys
from contextlib import contextmanager, redirect_stdout

from .cli_output import manifest_summary, print_fields, print_json, print_table
from .cli_parser import UsageError, parse_source_paths
from .command_parsers import data_parser
from .downloads import (
    abbreviate_version,
    acquisition_action,
    format_byte_size,
    format_species_label,
    list_data,
    needs_replacement,
    normalize_status,
    refuse_unreadable,
)


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
    headings = ("KIND", "NAME", "SPECIES", "VERSION", "SIZE", "STATUS", "ACQUIRE")
    preview, hidden = _preview(rows, limit)
    values = [
        (
            row.kind,
            row.name,
            format_species_label(row.species),
            abbreviate_version(row.version),
            format_byte_size(row.download_size),
            row.status,
            row.acquire,
        )
        for row in preview
    ]
    print_table(
        headings,
        values,
        empty="No matching resources",
        limit=0,
        total=len(rows),
        hidden=hidden or None,
        nowrap=("KIND", "NAME", "VERSION", "SIZE", "STATUS", "ACQUIRE"),
    )
    _repair_hints(rows)
    print("\nDetails: orfdata data info NAME\nInstall: orfdata data fetch NAME\n"
          "Local path: orfdata data path NAME")
    print("Filter with --species, --kind or --all-versions; use --json for paths and source metadata.")


def _repair_hints(rows):
    """The exact command for each damaged or unreadable row, which a bare NAME can't give."""
    damaged = [row for row in rows if row.command and row.status == "damaged"]
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
            if result.get("download_size") is not None:
                summary["download_size"] = format_byte_size(result["download_size"])
            if result.get("installed_size") is not None:
                summary["installed_size"] = format_byte_size(result["installed_size"])
            summary["distribution"] = {
                key: value for key, value in result.get("distribution", {}).items()
                if key in ("access_terms", "redistributable")
            }
            manifest = manifest_summary(result["manifest"])
            summary.update((key, manifest[key]) for key in ("references", "source_assets", "coverage"))
            print_fields(summary)
            print("\nUse --json for complete source pins and metadata.")
            return
        if "acquired" in result:  # What a fetch did, whatever the kind.
            result = {key: result[key] for key in _FETCH_FIELDS if result.get(key) not in (None, (), [], False)}
        else:
            display = dict(result)
            if display.get("download_size") is not None:
                display["download_size"] = format_byte_size(display["download_size"])
            if display.get("installed_size") is not None:
                display["installed_size"] = format_byte_size(display["installed_size"])
            display.pop("datacache_status", None)
            result = display
        print_fields(result)


# What every fetch summary shows; empty contexts and an unrestricted resource are left out.
_FETCH_FIELDS = ("kind", "name", "version", "title", "species", "contexts", "acquired", "status", "path",
                 "restricted")
# What `data viral` shows of the collection or selection it installed.
_COLLECTION_FIELDS = ("kind", "status", "collection", "path", "bundle_version", "genome_count",
                      "transcript_count", "protein_count", "target_species", "contexts", "subsets")


@contextmanager
def _offline_override(offline):
    if not offline:
        yield
        return
    previous = os.environ.get("DATACACHE_OFFLINE")
    os.environ["DATACACHE_OFFLINE"] = "1"
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("DATACACHE_OFFLINE", None)
        else:
            os.environ["DATACACHE_OFFLINE"] = previous


def _dir_stats(path):
    if not os.path.exists(path):
        return False, 0, 0
    total_bytes = 0
    file_count = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        for name in files:
            try:
                st = os.lstat(os.path.join(root, name))
            except OSError:
                continue
            total_bytes += st.st_size
            file_count += 1
    return True, total_bytes, file_count


def _show_dirs(cache_root, json_output):
    import datacache

    from .cache import resolve_cache_location

    loc = resolve_cache_location(override=cache_root)
    exists, total_bytes, file_count = _dir_stats(loc.path)
    offline = datacache.offline_source()
    entries = []
    if os.path.isdir(loc.path):
        try:
            children = sorted(os.scandir(loc.path), key=lambda entry: entry.name)
        except OSError:
            children = []
        for entry in children:
            if not entry.is_dir(follow_symlinks=False):
                continue
            _, sub_bytes, sub_files = _dir_stats(entry.path)
            entries.append(dict(
                name=entry.name,
                path=entry.path,
                size_bytes=sub_bytes,
                file_count=sub_files,
            ))
    payload = dict(
        cache_root=loc.path,
        origin=loc.origin,
        detail=loc.detail,
        exists=exists,
        size_bytes=total_bytes,
        file_count=file_count,
        offline=offline is not None,
        offline_source=offline,
        subdirectories=entries,
    )
    if json_output:
        print_json(payload)
        return
    fields = {
        "cache_root": loc.path,
        "origin": f"{loc.detail} ({loc.origin})",
        "size": f"{datacache.format_bytes(total_bytes)} ({file_count} file(s))" if exists else "not created yet",
    }
    if offline:
        fields["offline"] = f"enabled ({offline})"
    print_fields(fields)
    if entries:
        print()
        print_table(
            ("AREA", "SIZE", "FILES", "PATH"),
            [(e["name"], datacache.format_bytes(e["size_bytes"]), str(e["file_count"]), e["path"])
             for e in entries],
            limit=0,
        )


def _show_prune(cache_root, keep, apply, json_output):
    import datacache
    from datacache.cli import prune_root

    from .cache import get_data_dir

    data_dir = str(get_data_dir(override=cache_root))
    results = prune_root(data_dir, keep=keep, apply=apply)
    removed_count = sum(r["removed_count"] for r in results)
    reclaimable = sum(r["reclaimable_bytes"] or 0 for r in results)
    payload = dict(
        cache_root=data_dir,
        dry_run=not apply,
        keep=keep,
        removed_count=removed_count,
        reclaimable_bytes=reclaimable,
        stores=results,
    )
    if json_output:
        print_json(payload)
        return
    if removed_count == 0:
        print(f"Nothing to prune under {data_dir}.")
        return
    verb = "Removed" if apply else "Would remove"
    for r in results:
        if r["removed_count"]:
            size_text = datacache.format_bytes(r["reclaimable_bytes"] or 0)
            print(f"{verb} {r['removed_count']} stale generation(s) ({size_text}) from {r['store']}")
    if not apply:
        print("\nDry run only. Pass --yes to delete these stale generations.")


def run(argv=None, prog="orfdata download"):
    parser = data_parser(prog)
    args = parser.parse_args(argv)
    try:
        with _offline_override(getattr(args, "offline", False)):
            if args.command == "dirs":
                _show_dirs(args.cache_root, args.json)
                return
            if args.command == "prune":
                _show_prune(args.cache_root, args.keep, args.yes, args.json)
                return
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
                    paths = parse_source_paths(args.source_path)
                    # Check every name, restricted terms, source overrides and permissions before downloading.
                    for name in names:
                        check_fetch(name, **options, force=args.force, source_paths=paths,
                                    accept_terms=args.accept_terms)
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
                    all_versions=args.all_versions,
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
                    paths = parse_source_paths(args.source_path)
                    before = normalize_status(
                        registry.inspect(args.name, args.version, verify_files=True).status
                    )
                    refuse_unreadable(before, force=args.force, name=args.name)
                    replace_existing = needs_replacement(before, force=args.force)
                    published = registry.published(args.name, args.version)
                    acquire = "download" if published else "build"
                    if published:
                        registry.download(args.name, args.version, force=replace_existing, source_paths=paths,
                                          timeout=args.timeout)
                        inspection = registry.inspect(args.name, args.version)
                    else:
                        registry.build(args.name, args.version, force=replace_existing, source_paths=paths,
                                       download_options={"timeout": args.timeout})
                        inspection = registry.inspect_build(args.name, args.version)
                    manifest = registry.info(args.name, args.version)
                    output = dict(
                        kind="catalog",
                        name=args.name,
                        version=manifest.version,
                        title=manifest.title,
                        acquired=acquisition_action(before, acquire=acquire, force=args.force),
                        status=normalize_status(inspection.status),
                        datacache_status=inspection.status,
                        path=str(inspection.path),
                    )
                elif args.command == "dataset":
                    from .datasets import DatasetRegistry

                    registry = DatasetRegistry(cache_root=args.cache_root)
                    paths = parse_source_paths(args.source_path)
                    options = {"timeout": args.timeout, "show_progress": not args.quiet}
                    before = normalize_status(
                        registry.inspect(args.name, args.version, verify_files=True).status
                    )
                    refuse_unreadable(before, force=args.force, name=args.name)
                    replace_existing = needs_replacement(before, force=args.force)
                    published = registry.published(args.name, args.version)
                    acquire = "download" if published else "build"
                    if published:
                        registry.download(args.name, args.version, force=replace_existing, source_paths=paths, **options)
                        inspection = registry.inspect(args.name, args.version)
                    else:
                        registry.build(args.name, args.version, force=replace_existing, source_paths=paths,
                                       download_options=options)
                        inspection = registry.inspect_build(args.name, args.version)
                    manifest = registry.info(args.name, args.version)
                    output = dict(
                        kind="dataset",
                        name=args.name,
                        version=manifest.version,
                        title=manifest.title,
                        acquired=acquisition_action(before, acquire=acquire, force=args.force),
                        status=normalize_status(inspection.status),
                        datacache_status=inspection.status,
                        path=str(inspection.path),
                    )
                elif args.command == "reference":
                    from .dataset_cli import validate_search_reference_options
                    from .datasets import DatasetRegistry

                    registry = DatasetRegistry(cache_root=args.cache_root)
                    validate_search_reference_options(
                        registry, args.dataset, args.version, reference=args.reference, asset_id=args.asset,
                    )
                    options = {"timeout": args.timeout, "show_progress": not args.quiet}
                    paths = registry.download_reference(
                        args.dataset, args.version, reference=args.reference, asset_id=args.asset,
                        force=args.force, source_path=args.source_path, **options)
                    inspection = registry.inspect_reference(
                        args.dataset, args.version, reference=args.reference, asset_id=args.asset)
                    output = dict(
                        kind="reference",
                        dataset=args.dataset,
                        version=registry.info(args.dataset, args.version).version,
                        status=normalize_status(inspection.status),
                        datacache_status=inspection.status,
                        path=str(inspection.path),
                        files=paths,
                    )
                elif args.command == "genome":
                    from .genomes import install_genome

                    paths = parse_source_paths(args.source_path)
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
                        path = download_chain(
                            args.source_build, args.target_build, species=args.species,
                            accept_terms=args.accept_terms, force=args.force, source_path=args.source_path,
                            cache_root=args.cache_root, download_options=dict(timeout=args.timeout))
                        output = dict(path=str(path))
    except (ImportError, RuntimeError, ValueError, OSError) as error:
        parser.fail(error)
    _show_result(output, args.json)
