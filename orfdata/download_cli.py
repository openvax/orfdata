"""Discover and acquire ORFdata resources through one command."""

import os
import sys
from contextlib import contextmanager, redirect_stdout

from .cli_output import manifest_summary, print_fields, print_json, print_table
from .cli_parser import UsageError, parse_source_paths
from .command_parsers import data_parser
from .downloads import (
    abbreviate_version,
    format_byte_size,
    format_species_label,
    list_data,
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


def _show_prune(cache_root, keep, apply, sources, json_output):
    import datacache
    from datacache.cli import prune_root

    from .cache import get_data_dir
    from .source_store import prune_source_store

    data_dir = str(get_data_dir(override=cache_root))
    results = prune_root(data_dir, keep=keep, apply=apply)
    removed_count = sum(r["removed_count"] for r in results)
    reclaimable = sum(r["reclaimable_bytes"] or 0 for r in results)
    source_summary = prune_source_store(cache_root=cache_root, apply=apply) if sources else None
    payload = dict(
        cache_root=data_dir,
        dry_run=not apply,
        keep=keep,
        removed_count=removed_count,
        reclaimable_bytes=reclaimable,
        stores=results,
    )
    if source_summary is not None:
        payload["sources"] = source_summary
    if json_output:
        print_json(payload)
        return
    source_removed = 0 if source_summary is None else source_summary["removed_count"]
    if removed_count == 0 and source_removed == 0:
        print(f"Nothing to prune under {data_dir}.")
        return
    verb = "Removed" if apply else "Would remove"
    for r in results:
        if r["removed_count"]:
            size_text = datacache.format_bytes(r["reclaimable_bytes"] or 0)
            print(f"{verb} {r['removed_count']} stale generation(s) ({size_text}) from {r['store']}")
    if source_summary is not None and source_removed:
        size_text = datacache.format_bytes(source_summary["reclaimable_bytes"])
        print(f"{verb} {source_removed} cached build source(s) ({size_text}) from {source_summary['path']}")
    if not apply:
        print("\nDry run only. Pass --yes to delete these stale generations.")


def _fetch_chains(args):
    """Compatibility helper shared by ``data chains`` and ``chains download``."""
    from .data_resources import fetch_data
    from .liftover import find_pinned_chain, pinned_chains

    paths = parse_source_paths(args.source_path)
    if args.source_build is None:
        if paths:
            raise UsageError("--source-path installs one chain; name it with FROM TO")
        pins = pinned_chains(species=args.species)
        if not pins:
            raise ValueError(
                f"No chain files are pinned for {args.species}; "
                f"pinned species: {', '.join(sorted({c.species for c in pinned_chains()}))}"
            )
        downloaded, skipped = {}, {}
        for chain in pins:
            if chain.restricted and not args.accept_terms:
                skipped[chain.name] = (
                    f"{chain.provider}'s terms need accept_terms=True (--accept-terms): {chain.terms}"
                )
                continue
            result = fetch_data(
                chain.name,
                kind="chain",
                cache_root=args.cache_root,
                force=args.force,
                timeout=args.timeout,
                accept_terms=args.accept_terms,
                resume=args.resume,
            )
            downloaded[chain.name] = result["path"]
        return dict(downloaded=downloaded, skipped=skipped)
    if args.target_build is None:
        raise UsageError("Name both assemblies, FROM and TO, or neither")
    chain = find_pinned_chain(args.source_build, args.target_build, species=args.species)
    result = fetch_data(
        chain.name,
        kind="chain",
        cache_root=args.cache_root,
        force=args.force,
        source_paths=paths,
        timeout=args.timeout,
        accept_terms=args.accept_terms,
        resume=args.resume,
    )
    from .liftover import chain_registry

    try:
        file_path = str(
            chain_registry(chain, args.cache_root).local_path(
                chain.name, "1", asset=chain.filename, verify_files=False,
            )
        )
    except Exception:
        file_path = result["path"]
    return dict(result, path=file_path)


def _fetch_reference(args, *, registry=None):
    """Compatibility helper shared by ``data reference`` and ``datasets download-reference``."""
    from .data_resources import fetch_data
    from .dataset_cli import validate_search_reference_options
    from .datasets import DatasetRegistry

    reg = DatasetRegistry(cache_root=args.cache_root) if registry is None else registry
    dataset_name = getattr(args, "dataset", None) or args.name
    validate_search_reference_options(
        reg, dataset_name, args.version, reference=args.reference, asset_id=args.asset,
    )
    specs = [
        spec for spec in reg.search_references(dataset_name, args.version)
        if (args.reference is None or spec.reference.key == args.reference)
        and (args.asset is None or spec.source.asset_id == args.asset)
    ]
    if not specs:
        raise ValueError("No declared native search reference matches this dataset/context/asset")
    if len(specs) != 1:
        raise ValueError("Select an exact reference context and source asset")
    spec = specs[0]
    paths = parse_source_paths(args.source_path)
    return fetch_data(
        f"{dataset_name}:{spec.source.asset_id}",
        kind="reference",
        version=args.version,
        cache_root=args.cache_root,
        force=args.force,
        source_paths=paths,
        timeout=args.timeout,
        show_progress=not args.quiet,
        resume=args.resume,
        options={"registry": reg} if registry is not None else None,
    )


def run(argv=None, prog="orfdata download"):
    parser = data_parser(prog)
    args = parser.parse_args(argv)
    try:
        with _offline_override(getattr(args, "offline", False)):
            if args.command == "dirs":
                _show_dirs(args.cache_root, args.json)
                return
            if args.command == "prune":
                _show_prune(args.cache_root, args.keep, args.yes, args.sources, args.json)
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
                    fetch_extras = dict(
                        genome_fasta_path_or_url=args.genome_fasta,
                        prefer_reference_molecule=args.prefer_reference_molecule,
                        annotate_from_alternate=args.annotate_from_alternate,
                        prefer_alternate_annotation=args.prefer_alternate_annotation,
                    )
                    # Check every name, restricted terms, source overrides and permissions before downloading.
                    for name in names:
                        check_fetch(
                            name,
                            **options,
                            force=args.force,
                            source_paths=paths,
                            accept_terms=args.accept_terms,
                            **fetch_extras,
                        )
                    # Source processors may print progress; keep stdout machine-readable.
                    with redirect_stdout(sys.stderr):
                        output = [
                            fetch_data(
                                name,
                                **options,
                                force=args.force,
                                source_paths=paths,
                                timeout=args.timeout,
                                show_progress=not args.quiet and not args.json,
                                accept_terms=args.accept_terms,
                                resume=args.resume,
                                **fetch_extras,
                            )
                            for name in names
                        ]
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

            from .data_resources import fetch_data

            with redirect_stdout(sys.stderr):
                if args.command == "viral":
                    from .viral_panel import open_viral_references, select_viral_references

                    paths = parse_source_paths(args.source_path)
                    fetched = fetch_data(
                        args.collection,
                        kind="viral",
                        cache_root=args.cache_root,
                        force=args.force,
                        source_paths=paths,
                        timeout=args.timeout,
                        resume=args.resume,
                        prefer_reference_molecule=args.prefer_reference_molecule,
                        annotate_from_alternate=args.annotate_from_alternate,
                        prefer_alternate_annotation=args.prefer_alternate_annotation,
                        options={"species": args.species},
                    )
                    if args.context or args.exclude_context or args.subset:
                        indexed = select_viral_references(
                            contexts=args.context, exclude_contexts=args.exclude_context,
                            context_match=args.context_match, subsets=args.subset,
                            species=args.species, collection=args.collection,
                            cache_root=args.cache_root, force=args.force, timeout=args.timeout)
                        output = dict(indexed.to_dict(), kind="viral-selection", status="available")
                    else:
                        try:
                            bundle = open_viral_references(
                                species=args.species,
                                collection=args.collection,
                                cache_root=args.cache_root,
                            )
                            output = dict(bundle.to_dict(), kind="viral", status=fetched["status"])
                        except Exception:
                            output = fetched
                    if not args.json and "collection" in output:
                        output = {key: value for key, value in output.items() if key in _COLLECTION_FIELDS}
                elif args.command == "catalog":
                    output = fetch_data(
                        args.name,
                        kind="catalog",
                        version=args.version,
                        cache_root=args.cache_root,
                        force=args.force,
                        source_paths=parse_source_paths(args.source_path),
                        timeout=args.timeout,
                        resume=args.resume,
                    )
                elif args.command == "dataset":
                    output = fetch_data(
                        args.name,
                        kind="dataset",
                        version=args.version,
                        cache_root=args.cache_root,
                        force=args.force,
                        source_paths=parse_source_paths(args.source_path),
                        timeout=args.timeout,
                        show_progress=not args.quiet,
                        resume=args.resume,
                    )
                elif args.command == "reference":
                    output = _fetch_reference(args)
                elif args.command == "genome":
                    output = fetch_data(
                        args.name,
                        kind="genome",
                        cache_root=args.cache_root,
                        force=args.force,
                        source_paths=parse_source_paths(args.source_path),
                        timeout=args.timeout,
                        resume=args.resume,
                        genome_fasta_path_or_url=args.genome_fasta,
                    )
                else:
                    output = _fetch_chains(args)
    except (ImportError, RuntimeError, ValueError, OSError) as error:
        parser.fail(error)
    _show_result(output, args.json)
