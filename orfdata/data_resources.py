"""Named access to optional data, independent of source and species."""

from dataclasses import replace
from pathlib import Path

from .cache import DOWNLOAD_TIMEOUT
from .downloads import (
    _fetch_command,
    acquisition_action,
    list_data,
    needs_replacement,
    refuse_unreadable,
)
from .errors import NotInstalledError, TermsNotAccepted
from .resource_kinds import get_resource_kind


def _kind_options(
    *,
    genome_fasta_path_or_url=None,
    prefer_reference_molecule=True,
    annotate_from_alternate=True,
    prefer_alternate_annotation=False,
    extra_options=None,
):
    """Collect any non-default kind-specific acquisition options into a dict."""
    options = dict(extra_options or {})
    if genome_fasta_path_or_url is not None:
        options["genome_fasta_path_or_url"] = genome_fasta_path_or_url
    if prefer_reference_molecule is not True:
        options["prefer_reference_molecule"] = prefer_reference_molecule
    if annotate_from_alternate is not True:
        options["annotate_from_alternate"] = annotate_from_alternate
    if prefer_alternate_annotation is not False:
        options["prefer_alternate_annotation"] = prefer_alternate_annotation
    return options


def _resolve_kind(name, kind, cache_root, options=None):
    """Determine the single registered kind owning ``name`` without per-kind branching."""
    if kind is not None:
        get_resource_kind(kind)  # Validates kind name immediately.
        if options and options.get("registry") is not None:
            return kind
    rows = [
        row for row in list_data(kind=kind, cache_root=cache_root, verify_files=False, all_versions=True)
        if row.name == name
    ]
    if rows:
        kinds = {row.kind for row in rows}
        if len(kinds) != 1:
            raise ValueError(f"Ambiguous data resource {name!r}; select --kind")
        return next(iter(kinds))
    if kind is not None:
        # Allow kind-specific alias resolution (e.g. genome aliases such as 'hs1' or custom names).
        try:
            get_resource_kind(kind).inspect(name, None, cache_root=cache_root, verify_files=False, options=options)
            return kind
        except ValueError:
            pass
    raise ValueError(f"Unknown data resource {name!r}; see 'orfdata data list'")


def _resolve(name, kind, version, cache_root, verify_files, options=None):
    resolved_kind = _resolve_kind(name, kind, cache_root, options=options)
    handler = get_resource_kind(resolved_kind)
    row = handler.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
    return replace(row, command=_fetch_command(row, cache_root))


def data_info(name, *, kind=None, version=None, cache_root=None, verify_files=False):
    """Describe a named resource, its pins and acquisition route, entirely offline."""
    resolved_kind = _resolve_kind(name, kind, cache_root)
    handler = get_resource_kind(resolved_kind)
    row = handler.inspect(name, version, cache_root=cache_root, verify_files=verify_files)
    result = handler.info(name, row.version, cache_root=cache_root, verify_files=verify_files)
    result["command"] = _fetch_command(row, cache_root)
    return result


def data_path(name, *, kind=None, version=None, cache_root=None, verify_files=False):
    """Return an installed file/directory; missing data never triggers acquisition."""
    info = data_info(name, kind=kind, version=version, cache_root=cache_root,
                     verify_files=verify_files)
    if not info["installed"] or info["path"] is None:
        if info["status"] == "missing":
            raise NotInstalledError(f"{name!r} is not installed; install it with: {info['command']}")
        if info["status"] == "inaccessible":
            refuse_unreadable(info["status"], name=f"{name!r} at {info['path']}", command=info["command"])
        raise ValueError(f"{name!r} is {info['status']}; repair it with: {info['command']}")
    return Path(info["path"])


def check_fetch(
    name,
    *,
    kind=None,
    version=None,
    cache_root=None,
    force=False,
    source_paths=None,
    accept_terms=False,
    genome_fasta_path_or_url=None,
    prefer_reference_molecule=True,
    annotate_from_alternate=True,
    prefer_alternate_annotation=False,
    options=None,
):
    """Verify a resource and raise whatever would stop :func:`fetch_data`, before anything downloads.

    Returns the verified status row. Batch callers check every name first.
    Restricted resources raise :class:`~orfdata.errors.TermsNotAccepted` before
    file verification when ``accept_terms`` is false. Kind-specific options and
    ``--source-path`` overrides are validated against the resource's handler.
    """
    kind_opts = _kind_options(
        genome_fasta_path_or_url=genome_fasta_path_or_url,
        prefer_reference_molecule=prefer_reference_molecule,
        annotate_from_alternate=annotate_from_alternate,
        prefer_alternate_annotation=prefer_alternate_annotation,
        extra_options=options,
    )
    unverified = _resolve(name, kind, version, cache_root, False, options=kind_opts or None)
    if unverified.restricted and not accept_terms:
        raise TermsNotAccepted(f"{name!r} requires --accept-terms; see 'orfdata data info {name}'")
    handler = get_resource_kind(unverified.kind)
    handler.validate_options(
        unverified, cache_root=cache_root, source_paths=source_paths, options=kind_opts or None,
    )
    row = _resolve(name, unverified.kind, unverified.version, cache_root, True, options=kind_opts or None)
    if row.acquire is not None and not (row.installed and not force):
        refuse_unreadable(row.status, force=force, name=f"{name!r} at {row.path}", command=row.command)
    return row


def fetch_data(
    name,
    *,
    kind=None,
    version=None,
    cache_root=None,
    force=False,
    source_paths=None,
    timeout=DOWNLOAD_TIMEOUT,
    show_progress=True,
    accept_terms=False,
    resume=False,
    genome_fasta_path_or_url=None,
    prefer_reference_molecule=True,
    annotate_from_alternate=True,
    prefer_alternate_annotation=False,
    options=None,
):
    """Install one resource and return its status, downloading or building as declared.

    Downloads retain their checksum validation. Builds use the same pinned source
    adapters as the focused commands. Existing valid files are reused
    (``acquired="already installed"``), files that fail verification are
    replaced (``acquired="repaired"``), and ``force=True`` replaces valid files
    too (``acquired="replaced"``). Discovery and ordinary queries never
    invoke this function implicitly.

    ``source_paths`` accepts either a mapping of ``{asset: path}`` or a single
    bare path when the resource declares exactly one asset.

    ``resume=True`` keeps partial downloads so an interrupted one continues
    where it stopped (DataCache's ``resume="auto"``). Partials and staging then
    stay on disk after a failure, and an archive keeps every part until it is
    installed.
    """
    kind_opts = _kind_options(
        genome_fasta_path_or_url=genome_fasta_path_or_url,
        prefer_reference_molecule=prefer_reference_molecule,
        annotate_from_alternate=annotate_from_alternate,
        prefer_alternate_annotation=prefer_alternate_annotation,
        extra_options=options,
    )
    row = check_fetch(
        name,
        kind=kind,
        version=version,
        cache_root=cache_root,
        force=force,
        source_paths=source_paths,
        accept_terms=accept_terms,
        options=kind_opts or None,
    )
    needs_extra_action = bool(
        kind_opts.get("genome_fasta_path_or_url") is not None
        or kind_opts.get("prefer_reference_molecule") is False
        or kind_opts.get("annotate_from_alternate") is False
        or kind_opts.get("prefer_alternate_annotation") is True
    )
    if row.acquire is None or (row.installed and not force and not needs_extra_action):
        return dict(row.to_dict(), acquired="already installed")
    acquired = acquisition_action(row.status, acquire=row.acquire, force=force)
    replace_existing = needs_replacement(row.status, force=force)
    if type(resume) is not bool:
        raise ValueError("resume must be True or False")
    # Nothing to resume when local files stand in for the downloads.
    resumable = resume and not source_paths
    download_opts = dict(
        timeout=timeout,
        show_progress=show_progress,
        **({"resume": "auto"} if resumable else {}),
    )
    handler = get_resource_kind(row.kind)
    handler.acquire(
        row.name,
        row.version,
        cache_root=cache_root,
        force=replace_existing,
        source_paths=source_paths,
        accept_terms=accept_terms,
        download_options=download_opts,
        options=kind_opts or None,
    )
    result = _resolve(row.name, row.kind, row.version, cache_root, False, options=kind_opts or None)
    if not result.installed:
        raise RuntimeError(f"{name!r} installation finished with status {result.status}")
    return dict(result.to_dict(), acquired=acquired)
