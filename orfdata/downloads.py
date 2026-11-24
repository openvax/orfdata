"""Discover ORFdata resources and their local installation state."""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from os import PathLike
from shlex import join
from typing import Optional

from .cache import DOWNLOAD_TIMEOUT
from .reference import canonical_species, species_taxon

# Unified local resource states shared by data commands, catalogues, datasets and setup.
RESOURCE_STATES = ("available", "bundled", "missing", "damaged", "inaccessible")

# Unified acquisition actions shared by data fetch, focused data commands and setup.
ACQUISITION_ACTIONS = (
    "downloaded",
    "built",
    "repaired",
    "replaced",
    "already installed",
    "skipped",
    "failed",
)

_STATUS_MAP = {
    "available": "available",
    "bundled": "bundled",
    "missing": "missing",
    "invalid": "damaged",
    "corrupt": "damaged",
    "damaged": "damaged",
    "inaccessible": "inaccessible",
}

_SHORT_SPECIES = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "canis_lupus_familiaris": "dog",
    "danio_rerio": "zebrafish",
    "caenorhabditis_elegans": "worm",
    "saccharomyces_cerevisiae": "yeast",
}


def normalize_status(status: str, *, bundled: bool = False) -> str:
    """Map a local or DataCache inspection status to ORFdata's unified status vocabulary.

    ORFdata uses one documented set of local states across ``data`` commands,
    catalogues, datasets, viral collections and ``setup``: ``"available"``,
    ``"bundled"``, ``"missing"``, ``"damaged"``, and ``"inaccessible"``.
    DataCache's ``"invalid"`` and legacy ``"corrupt"`` both map to
    ``"damaged"``.

    Parameters
    ----------
    status : str
        Raw status string from DataCache or an ORFdata inspector.
    bundled : bool, optional
        When ``True``, returns ``"bundled"`` regardless of ``status``.

    Returns
    -------
    str
        One of :data:`RESOURCE_STATES`.

    Raises
    ------
    ValueError
        If ``status`` is not a recognized state.
    """
    if bundled:
        return "bundled"
    normalized = _STATUS_MAP.get(status)
    if normalized is None:
        raise ValueError(
            f"Unknown resource status {status!r}; expected one of {', '.join(RESOURCE_STATES)}"
        )
    return normalized


def acquisition_action(before_status: str, *, acquire: str | None, force: bool = False) -> str:
    """Return the unified acquisition action for installing over ``before_status``.

    Shared by :func:`~orfdata.data_resources.fetch_data`, focused ``data``
    subcommands, and :func:`~orfdata.install.install_species` so every command
    reports the same action vocabulary (:data:`ACQUISITION_ACTIONS`).

    Parameters
    ----------
    before_status : str
        Local status observed before acquisition (raw or normalized).
    acquire : str or None
        Declared acquisition route: ``"download"``, ``"build"``, or ``None``
        for bundled resources.
    force : bool, optional
        Whether forced replacement was requested.

    Returns
    -------
    str
        One of ``"already installed"``, ``"downloaded"``, ``"built"``,
        ``"repaired"``, or ``"replaced"``.
    """
    state = normalize_status(before_status)
    if acquire is None or (state in ("available", "bundled") and not force):
        return "already installed"
    if state in ("damaged", "inaccessible"):
        return "repaired"
    if state in ("available", "bundled") and force:
        return "replaced"
    return "downloaded" if acquire == "download" else "built"


def abbreviate_version(version: Optional[str], max_length: int = 14) -> str:
    """Shorten long version identifiers for human tables while leaving short ones intact.

    Parameters
    ----------
    version : str or None
        Full version string (such as a 67-character viral bundle hash).
    max_length : int, optional
        Maximum rendered length before truncating with ``"..."`` (default 14).

    Returns
    -------
    str
        Abbreviated version for table display, or ``"-"`` when ``version`` is
        empty or ``None``.
    """
    if not version:
        return "-"
    text = str(version)
    if len(text) <= max_length:
        return text
    if max_length <= 3:
        return text[:max_length]
    return text[: max_length - 3] + "..."


def format_species_label(species: Iterable[str]) -> str:
    """Format canonical species names compactly for human discovery tables.

    Parameters
    ----------
    species : iterable of str
        Canonical species names such as ``("homo_sapiens",)``.

    Returns
    -------
    str
        Comma-separated common names (e.g. ``"human"``, ``"mouse"``), or
        ``"-"`` when empty.
    """
    labels = [_SHORT_SPECIES.get(name, name) for name in species]
    return ", ".join(labels) if labels else "-"


def format_byte_size(size: Optional[int]) -> str:
    """Format an optional byte count with DataCache's human-readable units.

    Parameters
    ----------
    size : int or None
        Size in bytes, or ``None`` when unknown.

    Returns
    -------
    str
        Human-readable size string (e.g. ``"87.3 MiB"``) or ``"-"`` when
        ``size`` is ``None``.
    """
    if size is None:
        return "-"
    import datacache

    return datacache.format_bytes(size)


def release_sizes(registry, name, version, *, published=None):
    """Return ``(download_size, installed_size)`` in bytes for one registry release.

    Uses explicit ``distribution.download_size`` and
    ``distribution.installed_size`` when declared; otherwise derives download
    and installed sizes from pinned assets, archive parts, or build sources.

    Parameters
    ----------
    registry : ReleaseRegistry
        Catalogue or dataset registry owning ``name``.
    name : str
        Registered release name.
    version : str
        Resolved pinned release version.
    published : bool, optional
        Whether the release has a published artifact; queried from ``registry``
        when omitted.

    Returns
    -------
    tuple of (int or None, int or None)
        ``(download_size, installed_size)`` in bytes.
    """
    meta = registry.download_metadata(name, version)
    distribution = meta.get("distribution", {})
    is_published = registry.published(name, version) if published is None else published
    download_size = distribution.get("download_size")
    installed_size = distribution.get("installed_size")

    if is_published:
        assets = meta.get("assets", {})
        archive = meta.get("archive")
        if download_size is None:
            if archive and archive.get("parts"):
                part_sizes = [part.get("size") for part in archive["parts"]]
                if part_sizes and all(type(value) is int for value in part_sizes):
                    download_size = sum(part_sizes)
            elif assets and not any(spec.get("decompress") for spec in assets.values()):
                asset_sizes = [spec.get("size") for spec in assets.values()]
                if asset_sizes and all(type(value) is int for value in asset_sizes):
                    download_size = sum(asset_sizes)
        if installed_size is None and assets:
            asset_sizes = [spec.get("size") for spec in assets.values()]
            if asset_sizes and all(type(value) is int for value in asset_sizes):
                installed_size = sum(asset_sizes)
    else:
        if download_size is None:
            manifest = registry.info(name, version)
            source_sizes = [source.size for source in manifest.sources]
            if source_sizes and all(type(value) is int for value in source_sizes):
                download_size = sum(source_sizes)

    return download_size, installed_size


@dataclass(frozen=True)
class DataStatus:
    """Identity, availability, sizes and acquisition route for one ORFdata resource."""

    kind: str
    name: str
    version: Optional[str]
    title: str
    species: tuple[str, ...]
    taxon_ids: tuple[int, ...]
    contexts: tuple[str, ...]
    status: str
    path: Optional[str]
    origin: Optional[str]
    acquire: Optional[str]
    command: Optional[str]
    restricted: bool = False
    default: bool = True
    download_size: Optional[int] = None
    installed_size: Optional[int] = None
    datacache_status: Optional[str] = None

    @property
    def installed(self) -> bool:
        """Whether the resource is usable from local files."""
        return self.status in ("available", "bundled")

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-ready representation, as other OpenVax status objects do."""
        return dict(asdict(self), installed=self.installed)


def default_collections(species: str | int) -> tuple[str, ...]:
    """The reference collections ``install_species`` installs for one species.

    A collection is a versioned, downloadable reference set that is not an ORF
    catalogue. Which ones a species has is read from the registry, so adding a
    collection for another species adds it here too.
    """
    from .viral_registry import viral_collections

    return tuple(collection.name for collection in viral_collections(species=species))


# The command group that manages each kind a registry defines.
_GROUPS = {"catalog": "catalogs", "dataset": "datasets", "reference": "datasets"}


def install_command(name, kind, version=None, *, cache_root=None, restricted=False, force=False, verb=None,
                    registry_file=None):
    """The command that installs one resource in ``cache_root``, as a shell string.

    By default it is ``orfdata data fetch``, which downloads or builds as the
    release declares. ``verb`` names a command of the kind's group instead
    (``download``, ``build`` or ``download-reference``), as an explicit build
    or a ``registry_file`` of custom definitions needs.
    """
    if verb is None:
        args = ["orfdata", "data", "fetch", name, "--kind", kind]
        if cache_root is not None:
            args += ["--cache-root", str(cache_root)]
    else:  # The group's options come before its subcommand.
        args = ["orfdata", _GROUPS[kind]] + ([] if registry_file is None else ["--registry", str(registry_file)])
        args += ([] if cache_root is None else ["--cache-root", str(cache_root)]) + [verb]
        dataset, _, asset = name.partition(":")
        args += [dataset, "--asset", asset] if kind == "reference" else [name]
    if version is not None:
        args += ["--version", version]
    if restricted:
        args.append("--accept-terms")
    if force:
        args.append("--force")
    return join(args)


def install_hint(registry, call, name, kind, version, *, force=False, build=False):
    """``"<Python call>, or on the command line: <command>"`` for one resource of ``registry``, in its cache.

    ``call`` is a method call on the registry, such as ``download('x', '1')``.
    ``build=True`` asks for an explicit source build of a published release. A
    registry made from custom definitions is not what ``orfdata data fetch``
    reads, so its command is the group's own, given the same ``--registry``
    file. The command comes last and unquoted, so it can be copied as it is.
    """
    root = getattr(registry, "_cache_root", None)
    if getattr(registry, "_packaged", True):
        owner = f"{type(registry).__name__}({'' if root is None else f'cache_root={str(root)!r}'})"
        command = install_command(name, kind, version, cache_root=root, force=force,
                                  verb="build" if build else None)
        return f"{owner}.{call}, or on the command line: {command}"
    command = install_command(name, kind, version, cache_root=root, force=force,
                              verb=call.split("(", 1)[0].replace("_", "-"))
    return (f"registry.{call}, or on the command line, adding the same --registry file after "
            f"'orfdata {_GROUPS[kind]}': {command}")


def _fetch_command(row, cache_root=None):
    if row.acquire is None:
        return None
    return install_command(row.name, row.kind, row.version, cache_root=cache_root, restricted=row.restricted)


def needs_replacement(status: str, *, force: bool = False) -> bool:
    """Whether acquiring over a resource's local state replaces its files; ``data fetch`` and setup share this.

    Files that fail verification are replaced; valid ones only with ``force``.
    Check :func:`refuse_unreadable` first: unreadable files aren't replaced unless forced.
    """
    return force or normalize_status(status) not in ("available", "bundled", "missing")


def refuse_unreadable(status: str, *, force: bool = False, name: str = "", command: str | None = None) -> None:
    """Raise for unreadable files unless ``force``: replacing them could hide a permissions problem, or fail on
    it halfway through a download. ``command`` is the command line that would replace them with ``--force``."""
    if normalize_status(status) == "inaccessible" and not force:
        retry = f"force=True, or on the command line: {command} --force" if command else "force=True"
        raise ValueError(f"{name or 'The resource'} can't be read; fix its permissions, choose another "
                         f"cache root, or replace it with {retry}")


def collection_state(
    name: str,
    *,
    species: str | int | None = None,
    cache_root: str | PathLike[str] | None = None,
    verify_files: bool = True,
) -> str:
    """Local state of one reference collection: ``available``, ``missing``, ``damaged`` or ``inaccessible``.

    ``verify_files=False`` checks only that the expected files are present,
    without hashing them.
    """
    from .viral_panel import inspect_viral_collection

    return inspect_viral_collection(
        name, species=species, cache_root=cache_root, verify_files=verify_files,
    ).status


def install_collection(
    name: str,
    *,
    species: str | int | None = None,
    force: bool = False,
    cache_root: str | PathLike[str] | None = None,
    timeout: float | None = DOWNLOAD_TIMEOUT,
    source_paths: dict[str, str | PathLike[str]] | None = None,
    download_options: dict[str, object] | None = None,
):
    """Download or build one reference collection, as ``orfdata data fetch`` does."""
    from .viral_panel import download_viral_references
    from .viral_registry import viral_collection

    collection = viral_collection(name, species=species)
    return download_viral_references(
        collection=collection,
        cache_root=cache_root,
        force=force,
        timeout=timeout,
        source_paths=source_paths,
        download_options=download_options,
    )


def _viral(cache_root, verify_files):
    """Every published viral collection, whichever species it targets."""
    from .viral_panel import bundle_version, inspect_viral_collection, viral_store_path
    from .viral_registry import viral_collections

    rows = []
    for collection in viral_collections():
        species = canonical_species(collection.target_species)
        inspected = inspect_viral_collection(
            collection, species=species, cache_root=cache_root, verify_files=verify_files,
        )
        acquire = "download" if collection.published else "build"
        area = "downloads" if collection.published else "collections"
        download_size = collection.download_size
        installed_size = collection.installed_size
        if collection.published and download_size is None:
            sizes = [spec.get("size") for spec in collection.assets.values()]
            if sizes and all(type(value) is int for value in sizes):
                download_size = sum(sizes)
        if collection.published and installed_size is None:
            sizes = [spec.get("size") for spec in collection.assets.values()]
            if sizes and all(type(value) is int for value in sizes):
                installed_size = sum(sizes)
        rows.append(DataStatus(
            kind="viral",
            name=collection.name,
            version=bundle_version(collection),
            title=f"{collection.title}: genome, transcriptome and proteome layers",
            species=(species,),
            taxon_ids=(collection.target_taxon_id,),
            contexts=collection.contexts,
            status=inspected.status,
            path=str(inspected.path or viral_store_path(cache_root, species, collection, area=area)),
            origin=inspected.origin if inspected.status == "available" else None,
            acquire=acquire,
            command=None,
            default=True,
            download_size=download_size,
            installed_size=installed_size,
            datacache_status=inspected.datacache_status,
        ))
    return tuple(rows)


def _catalogs(cache_root, verify_files):
    from importlib.resources import files

    from .catalog import _BUNDLED, catalog_status
    from .catalogs import CatalogRegistry

    registry = CatalogRegistry(cache_root=cache_root)
    rows = []
    for item in catalog_status(registry=registry):
        inspection = None
        if item.acquire is None:
            bundled_resource = files("orfdata").joinpath(
                f"data/{_BUNDLED[item.species, item.genome_build]}/catalog.json.gz"
            )
            status = "bundled"
            datacache_status = None
            path = str(bundled_resource)
            origin = "bundled"
            is_default = True
            download_size = None
            try:
                installed_size = bundled_resource.stat().st_size
            except (AttributeError, OSError):
                installed_size = len(bundled_resource.read_bytes())
        else:
            inspection = registry.inspect(item.name, item.version, verify_files=verify_files)
            datacache_status = inspection.status
            status = normalize_status(datacache_status)
            path = str(inspection.path)
            origin = item.installed if status == "available" else None
            is_default = item.version == registry.resolve_version(item.name)
            download_size, installed_size = release_sizes(registry, item.name, item.version)
        rows.append(DataStatus(
            kind="catalog",
            name=item.name,
            version=item.version,
            title=item.title,
            species=(item.species,),
            taxon_ids=(species_taxon(item.species),),
            contexts=(),
            status=status,
            path=path,
            origin=origin,
            acquire=item.acquire,
            command=None,
            default=is_default,
            download_size=download_size,
            installed_size=installed_size,
            datacache_status=datacache_status,
        ))
    return tuple(rows)


def _datasets(cache_root, verify_files):
    from .datasets import DatasetRegistry

    registry = DatasetRegistry(cache_root=cache_root)
    rows = []
    for manifest in registry.list():
        name, version = manifest.dataset_id, manifest.version
        published = registry.published(name, version)
        inspection = registry.inspect(name, version, verify_files=verify_files)
        status = normalize_status(inspection.status)
        species = tuple(sorted({canonical_species(reference.species)
                                for reference in manifest.references}))
        taxa = tuple(sorted({reference.taxon_id for reference in manifest.references}))
        acquire = "download" if published else "build"
        download_size, installed_size = release_sizes(
            registry, name, version, published=published,
        )
        rows.append(DataStatus(
            kind="dataset",
            name=name,
            version=version,
            title=manifest.title,
            species=species,
            taxon_ids=taxa,
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin=("download" if published else "build") if status == "available" else None,
            acquire=acquire,
            command=None,
            default=version == registry.resolve_version(name),
            download_size=download_size,
            installed_size=installed_size,
            datacache_status=inspection.status,
        ))
    return tuple(rows)


def _references(cache_root, verify_files):
    from .datasets import DatasetRegistry

    registry = DatasetRegistry(cache_root=cache_root)
    rows = []
    for manifest in registry.list():
        is_default = manifest.version == registry.resolve_version(manifest.dataset_id)
        for spec in registry.search_references(manifest.dataset_id, manifest.version):
            inspection = registry.inspect_reference(
                manifest.dataset_id, manifest.version, reference=spec.reference.key,
                asset_id=spec.source.asset_id, verify_files=verify_files)
            status = normalize_status(inspection.status)
            size = spec.source.size if type(spec.source.size) is int else None
            rows.append(DataStatus(
                kind="reference",
                name=f"{manifest.dataset_id}:{spec.source.asset_id}",
                version=manifest.version,
                title=f"{manifest.title}: native {spec.namespace} protein search reference",
                species=(canonical_species(spec.reference.species),),
                taxon_ids=(spec.reference.taxon_id,),
                contexts=(),
                status=status,
                path=str(inspection.path),
                origin="download" if status == "available" else None,
                acquire="download",
                command=None,
                default=is_default,
                download_size=size,
                installed_size=size,
                datacache_status=inspection.status,
            ))
    return tuple(rows)


def _genomes(cache_root, verify_files):
    from .genomes import genome_registry, pinned_genomes

    rows = []
    for genome in pinned_genomes():
        inspection = genome_registry(genome, cache_root).inspect(
            genome.name, "1", verify_files=verify_files)
        status = normalize_status(inspection.status)
        asset_sizes = [asset.get("size") for asset in genome.assets.values()]
        download_size = (sum(asset_sizes)
                         if asset_sizes and all(type(value) is int for value in asset_sizes)
                         else None)
        rows.append(DataStatus(
            kind="genome",
            name=genome.name,
            version=genome.annotation_version,
            title=genome.title,
            species=(genome.species,),
            taxon_ids=(species_taxon(genome.species),),
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin="download" if status == "available" else None,
            acquire="download",
            command=None,
            default=True,
            download_size=download_size,
            installed_size=None,
            datacache_status=inspection.status,
        ))
    return tuple(rows)


def _chains(cache_root, verify_files):
    from .liftover import chain_registry, pinned_chains

    rows = []
    for chain in pinned_chains():
        inspection = chain_registry(chain, cache_root).inspect(
            chain.name, "1", verify_files=verify_files)
        status = normalize_status(inspection.status)
        size = chain.size if type(chain.size) is int else None
        rows.append(DataStatus(
            kind="chain",
            name=chain.name,
            version="1",
            title=f"{chain.source_build} to {chain.target_build} chain ({chain.provider})",
            species=(chain.species,),
            taxon_ids=(species_taxon(chain.species),),
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin="download" if status == "available" else None,
            acquire="download",
            command=None,
            restricted=chain.restricted,
            default=True,
            download_size=size,
            installed_size=size,
            datacache_status=inspection.status,
        ))
    return tuple(rows)


def _records(cache_root, verify_files):
    from .source_records import pinned_records, records_registry

    rows = []
    for pin in pinned_records():
        inspection = records_registry(pin, cache_root).inspect(
            pin.name, "1", verify_files=verify_files)
        status = normalize_status(inspection.status)
        size = pin.size if type(pin.size) is int else None
        rows.append(DataStatus(
            kind="records",
            name=pin.name,
            version="1",
            title=pin.title,
            species=(pin.species,),
            taxon_ids=(species_taxon(pin.species),),
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin="download" if status == "available" else None,
            acquire="download",
            command=None,
            default=True,
            download_size=size,
            installed_size=size,
            datacache_status=inspection.status,
        ))
    return tuple(rows)


def list_data(
    *,
    kind: str | None = None,
    species: str | int | None = None,
    contexts: Iterable[str] | None = None,
    context_match: str = "any",
    installed: bool | None = None,
    local: bool = False,
    all_versions: bool = False,
    cache_root: str | PathLike[str] | None = None,
    verify_files: bool = False,
) -> tuple[DataStatus, ...]:
    """List known data resources and local state without acquiring anything.

    ``kind`` may be ``viral``, ``catalog``, ``dataset``, ``reference``,
    ``genome``, ``chain`` or ``records``. ``species`` accepts ORFdata's species aliases.
    ``installed`` can be true or false to select usable resources or
    unavailable ones; ``local=True`` selects every resource with local files,
    usable, damaged or unreadable, as ``orfdata data status`` shows them.
    By default only each resource's default version is listed; pass
    ``all_versions=True`` to include superseded and historical versions.
    ``verify_files=True`` hashes pinned files and validates the viral bundle.
    """
    kinds = ("viral", "catalog", "dataset", "reference", "genome", "chain", "records")
    if kind == "all":
        kind = None
    if kind is not None and kind not in kinds:
        raise ValueError(f"Unknown data kind {kind!r}; choose from {', '.join(kinds)}")
    if installed is not None and type(installed) is not bool:
        raise ValueError("installed must be True, False or None")
    if type(local) is not bool:
        raise ValueError("local must be True or False")
    if type(all_versions) is not bool:
        raise ValueError("all_versions must be True or False")
    if contexts is None:
        selected_contexts = ()
    else:
        from .viral_coverage import normalize_viral_contexts

        # A misspelled context is an error here too, not an empty listing.
        selected_contexts = normalize_viral_contexts(contexts)
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    selected_species = None
    if species is not None:
        if species_taxon(species) is None:
            raise ValueError(f"Unknown species alias {species!r}")
        selected_species = canonical_species(species)

    finders = {
        "viral": _viral,
        "catalog": _catalogs,
        "dataset": _datasets,
        "reference": _references,
        "genome": _genomes,
        "chain": _chains,
        "records": _records,
    }
    rows = []
    for current in kinds if kind is None else (kind,):
        rows.extend(finders[current](cache_root, verify_files))
    if not all_versions and not local:
        rows = [row for row in rows if row.default]
    if selected_species is not None:
        rows = [row for row in rows if selected_species in row.species]
    if selected_contexts:
        predicate = any if context_match == "any" else all
        rows = [row for row in rows
                if predicate(context in row.contexts for context in selected_contexts)]
    if installed is not None:
        rows = [row for row in rows if row.installed is installed]
    if local:  # Bundled, installed, damaged or unreadable: anything not simply absent.
        rows = [row for row in rows if row.status != "missing"]
    return tuple(replace(row, command=_fetch_command(row, cache_root)) for row in
                 sorted(rows, key=lambda row: (row.kind, row.name, row.version or "")))
