"""Install everything orfdata uses for one species' extra ORFs, in one explicit step."""

from dataclasses import dataclass
from shlex import join

from .cache import DOWNLOAD_TIMEOUT


class InstallError(ValueError):
    """Some of a species' catalogues or chain files failed to install; ``result`` says what did and didn't."""

    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


@dataclass(frozen=True)
class _SetupItem:
    """One planned resource acquisition item in :func:`install_species`."""

    bucket: str
    label: str
    kind: str
    name: str
    version: str | None
    restricted: bool = False
    skip_reason: str | None = None
    options: dict | None = None
    state_override: str | None = None


def install_species(species, *, accept_terms=False, force=False, cache_root=None, registry=None, timeout=DOWNLOAD_TIMEOUT,
                    catalogs=None, collections=None, resume=False):
    """Install a species' registered catalogues, chains and default reference collections.

    Afterwards ``load_catalog``, ``extra_orfs`` and ``Proteome.build`` work
    offline on every assembly of the species the chains reach: every release
    one of those assemblies uses is installed. One already installed is kept
    unless ``force``, a damaged one is repaired, and an unreadable one fails
    (and is reported) unless ``force``, as ``orfdata data fetch`` treats them. Chains with usage
    restrictions, such as UCSC's non-commercial mouse chains, need
    ``accept_terms=True``. A download that stalls for ``timeout`` seconds is
    retried, then fails, rather than waiting forever (``None`` waits). A
    source file several catalogues pin, such as an Ensembl GTF, is downloaded
    once into ORFdata's content-addressed SHA-256 source store and reused across
    builds and setup runs. ``catalogs`` chooses which of the species' catalogues,
    as in ``load_catalog``: its default set unless ``"all"`` or names are given.

    ``resume=True`` lets an interrupted catalogue, chain or build-source
    download continue where it stopped, as in :func:`~orfdata.fetch_data`.

    Human setup also installs the default GDC viral reference collection.

    The returned mapping adds ``"collections"`` when reference collections are
    selected. Every item is tried; if any fails, :class:`InstallError`
    (a ``ValueError``) names each failure after the rest are installed, and
    carries the same result.
    """
    from .catalog import catalog_choices, catalog_install_state, plan_catalogs, select_catalogs
    from .catalogs import CatalogRegistry
    from .downloads import (
        acquisition_action,
        default_collections,
        needs_replacement,
        refuse_unreadable,
    )
    from .liftover import pinned_chains
    from .reference import canonical_species
    from .resource_kinds import get_resource_kind

    species = canonical_species(species)
    collections = _selected_collections(default_collections(species), collections, species)
    if registry is not None and cache_root is not None:
        raise ValueError("Pass a registry or a cache_root, not both; the registry has its own cache root")
    registry = CatalogRegistry(cache_root=cache_root) if registry is None else registry
    root = getattr(registry, "_cache_root", cache_root)
    found = select_catalogs([c for c in catalog_choices(registry) if c.species == species], catalogs, species)
    chains = pinned_chains(species=species)
    if not any(c.registry is not None for c in found) and not chains:
        ready = sorted({c.species for c in catalog_choices(registry) if c.registry is not None}
                       | {c.species for c in pinned_chains()})
        raise ValueError(f"orfdata has no ORF catalogue or chain file to install for {species}; "
                         f"it has them for {', '.join(ready)}")
    if type(resume) is not bool:
        raise ValueError("resume must be True or False")

    # Build one unified plan across catalogues, chains and reference collections.
    builds = {c.genome_build for c in found} | {b for c in chains for b in (c.source_build, c.target_build)}
    needed_catalogs = {}
    for build in sorted(builds):
        for choice, _ in plan_catalogs(found, species, build):
            if choice.registry is not None:
                needed_catalogs.setdefault((choice.name, choice.version), choice)

    plan: list[_SetupItem] = []
    for (name, version), choice in sorted(needed_catalogs.items()):
        state = catalog_install_state(choice)[0]
        if state in ("damaged", "inaccessible"):
            state = _damage(registry, name, version)
        plan.append(
            _SetupItem(
                bucket="catalogs",
                label=f"{name} {version}",
                kind="catalog",
                name=name,
                version=version,
                options={"registry": registry},
                state_override=state,
            )
        )
    for chain in chains:
        skip = (
            f"skipped: {chain.provider}'s terms need accept_terms=True (--accept-terms): {chain.terms}"
            if chain.restricted and accept_terms is not True
            else None
        )
        plan.append(
            _SetupItem(
                bucket="chains",
                label=chain.name,
                kind="chain",
                name=chain.name,
                version="1",
                restricted=chain.restricted,
                skip_reason=skip,
            )
        )
    for collection_name in collections:
        plan.append(
            _SetupItem(
                bucket="collections",
                label=collection_name,
                kind="viral",
                name=collection_name,
                version=None,
                options={"species": species},
            )
        )

    result = dict(species=species, catalogs={}, chains={})
    if collections:
        result["collections"] = {}
    failed = []
    retry = join(["orfdata", "setup", "--species", species] + ([] if root is None else ["--cache-root", str(root)]))
    download_options = ({} if timeout is None else dict(timeout=timeout)) | ({"resume": "auto"} if resume else {})

    # Generic acquisition loop over the resource-kind registry: no per-kind branches.
    for item in plan:
        if item.skip_reason is not None:
            result[item.bucket][item.label] = item.skip_reason
            continue
        handler = get_resource_kind(item.kind)
        try:
            row = handler.inspect(
                item.name,
                item.version,
                cache_root=root,
                verify_files=False,
                options=item.options,
            )
            state = item.state_override if item.state_override is not None else row.status
            if state == "available" and not force:
                result[item.bucket][item.label] = "already installed"
                continue
            refuse_unreadable(state, force=force, name=item.label, command=retry)
            repair = needs_replacement(state, force=force)
            handler.acquire(
                item.name,
                item.version,
                cache_root=root,
                force=repair,
                accept_terms=accept_terms,
                download_options=download_options,
                options=item.options,
            )
            result[item.bucket][item.label] = acquisition_action(state, acquire=row.acquire, force=force)
        except Exception as error:  # Reported with the rest, after everything else is tried.
            result[item.bucket][item.label] = f"failed: {error}"
            failed.append((item.label, str(error)))

    if failed:
        raise InstallError(
            f"Could not install {', '.join(name for name, _ in failed)} for {species}: "
            + "; ".join(f"{name}: {message}" for name, message in failed),
            result,
        )
    return result


def _damage(registry, name, version):
    """The unified status of whichever of a release's stores is damaged or inaccessible."""
    from .downloads import normalize_status

    statuses = [normalize_status(registry.inspect_build(name, version, verify_files=False).status)]
    if registry.published(name, version):
        statuses.insert(0, normalize_status(registry.inspect(name, version, verify_files=False).status))
    return next((status for status in statuses if status not in ("available", "missing")), "missing")


def _selected_collections(available, chosen, species):
    """``collections``, as ``catalogs`` works: the default set, "all", "none", or names."""
    if chosen is None:
        return tuple(available)
    names = [chosen] if isinstance(chosen, str) else list(chosen)
    if "all" in names:
        return tuple(available)
    if not names or names == ["none"]:
        return ()
    unknown = [name for name in names if name not in available]
    if unknown:
        raise ValueError(
            f"orfdata has no reference collection {', '.join(map(repr, unknown))} for {species}; "
            f"it has {', '.join(available) or 'none'}"
        )
    return tuple(dict.fromkeys(names))
