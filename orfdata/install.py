"""Install everything orfdata uses for one species' extra ORFs, in one explicit step."""

from shlex import join

from .net import DOWNLOAD_TIMEOUT


class InstallError(ValueError):
    """Some of a species' catalogues or chain files failed to install; ``result`` says what did and didn't."""

    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def install_species(species, *, accept_terms=False, force=False, cache_root=None, registry=None, timeout=DOWNLOAD_TIMEOUT,
                    catalogs=None, collections=None):
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
    once. ``catalogs`` chooses which of the species' catalogues, as in
    ``load_catalog``: its default set unless ``"all"`` or names are given.

    Human setup also installs the default GDC viral reference collection.

    The returned mapping adds ``"viral"`` for a human setup. Every item is
    tried; if any fails, :class:`InstallError`
    (a ``ValueError``) names each failure after the rest are installed, and
    carries the same result.
    """
    import tempfile

    from .catalog import _catalog_choices, _installed, _plan, _selected, _species
    from .catalogs import CatalogRegistry
    from .downloads import default_collections, needs_replacement, refuse_unreadable
    from .liftover import _chain_registry, download_chain, pinned_chains

    species = _species(species)
    collections = _selected_collections(default_collections(species), collections, species)
    if registry is not None and cache_root is not None:
        raise ValueError("Pass a registry or a cache_root, not both; the registry has its own cache root")
    registry = CatalogRegistry(cache_root=cache_root) if registry is None else registry
    root = getattr(registry, "_cache_root", cache_root)
    found = _selected([c for c in _catalog_choices(registry) if c.species == species], catalogs, species)
    chains = pinned_chains(species=species)
    if not any(c.registry is not None for c in found) and not chains:
        ready = sorted({c.species for c in _catalog_choices(registry) if c.registry is not None}
                       | {c.species for c in pinned_chains()})
        raise ValueError(f"orfdata has no ORF catalogue or chain file to install for {species}; "
                         f"it has them for {', '.join(ready)}")
    # The releases load_catalog uses on each assembly a catalogue is on or the chains reach.
    builds = {c.genome_build for c in found} | {b for c in chains for b in (c.source_build, c.target_build)}
    needed = {}
    for build in sorted(builds):
        for choice, _ in _plan(found, species, build):
            if choice.registry is not None:
                needed.setdefault((choice.name, choice.version), choice)
    result, failed = dict(species=species, catalogs={}, chains={}), []
    retry = join(["orfdata", "setup", "--species", species] + ([] if root is None else
                                                                ["--cache-root", str(root)]))
    options = {} if timeout is None else dict(timeout=timeout)
    states = {key: _installed(choice) for key, choice in needed.items()}
    for key in sorted(needed):
        if states[key] in ("downloaded", "built") and not force:
            result["catalogs"][f"{key[0]} {key[1]}"] = "already installed"
    builds = [key for key in sorted(needed) if f"{key[0]} {key[1]}" not in result["catalogs"]
              and not registry.published(*key)]
    with tempfile.TemporaryDirectory(prefix="orfdata-setup-") as staging:
        shared = _shared_sources([registry.info(*key) for key in builds], staging, options)
        for (name, version) in sorted(needed):
            label = f"{name} {version}"
            if label in result["catalogs"]:
                continue
            state = states[name, version]
            try:
                # The shared policy needs the damaged store's own status: replaced if invalid, refused if unreadable.
                status = _damage(registry, name, version) if state == "damaged" else "missing"
                refuse_unreadable(status, force=force, name=label, command=retry)
                repair = needs_replacement(status, force=force)
                if (name, version) in builds:
                    paths = {s.asset_id: shared[s.sha256] for s in registry.info(name, version).sources
                             if s.sha256 in shared}
                    registry.build(name, version, force=repair, source_paths=paths or None, download_options=options)
                    verb = "built"
                else:
                    registry.download(name, version, force=repair, **options)
                    verb = "downloaded"
                result["catalogs"][label] = "repaired" if state == "damaged" else verb
            except Exception as error:  # Reported with the rest, after everything else is tried.
                result["catalogs"][label] = f"failed: {error}"
                failed.append((label, str(error)))
    for chain in chains:
        if chain.restricted and accept_terms is not True:
            result["chains"][chain.name] = (f"skipped: {chain.provider}'s terms need accept_terms=True "
                                            f"(--accept-terms): {chain.terms}")
            continue
        try:
            state = _chain_registry(chain, root).inspect(chain.name, "1", verify_files=False).status
            if state == "available" and not force:
                result["chains"][chain.name] = "already downloaded"
                continue
            refuse_unreadable(state, force=force, name=chain.name, command=retry)
            download_chain(chain.source_build, chain.target_build, species=chain.species, accept_terms=accept_terms,
                           force=needs_replacement(state, force=force), cache_root=root,
                           download_options=options)
            result["chains"][chain.name] = {"missing": "downloaded", "available": "reinstalled"}.get(state, "repaired")
        except Exception as error:
            result["chains"][chain.name] = f"failed: {error}"
            failed.append((chain.name, str(error)))
    if collections:
        result["collections"] = _install_collections(
            collections, species, force=force, cache_root=root, timeout=timeout, failed=failed, retry=retry)
    if failed:
        raise InstallError(
            f"Could not install {', '.join(name for name, _ in failed)} for {species}: "
            + "; ".join(f"{name}: {message}" for name, message in failed), result)
    return result


def _install_collections(names, species, *, force, cache_root, timeout, failed, retry=None):
    """Install each reference collection, reporting what every one needed."""
    from .downloads import (
        collection_state,
        install_collection,
        needs_replacement,
        refuse_unreadable,
    )

    installed = {}
    for name in names:
        try:
            state = collection_state(name, species=species, cache_root=cache_root)
            if state == "available" and not force:
                installed[name] = "already installed"
                continue
            refuse_unreadable(state, force=force, name=name, command=retry)
            install_collection(name, species=species, force=needs_replacement(state, force=force),
                               cache_root=cache_root, timeout=timeout)
            installed[name] = {"missing": "downloaded",
                               "available": "reinstalled"}.get(state, "repaired")
        except Exception as error:  # Reported with the rest, after everything else is tried.
            installed[name] = f"failed: {error}"
            failed.append((name, str(error)))
    return installed


def _damage(registry, name, version):
    """The status of whichever of a release's stores is damaged: its download, else its local build."""
    statuses = [registry.inspect_build(name, version, verify_files=False).status]
    if registry.published(name, version):
        statuses.insert(0, registry.inspect(name, version, verify_files=False).status)
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
            f"it has {', '.join(available) or 'none'}")
    return tuple(dict.fromkeys(names))


def _shared_sources(manifests, staging, options):
    """Download once each source two or more of these builds pin; ``{sha256: path}``.

    A source that fails here is left for each build to download, and to
    report, itself.
    """
    from collections import Counter

    uses = Counter(key for m in manifests for key in {(s.url, s.sha256, s.size) for s in m.sources})
    shared = {}
    for (url, sha256, size), count in sorted(uses.items()):
        if count > 1:
            try:
                shared[sha256] = _fetch(url, sha256, size, staging, options)
            except Exception:
                continue
    return shared


def _fetch(url, sha256, size, staging, options):
    """One pinned source file, downloaded and checked into ``staging``."""
    from pathlib import Path

    from .datasets import _datacache

    name = url.rsplit("/", 1)[-1]
    registry = _datacache().VersionedDatasetRegistry(
        {"source": dict(default_version="1", versions={"1": {name: dict(url=url, sha256=sha256, size=size)}})},
        store_path=lambda dataset, version: Path(staging) / sha256)
    registry.download("source", "1", **options)
    return str(registry.local_path("source", "1", asset=name))
