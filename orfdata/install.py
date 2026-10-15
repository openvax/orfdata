"""Install everything orfdata uses for one species' extra ORFs, in one explicit step."""


class InstallError(ValueError):
    """Some of a species' catalogues or chain files failed to install; ``result`` says what did and didn't."""

    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def install_species(species, *, accept_terms=False, force=False, cache_root=None, registry=None):
    """Download or build each registered catalogue for ``species`` and download its pinned chain files.

    Afterwards ``load_catalog``, ``extra_orfs`` and ``Proteome.build`` work
    offline on every assembly of the species the chains reach: every release
    one of those assemblies uses is installed. One already installed is kept
    unless ``force``, and a damaged one is repaired. Chains with usage
    restrictions, such as UCSC's non-commercial mouse chains, need
    ``accept_terms=True``.

    Returns ``{"species", "catalogs", "chains"}``, each item mapped to what
    happened to it. Every item is tried; if any fails, :class:`InstallError`
    (a ``ValueError``) names each failure after the rest are installed, and
    carries the same result.
    """
    from .catalog import _catalog_choices, _installed, _plan, _species
    from .catalogs import CatalogRegistry
    from .liftover import _chain_registry, download_chain, pinned_chains

    species = _species(species)
    if registry is not None and cache_root is not None:
        raise ValueError("Pass a registry or a cache_root, not both; the registry has its own cache root")
    registry = CatalogRegistry(cache_root=cache_root) if registry is None else registry
    root = getattr(registry, "_cache_root", cache_root)
    found = [c for c in _catalog_choices(registry) if c.species == species]
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
    for (name, version), choice in sorted(needed.items()):
        label = f"{name} {version}"
        state = _installed(choice)
        try:
            if state in ("downloaded", "built") and not force:
                result["catalogs"][label] = "already installed"
                continue
            repair = force or state == "damaged"
            verb = "downloaded" if registry.published(name, version) else "built"
            install = registry.download if verb == "downloaded" else registry.build
            install(name, version, force=repair)
            result["catalogs"][label] = "repaired" if state == "damaged" else verb
        except Exception as error:  # Reported with the rest, after everything else is tried.
            result["catalogs"][label] = f"failed: {error}"
            failed.append(label)
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
            download_chain(chain.source_build, chain.target_build, species=chain.species, accept_terms=accept_terms,
                           force=force or state != "missing", cache_root=root)
            result["chains"][chain.name] = "downloaded" if state == "missing" else "repaired"
        except Exception as error:
            result["chains"][chain.name] = f"failed: {error}"
            failed.append(chain.name)
    if failed:
        raise InstallError(f"Could not install {', '.join(failed)} for {species}: "
                           + "; ".join(f"{item}: {(result['catalogs'].get(item) or result['chains'][item])[8:]}"
                                       for item in failed), result)
    return result
