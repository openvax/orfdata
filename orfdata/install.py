"""Install everything orfdata uses for one species' extra ORFs, in one explicit step."""


def install_species(species, *, accept_terms=False, force=False, cache_root=None, registry=None):
    """Download or build each registered catalogue for ``species`` and download its pinned chain files.

    Afterwards ``load_catalog``, ``extra_orfs`` and ``Proteome.build`` work
    offline on every assembly of the species the chains reach. Each
    catalogue's default release is installed; one already installed is
    kept unless ``force``. Chains with usage restrictions, such as UCSC's
    non-commercial mouse chains, need ``accept_terms=True``. Returns what was
    installed, kept or skipped.
    """
    from .catalog import _Choice, _installed, _species
    from .catalogs import CatalogRegistry
    from .liftover import download_chains, pinned_chains

    species = _species(species)
    if registry is not None and cache_root is not None:
        raise ValueError("Pass a registry or a cache_root, not both; the registry has its own cache root")
    registry = CatalogRegistry(cache_root=cache_root) if registry is None else registry
    result = dict(species=species, catalogs={}, chains=None)
    names = dict.fromkeys(m.catalog_id for m in registry.list(species=species))
    if not names and not any(c.species == species for c in pinned_chains()):
        raise ValueError(f"orfdata has no ORF catalogue or chain file to install for {species}")
    for name in names:
        manifest = registry.info(name)
        choice = _Choice(species, manifest.genome_build, name, registry, version=manifest.version)
        label = f"{name} {manifest.version}"
        state = _installed(choice)
        repair = force or state == "damaged"
        if state in ("downloaded", "built") and not force:
            result["catalogs"][label] = "already installed"
        elif registry.published(name):
            registry.download(name, force=repair)
            result["catalogs"][label] = "repaired" if state == "damaged" else "downloaded"
        else:
            registry.build(name, force=repair)
            result["catalogs"][label] = "repaired" if state == "damaged" else "built"
    if any(c.species == species for c in pinned_chains()):
        result["chains"] = download_chains(species=species, accept_terms=accept_terms, force=force,
                                           cache_root=getattr(registry, "_cache_root", cache_root))
    return result
