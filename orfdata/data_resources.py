"""Named access to optional data, independent of source and species."""

from dataclasses import replace
from pathlib import Path

from .downloads import _fetch_command, list_data
from .errors import NotInstalledError


def _registry(kind, cache_root):
    from .catalogs import CatalogRegistry
    from .datasets import DatasetRegistry

    return {"catalog": CatalogRegistry, "dataset": DatasetRegistry}[kind](cache_root=cache_root)


def _resolve(name, kind, version, cache_root, verify_files):
    rows = [r for r in list_data(kind=kind, cache_root=cache_root, verify_files=False)
            if r.name == name]
    if not rows:
        raise ValueError(f"Unknown data resource {name!r}; see 'orfdata data list'")
    if len({r.kind for r in rows}) != 1:
        raise ValueError(f"Ambiguous data resource {name!r}; select --kind")
    row = rows[0]
    if row.kind in ("dataset", "catalog") and row.acquire is not None:
        registry = _registry(row.kind, cache_root)
        manifest = registry.info(name, version)
        inspection = registry.inspect(name, manifest.version, verify_files=verify_files)
        acquire = "download" if registry.published(name, manifest.version) else "build"
        row = replace(row, version=manifest.version, title=manifest.title,
                      status=inspection.status, path=str(inspection.path), acquire=acquire,
                      origin=acquire if inspection.status == "available" else None)
    elif row.kind == "reference":
        from .datasets import DatasetRegistry

        registry = DatasetRegistry(cache_root=cache_root)
        dataset, asset = name.split(":", 1)
        manifest = registry.info(dataset, version)
        inspection = registry.inspect_reference(dataset, manifest.version, asset_id=asset,
                                                verify_files=verify_files)
        row = replace(row, version=manifest.version, status=inspection.status,
                      path=str(inspection.path),
                      origin="download" if inspection.status == "available" else None)
    else:
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        if verify_files:
            row = next(r for r in list_data(kind=row.kind, cache_root=cache_root, verify_files=True)
                       if r.name == name and r.version == row.version)
    return replace(row, command=_fetch_command(row))


def data_info(name, *, kind=None, version=None, cache_root=None, verify_files=False):
    """Describe a named resource, its pins and acquisition route, entirely offline."""
    row = _resolve(name, kind, version, cache_root, verify_files)
    result = row.to_dict()
    if row.kind in ("catalog", "dataset") and row.acquire is not None:
        registry = _registry(row.kind, cache_root)
        result["manifest"] = registry.info(name, row.version).to_dict()
        result.update(registry.download_metadata(name, row.version))
    elif row.kind == "reference":
        from .datasets import DatasetRegistry

        dataset, asset = name.split(":", 1)
        registry = DatasetRegistry(cache_root=cache_root)
        result["sources"] = [s.source.to_dict() for s in registry.search_references(dataset, row.version)
                             if s.source.asset_id == asset]
    elif row.kind == "chain":
        from dataclasses import asdict

        from .liftover import pinned_chains

        result["source"] = asdict(next(c for c in pinned_chains() if c.name == name))
    elif row.kind == "genome":
        from dataclasses import asdict

        from .genomes import pinned_genomes

        result["source"] = asdict(next(g for g in pinned_genomes() if g.name == name))
    return result


def data_path(name, *, kind=None, version=None, cache_root=None, verify_files=False):
    """Return an installed file/directory; missing data never triggers acquisition."""
    info = data_info(name, kind=kind, version=version, cache_root=cache_root,
                     verify_files=verify_files)
    if not info["installed"] or info["path"] is None:
        if info["status"] == "missing":
            raise NotInstalledError(f"{name!r} is not installed; install it with {info['command']}")
        raise ValueError(f"{name!r} is {info['status']}; install it with {info['command']}")
    return Path(info["path"])


def fetch_data(name, *, kind=None, version=None, cache_root=None, force=False,
               source_paths=None, timeout=60, show_progress=True, accept_terms=False):
    """Install one resource and return its status, downloading or building as declared.

    Downloads retain their checksum validation. Builds use the same pinned source
    adapters as the focused commands. Existing valid files are reused; discovery
    and ordinary queries never invoke this function implicitly.
    """
    row = _resolve(name, kind, version, cache_root, True)
    if source_paths is not None and row.kind not in ("catalog", "dataset", "genome"):
        raise ValueError("Source overrides here require a catalogue, dataset or genome; use the focused command for other kinds")
    if row.restricted and not accept_terms:
        raise ValueError(f"{name!r} requires --accept-terms; see 'orfdata data info {name}'")
    if row.acquire is None or (row.installed and not force):
        return dict(row.to_dict(), acquired="already installed")
    if row.kind in ("catalog", "dataset"):
        registry = _registry(row.kind, cache_root)
        options = dict(timeout=timeout, show_progress=show_progress)
        if row.acquire == "download":
            registry.download(name, row.version, force=force, source_paths=source_paths, **options)
        else:
            registry.build(name, row.version, force=force, source_paths=source_paths,
                           download_options=options)
    elif row.kind == "reference":
        from .datasets import DatasetRegistry

        dataset, asset = name.split(":", 1)
        DatasetRegistry(cache_root=cache_root).download_reference(
            dataset, row.version, asset_id=asset, force=force, timeout=timeout,
            show_progress=show_progress)
    elif row.kind == "genome":
        from .genomes import install_genome

        install_genome(name, cache_root=cache_root, force=force, source_paths=source_paths,
                       download_options=dict(timeout=timeout, show_progress=show_progress))
    elif row.kind == "chain":
        from .liftover import download_chain, pinned_chains

        chain = next(c for c in pinned_chains() if c.name == name)
        download_chain(chain.source_build, chain.target_build, species=chain.species,
                       cache_root=cache_root, force=force, accept_terms=accept_terms,
                       download_options=dict(timeout=timeout, show_progress=show_progress))
    elif row.kind == "viral":
        from .downloads import install_collection

        install_collection(name, cache_root=cache_root, force=force, timeout=timeout)
    result = _resolve(name, row.kind, row.version, cache_root, False)
    if not result.installed:
        raise RuntimeError(f"{name!r} installation finished with status {result.status}")
    return dict(result.to_dict(), acquired="downloaded" if row.acquire == "download" else "built")
