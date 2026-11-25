"""Uniform resource-kind protocol and registry for ORFdata discovery, inspection and acquisition.

Every ORFdata resource kind (viral collections, viral inventories, catalogues,
quantitative datasets, native search references, pinned genomes, liftover
chains, and pinned source-record archives) implements :class:`ResourceKind` and
registers itself in :data:`RESOURCE_KINDS`.

``list_data``, ``data_info``, ``data_path``, ``check_fetch``, ``fetch_data``
and ``install_species`` operate generically through this registry without
per-kind ``if/elif`` dispatch.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import asdict
from importlib.resources import files
from os import PathLike
from typing import Protocol, runtime_checkable

from .cache import DOWNLOAD_TIMEOUT
from .downloads import DataStatus, normalize_status, release_sizes
from .reference import canonical_species, species_taxon


def resolve_source_paths(
    source_paths: str | PathLike[str] | Mapping[str, object] | None,
    assets: Iterable[str],
) -> dict[str, object] | None:
    """Normalize unified ``--source-path`` overrides against a resource's declared assets.

    Accepts either a mapping of ``{asset_name: path}`` or a single bare path
    when the resource declares exactly one asset.

    Parameters
    ----------
    source_paths : str, os.PathLike, mapping, or None
        Caller-supplied local override(s).
    assets : iterable of str
        Declared asset identifiers accepted by the resource.

    Returns
    -------
    dict or None
        Normalized ``{asset_name: path}`` mapping, or ``None`` when no override
        was provided.

    Raises
    ------
    ValueError
        If overrides are not supported for the resource, if a bare path is
        given for a multi-asset resource, or if an unknown asset name is used.
    TypeError
        If ``source_paths`` is neither a path nor a mapping.
    """
    if source_paths is None:
        return None
    declared = tuple(assets)
    if isinstance(source_paths, (str, PathLike)):
        if not declared:
            raise ValueError("Source overrides are not supported for this resource")
        if len(declared) != 1:
            raise ValueError(
                "Source overrides with a bare PATH require a resource with one asset; "
                f"choose from {', '.join(declared)} using ASSET=PATH"
            )
        return {declared[0]: os.fspath(source_paths)}
    if isinstance(source_paths, Mapping):
        resolved = {str(key): value for key, value in source_paths.items()}
        if not resolved:
            return None
        if not declared:
            raise ValueError("Source overrides are not supported for this resource")
        unknown = sorted(set(resolved) - set(declared))
        if unknown:
            raise ValueError(
                f"Source overrides name unknown asset(s) {', '.join(unknown)}; "
                f"available: {', '.join(declared)}"
            )
        return resolved
    raise TypeError("source_paths must be a mapping of ASSET -> PATH or a single PATH")


def _reject_extra_options(
    kind: str,
    name: str,
    options: Mapping[str, object] | None,
    allowed: Iterable[str] = (),
) -> dict[str, object]:
    """Validate that ``options`` contains only keys in ``allowed`` with non-default values."""
    if not options:
        return {}
    allowed_set = set(allowed)
    unknown = sorted(key for key, value in options.items() if key not in allowed_set and value is not None)
    if unknown:
        raise ValueError(
            f"Unsupported option(s) for {kind} resource {name!r}: {', '.join(unknown)}"
        )
    return {key: value for key, value in options.items() if key in allowed_set}


@runtime_checkable
class ResourceKind(Protocol):
    """Protocol implemented once per acquirable ORFdata resource kind."""

    kind: str

    def rows(
        self,
        *,
        cache_root: str | PathLike[str] | None = None,
        verify_files: bool = False,
    ) -> tuple[DataStatus, ...]:
        """Return every known :class:`~orfdata.downloads.DataStatus` row of this kind."""
        ...

    def inspect(
        self,
        name: str,
        version: str | None = None,
        *,
        cache_root: str | PathLike[str] | None = None,
        verify_files: bool = False,
        options: Mapping[str, object] | None = None,
    ) -> DataStatus:
        """Inspect one named resource of this kind offline."""
        ...

    def info(
        self,
        name: str,
        version: str | None = None,
        *,
        cache_root: str | PathLike[str] | None = None,
        verify_files: bool = False,
        options: Mapping[str, object] | None = None,
    ) -> dict[str, object]:
        """Return offline JSON-ready metadata and pins for one named resource."""
        ...

    def validate_options(
        self,
        row: DataStatus,
        *,
        cache_root: str | PathLike[str] | None = None,
        source_paths: str | PathLike[str] | Mapping[str, object] | None = None,
        options: Mapping[str, object] | None = None,
    ) -> tuple[dict[str, object] | None, dict[str, object]]:
        """Validate ``source_paths`` and kind-specific ``options`` before acquisition."""
        ...

    def acquire(
        self,
        name: str,
        version: str | None = None,
        *,
        cache_root: str | PathLike[str] | None = None,
        force: bool = False,
        source_paths: str | PathLike[str] | Mapping[str, object] | None = None,
        accept_terms: bool = False,
        download_options: Mapping[str, object] | None = None,
        options: Mapping[str, object] | None = None,
    ) -> object:
        """Download or build one named resource of this kind."""
        ...


class ViralResourceKind:
    """Published or locally built viral reference collections (e.g. ``gdc-viruses``)."""

    kind = "viral"
    _ALLOWED_OPTIONS = (
        "species",
        "prefer_reference_molecule",
        "annotate_from_alternate",
        "prefer_alternate_annotation",
    )

    def _row(self, collection, cache_root, verify_files) -> DataStatus:
        from .viral_panel import bundle_version, inspect_viral_collection, viral_store_path

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
        return DataStatus(
            kind=self.kind,
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
        )

    def rows(self, *, cache_root=None, verify_files=False) -> tuple[DataStatus, ...]:
        from .viral_registry import viral_collections

        return tuple(self._row(collection, cache_root, verify_files) for collection in viral_collections())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        from .viral_registry import viral_collection

        species = None if not options else options.get("species")
        collection = viral_collection(name, species=species)
        row = self._row(collection, cache_root, verify_files)
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        return row

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        return self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options).to_dict()

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        from .viral_registry import viral_collection

        cleaned = _reject_extra_options(self.kind, row.name, options, self._ALLOWED_OPTIONS)
        collection = viral_collection(row.name, species=cleaned.get("species"))
        prefer_ref = cleaned.get("prefer_reference_molecule", True)
        fill_alt = cleaned.get("annotate_from_alternate", True)
        prefer_alt = cleaned.get("prefer_alternate_annotation", False)
        default_content = bool(prefer_ref and fill_alt and not prefer_alt)
        assets = tuple(collection.assets) if (collection.published and default_content) else ()
        resolved = resolve_source_paths(source_paths, assets)
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        from .viral_panel import download_viral_references
        from .viral_registry import viral_collection

        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, cleaned = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        collection = viral_collection(name, species=cleaned.get("species"))
        opts = dict(download_options or {})
        timeout = opts.get("timeout", DOWNLOAD_TIMEOUT)
        return download_viral_references(
            species=cleaned.get("species"),
            collection=collection,
            cache_root=cache_root,
            force=force,
            timeout=timeout,
            prefer_reference_molecule=cleaned.get("prefer_reference_molecule", True),
            annotate_from_alternate=cleaned.get("annotate_from_alternate", True),
            prefer_alternate_annotation=cleaned.get("prefer_alternate_annotation", False),
            source_paths=resolved_paths,
            download_options=opts,
        )


class InventoryResourceKind:
    """Published viral screening inventory files (e.g. ``gdc-viruses-inventory``)."""

    kind = "inventory"

    @staticmethod
    def _collection_name(collection) -> str:
        return f"{collection.name}-inventory"

    def _find(self, name: str):
        from .viral_registry import viral_collections

        for collection in viral_collections():
            if self._collection_name(collection) == name:
                return collection
        raise ValueError(f"Unknown viral inventory {name!r}; see 'orfdata data list --kind inventory'")

    def _row(self, collection, cache_root, verify_files) -> DataStatus:
        from .viral_coverage import inspect_viral_inventory

        inspected = inspect_viral_inventory(
            collection=collection, cache_root=cache_root, verify_files=verify_files,
        )
        status = normalize_status(inspected.status)
        species = canonical_species(collection.target_species)
        return DataStatus(
            kind=self.kind,
            name=self._collection_name(collection),
            version=collection.panel,
            title=f"{collection.title}: published viral inventory ({collection.panel_filename})",
            species=(species,),
            taxon_ids=(collection.target_taxon_id,),
            contexts=(),
            status=status,
            path=str(inspected.path),
            origin="download" if status == "available" else None,
            acquire="download",
            command=None,
            default=True,
            download_size=collection.panel_size,
            installed_size=collection.panel_size,
            datacache_status=inspected.status,
        )

    def rows(self, *, cache_root=None, verify_files=False) -> tuple[DataStatus, ...]:
        from .viral_registry import viral_collections

        return tuple(self._row(collection, cache_root, verify_files) for collection in viral_collections())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        collection = self._find(name)
        row = self._row(collection, cache_root, verify_files)
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        return row

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        collection = self._find(name)
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        result = row.to_dict()
        result["source"] = {
            "collection": collection.name,
            "panel": collection.panel,
            "filename": collection.panel_filename,
            "url": collection.panel_url,
            "sha256": collection.panel_sha256,
            "size": collection.panel_size,
        }
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        cleaned = _reject_extra_options(self.kind, row.name, options)
        collection = self._find(row.name)
        resolved = resolve_source_paths(source_paths, (collection.panel_filename,))
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        from .viral_coverage import download_viral_inventory

        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, _ = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        collection = self._find(name)
        opts = dict(download_options or {})
        timeout = opts.get("timeout", DOWNLOAD_TIMEOUT)
        source_path = None if resolved_paths is None else resolved_paths[collection.panel_filename]
        return download_viral_inventory(
            collection=collection,
            cache_root=cache_root,
            force=force,
            timeout=timeout,
            source_path=source_path,
            download_options=opts,
        )


class CatalogResourceKind:
    """Bundled and registered ORF catalogues."""

    kind = "catalog"
    _ALLOWED_OPTIONS = ("registry", "build")

    @staticmethod
    def _registry(cache_root, options):
        from .catalogs import CatalogRegistry

        if options and options.get("registry") is not None:
            return options["registry"]
        return CatalogRegistry(cache_root=cache_root)

    def _bundled_row(self, item) -> DataStatus:
        from .catalog import bundled_catalogs

        resource_name = bundled_catalogs()[item.species, item.genome_build]
        bundled_resource = files("orfdata").joinpath(f"data/{resource_name}/catalog.json.gz")
        try:
            installed_size = bundled_resource.stat().st_size
        except (AttributeError, OSError):
            installed_size = len(bundled_resource.read_bytes())
        return DataStatus(
            kind=self.kind,
            name=item.name,
            version=item.version,
            title=item.title,
            species=(item.species,),
            taxon_ids=(species_taxon(item.species),),
            contexts=(),
            status="bundled",
            path=str(bundled_resource),
            origin="bundled",
            acquire=None,
            command=None,
            default=True,
            download_size=None,
            installed_size=installed_size,
            datacache_status=None,
        )

    def rows(self, *, cache_root=None, verify_files=False, options=None) -> tuple[DataStatus, ...]:
        from .catalog import catalog_status

        registry = self._registry(cache_root, options)
        rows = []
        for item in catalog_status(registry=registry):
            if item.acquire is None:
                rows.append(self._bundled_row(item))
            else:
                inspection = registry.inspect(item.name, item.version, verify_files=verify_files)
                datacache_status = inspection.status
                status = normalize_status(datacache_status)
                download_size, installed_size = release_sizes(registry, item.name, item.version)
                rows.append(
                    DataStatus(
                        kind=self.kind,
                        name=item.name,
                        version=item.version,
                        title=item.title,
                        species=(item.species,),
                        taxon_ids=(species_taxon(item.species),),
                        contexts=(),
                        status=status,
                        path=str(inspection.path),
                        origin=item.installed if status == "available" else None,
                        acquire=item.acquire,
                        command=None,
                        default=item.version == registry.resolve_version(item.name),
                        download_size=download_size,
                        installed_size=installed_size,
                        datacache_status=datacache_status,
                    )
                )
        return tuple(rows)

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        from .catalog import catalog_status

        registry = self._registry(cache_root, options)
        explicit_build = bool(options and options.get("build"))
        if name == "bundled":
            item = next(row for row in catalog_status(registry=registry) if row.name == "bundled")
            row = self._bundled_row(item)
            if version is not None and version != row.version:
                raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
            return row
        manifest = registry.info(name, version)
        published = registry.published(name, manifest.version)
        use_build = explicit_build or not published
        inspection = (
            registry.inspect_build(name, manifest.version, verify_files=verify_files)
            if use_build
            else registry.inspect(name, manifest.version, verify_files=verify_files)
        )
        status = normalize_status(inspection.status)
        acquire = "build" if use_build else "download"
        download_size, installed_size = release_sizes(
            registry, name, manifest.version, published=not use_build,
        )
        species = canonical_species(manifest.references[0].species)
        return DataStatus(
            kind=self.kind,
            name=name,
            version=manifest.version,
            title=manifest.title,
            species=(species,),
            taxon_ids=(species_taxon(species),),
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin=acquire if status == "available" else None,
            acquire=acquire,
            command=None,
            default=manifest.version == registry.resolve_version(name),
            download_size=download_size,
            installed_size=installed_size,
            datacache_status=inspection.status,
        )

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        result = row.to_dict()
        if row.acquire is not None:
            registry = self._registry(cache_root, options)
            result["manifest"] = registry.info(name, row.version).to_dict()
            result.update(registry.download_metadata(name, row.version))
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        from .catalogs import CATALOG_ASSET

        cleaned = _reject_extra_options(self.kind, row.name, options, self._ALLOWED_OPTIONS)
        if row.acquire is None:
            resolved = resolve_source_paths(source_paths, ())
            return resolved, cleaned
        registry = self._registry(cache_root, cleaned)
        explicit_build = bool(cleaned.get("build"))
        use_build = explicit_build or row.acquire == "build"
        if use_build:
            manifest = registry.info(row.name, row.version)
            assets = tuple(source.asset_id for source in manifest.sources)
        else:
            assets = (CATALOG_ASSET,)
        resolved = resolve_source_paths(source_paths, assets)
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, cleaned = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        if row.acquire is None:
            return None
        registry = self._registry(cache_root, cleaned)
        opts = dict(download_options or {})
        if row.acquire == "download":
            return registry.download(
                name, row.version, force=force, source_paths=resolved_paths, **opts,
            )
        return registry.build(
            name, row.version, force=force, source_paths=resolved_paths, download_options=opts,
        )


class DatasetResourceKind:
    """Registered quantitative datasets (single-file, archived, or locally built)."""

    kind = "dataset"
    _ALLOWED_OPTIONS = ("registry", "build")

    @staticmethod
    def _registry(cache_root, options):
        from .datasets import DatasetRegistry

        if options and options.get("registry") is not None:
            return options["registry"]
        return DatasetRegistry(cache_root=cache_root)

    def _row(self, registry, manifest, *, verify_files, explicit_build=False) -> DataStatus:
        name, version = manifest.dataset_id, manifest.version
        published = registry.published(name, version)
        use_build = explicit_build or not published
        inspection = (
            registry.inspect_build(name, version, verify_files=verify_files)
            if use_build
            else registry.inspect(name, version, verify_files=verify_files)
        )
        status = normalize_status(inspection.status)
        species = tuple(sorted({canonical_species(ref.species) for ref in manifest.references}))
        taxa = tuple(sorted({ref.taxon_id for ref in manifest.references}))
        acquire = "build" if use_build else "download"
        download_size, installed_size = release_sizes(
            registry, name, version, published=not use_build,
        )
        return DataStatus(
            kind=self.kind,
            name=name,
            version=version,
            title=manifest.title,
            species=species,
            taxon_ids=taxa,
            contexts=(),
            status=status,
            path=str(inspection.path),
            origin=acquire if status == "available" else None,
            acquire=acquire,
            command=None,
            default=version == registry.resolve_version(name),
            download_size=download_size,
            installed_size=installed_size,
            datacache_status=inspection.status,
        )

    def rows(self, *, cache_root=None, verify_files=False, options=None) -> tuple[DataStatus, ...]:
        registry = self._registry(cache_root, options)
        return tuple(self._row(registry, manifest, verify_files=verify_files) for manifest in registry.list())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        registry = self._registry(cache_root, options)
        manifest = registry.info(name, version)
        explicit_build = bool(options and options.get("build"))
        return self._row(registry, manifest, verify_files=verify_files, explicit_build=explicit_build)

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        registry = self._registry(cache_root, options)
        result = row.to_dict()
        result["manifest"] = registry.info(name, row.version).to_dict()
        result.update(registry.download_metadata(name, row.version))
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        from .datasets import DATASET_ASSET

        cleaned = _reject_extra_options(self.kind, row.name, options, self._ALLOWED_OPTIONS)
        registry = self._registry(cache_root, cleaned)
        explicit_build = bool(cleaned.get("build"))
        use_build = explicit_build or row.acquire == "build"
        if use_build:
            manifest = registry.info(row.name, row.version)
            assets = tuple(source.asset_id for source in manifest.sources)
        else:
            assets = (DATASET_ASSET,)
        resolved = resolve_source_paths(source_paths, assets)
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, cleaned = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        registry = self._registry(cache_root, cleaned)
        opts = dict(download_options or {})
        if row.acquire == "download":
            return registry.download(
                name, row.version, force=force, source_paths=resolved_paths, **opts,
            )
        return registry.build(
            name, row.version, force=force, source_paths=resolved_paths, download_options=opts,
        )


class ReferenceResourceKind:
    """Native protein search FASTAs declared by datasets (``DATASET:ASSET``)."""

    kind = "reference"
    _ALLOWED_OPTIONS = ("registry",)

    @staticmethod
    def _registry(cache_root, options):
        from .datasets import DatasetRegistry

        if options and options.get("registry") is not None:
            return options["registry"]
        return DatasetRegistry(cache_root=cache_root)

    @staticmethod
    def _split(name: str) -> tuple[str, str]:
        if ":" not in name:
            raise ValueError(f"Unknown search reference {name!r}; expected DATASET:ASSET")
        dataset, asset = name.split(":", 1)
        return dataset, asset

    def _row(self, registry, manifest, spec, *, verify_files) -> DataStatus:
        inspection = registry.inspect_reference(
            manifest.dataset_id,
            manifest.version,
            reference=spec.reference.key,
            asset_id=spec.source.asset_id,
            verify_files=verify_files,
        )
        status = normalize_status(inspection.status)
        size = spec.source.size if type(spec.source.size) is int else None
        return DataStatus(
            kind=self.kind,
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
            default=manifest.version == registry.resolve_version(manifest.dataset_id),
            download_size=size,
            installed_size=size,
            datacache_status=inspection.status,
        )

    def rows(self, *, cache_root=None, verify_files=False, options=None) -> tuple[DataStatus, ...]:
        registry = self._registry(cache_root, options)
        rows = []
        for manifest in registry.list():
            for spec in registry.search_references(manifest.dataset_id, manifest.version):
                rows.append(self._row(registry, manifest, spec, verify_files=verify_files))
        return tuple(rows)

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        registry = self._registry(cache_root, options)
        dataset, asset = self._split(name)
        manifest = registry.info(dataset, version)
        specs = [
            spec for spec in registry.search_references(dataset, manifest.version)
            if spec.source.asset_id == asset
        ]
        if not specs:
            raise ValueError(f"Unknown search reference {name!r}; see 'orfdata data list --kind reference'")
        return self._row(registry, manifest, specs[0], verify_files=verify_files)

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        registry = self._registry(cache_root, options)
        dataset, asset = self._split(name)
        result = row.to_dict()
        result["sources"] = [
            spec.source.to_dict()
            for spec in registry.search_references(dataset, row.version)
            if spec.source.asset_id == asset
        ]
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        cleaned = _reject_extra_options(self.kind, row.name, options, self._ALLOWED_OPTIONS)
        _, asset = self._split(row.name)
        resolved = resolve_source_paths(source_paths, (asset,))
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, cleaned = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        registry = self._registry(cache_root, cleaned)
        dataset, asset = self._split(name)
        source_path = None if resolved_paths is None else resolved_paths[asset]
        return registry.download_reference(
            dataset,
            row.version,
            asset_id=asset,
            force=force,
            source_path=source_path,
            **(download_options or {}),
        )


class GenomeResourceKind:
    """Pinned genome annotations (e.g. ``t2t-chm13``)."""

    kind = "genome"
    _ALLOWED_OPTIONS = ("genome_fasta_path_or_url",)

    def _row(self, genome, cache_root, verify_files) -> DataStatus:
        from .genomes import genome_registry

        inspection = genome_registry(genome, cache_root).inspect(
            genome.name, "1", verify_files=verify_files,
        )
        status = normalize_status(inspection.status)
        asset_sizes = [asset.get("size") for asset in genome.assets.values()]
        download_size = (
            sum(asset_sizes)
            if asset_sizes and all(type(value) is int for value in asset_sizes)
            else None
        )
        return DataStatus(
            kind=self.kind,
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
        )

    def rows(self, *, cache_root=None, verify_files=False) -> tuple[DataStatus, ...]:
        from .genomes import pinned_genomes

        return tuple(self._row(genome, cache_root, verify_files) for genome in pinned_genomes())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        from .genomes import resolve_pinned_genome

        genome = resolve_pinned_genome(name)
        row = self._row(genome, cache_root, verify_files)
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        return row

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        from .genomes import resolve_pinned_genome

        genome = resolve_pinned_genome(name)
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        result = row.to_dict()
        result["source"] = asdict(genome)
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        from .genomes import GENOME_ROLES

        cleaned = _reject_extra_options(self.kind, row.name, options, self._ALLOWED_OPTIONS)
        resolved = resolve_source_paths(source_paths, GENOME_ROLES)
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        from .genomes import install_genome, resolve_pinned_genome

        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, cleaned = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        genome = resolve_pinned_genome(name)
        return install_genome(
            genome.name,
            cache_root=cache_root,
            force=force,
            source_paths=resolved_paths,
            genome_fasta_path_or_url=cleaned.get("genome_fasta_path_or_url"),
            download_options=dict(download_options or {}) or None,
        )


class ChainResourceKind:
    """Pinned liftover chain files between assemblies."""

    kind = "chain"

    def _find(self, name: str):
        from .liftover import pinned_chains

        for chain in pinned_chains():
            if chain.name == name:
                return chain
        raise ValueError(f"Unknown chain file {name!r}; see 'orfdata data list --kind chain'")

    def _row(self, chain, cache_root, verify_files) -> DataStatus:
        from .liftover import chain_registry

        inspection = chain_registry(chain, cache_root).inspect(
            chain.name, "1", verify_files=verify_files,
        )
        status = normalize_status(inspection.status)
        size = chain.size if type(chain.size) is int else None
        inspected_path = getattr(inspection, "path", None)
        return DataStatus(
            kind=self.kind,
            name=chain.name,
            version="1",
            title=f"{chain.source_build} to {chain.target_build} chain ({chain.provider})",
            species=(chain.species,),
            taxon_ids=(species_taxon(chain.species),),
            contexts=(),
            status=status,
            path=None if inspected_path is None else str(inspected_path),
            origin="download" if status == "available" else None,
            acquire="download",
            command=None,
            restricted=chain.restricted,
            default=True,
            download_size=size,
            installed_size=size,
            datacache_status=inspection.status,
        )

    def rows(self, *, cache_root=None, verify_files=False) -> tuple[DataStatus, ...]:
        from .liftover import pinned_chains

        return tuple(self._row(chain, cache_root, verify_files) for chain in pinned_chains())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        chain = self._find(name)
        row = self._row(chain, cache_root, verify_files)
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        return row

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        chain = self._find(name)
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        result = row.to_dict()
        result["source"] = asdict(chain)
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        cleaned = _reject_extra_options(self.kind, row.name, options)
        chain = self._find(row.name)
        resolved = resolve_source_paths(source_paths, (chain.filename,))
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        from .liftover import download_chain

        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, _ = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        chain = self._find(name)
        source_path = None if resolved_paths is None else resolved_paths[chain.filename]
        return download_chain(
            chain.source_build,
            chain.target_build,
            species=chain.species,
            cache_root=cache_root,
            force=force,
            accept_terms=accept_terms,
            source_path=source_path,
            download_options=dict(download_options or {}) or None,
        )


class RecordsResourceKind:
    """Pinned source-record archives used to reproduce viral collections offline."""

    kind = "records"

    def _find(self, name: str):
        from .source_records import pinned_records

        for pin in pinned_records():
            if pin.name == name:
                return pin
        raise ValueError(f"Unknown source records {name!r}; see 'orfdata data list --kind records'")

    def _row(self, pin, cache_root, verify_files) -> DataStatus:
        from .source_records import records_registry

        inspection = records_registry(pin, cache_root).inspect(
            pin.name, "1", verify_files=verify_files,
        )
        status = normalize_status(inspection.status)
        size = pin.size if type(pin.size) is int else None
        return DataStatus(
            kind=self.kind,
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
        )

    def rows(self, *, cache_root=None, verify_files=False) -> tuple[DataStatus, ...]:
        from .source_records import pinned_records

        return tuple(self._row(pin, cache_root, verify_files) for pin in pinned_records())

    def inspect(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> DataStatus:
        pin = self._find(name)
        row = self._row(pin, cache_root, verify_files)
        if version is not None and version != row.version:
            raise ValueError(f"{name!r} has version {row.version!r}, not {version!r}")
        return row

    def info(self, name, version=None, *, cache_root=None, verify_files=False, options=None) -> dict[str, object]:
        pin = self._find(name)
        row = self.inspect(name, version, cache_root=cache_root, verify_files=verify_files, options=options)
        result = row.to_dict()
        result["source"] = asdict(pin)
        return result

    def validate_options(self, row, *, cache_root=None, source_paths=None, options=None):
        cleaned = _reject_extra_options(self.kind, row.name, options)
        pin = self._find(row.name)
        resolved = resolve_source_paths(source_paths, (pin.filename,))
        return resolved, cleaned

    def acquire(
        self,
        name,
        version=None,
        *,
        cache_root=None,
        force=False,
        source_paths=None,
        accept_terms=False,
        download_options=None,
        options=None,
    ):
        from .source_records import install_records

        row = self.inspect(name, version, cache_root=cache_root, verify_files=False, options=options)
        resolved_paths, _ = self.validate_options(
            row, cache_root=cache_root, source_paths=source_paths, options=options,
        )
        pin = self._find(name)
        source_path = None if resolved_paths is None else resolved_paths[pin.filename]
        return install_records(
            name,
            cache_root=cache_root,
            force=force,
            source_path=source_path,
            download_options=dict(download_options or {}) or None,
        )


class ResourceKindRegistry:
    """Ordered registry mapping resource-kind names to :class:`ResourceKind` handlers."""

    def __init__(self, handlers: Iterable[ResourceKind] = ()):
        self._handlers: dict[str, ResourceKind] = {}
        for handler in handlers:
            self.register(handler)

    def register(self, handler: ResourceKind) -> None:
        """Register or replace one :class:`ResourceKind` handler by its ``kind`` name."""
        kind = getattr(handler, "kind", None)
        if not isinstance(kind, str) or not kind.strip():
            raise ValueError("A ResourceKind handler requires a nonempty string 'kind'")
        self._handlers[kind] = handler

    def kinds(self) -> tuple[str, ...]:
        """Return every registered resource-kind name in registration order."""
        return tuple(self._handlers)

    def get(self, kind: str) -> ResourceKind:
        """Return the handler registered for ``kind``, raising ``ValueError`` if unknown."""
        handler = self._handlers.get(kind)
        if handler is None:
            raise ValueError(
                f"Unknown data kind {kind!r}; choose from {', '.join(self._handlers)}"
            )
        return handler

    def handlers(self) -> tuple[ResourceKind, ...]:
        """Return every registered handler in registration order."""
        return tuple(self._handlers.values())


RESOURCE_KINDS = ResourceKindRegistry(
    (
        ViralResourceKind(),
        InventoryResourceKind(),
        CatalogResourceKind(),
        DatasetResourceKind(),
        ReferenceResourceKind(),
        GenomeResourceKind(),
        ChainResourceKind(),
        RecordsResourceKind(),
    )
)


def registered_resource_kinds() -> tuple[str, ...]:
    """Return every registered ORFdata resource kind name."""
    return RESOURCE_KINDS.kinds()


def get_resource_kind(kind: str) -> ResourceKind:
    """Return the :class:`ResourceKind` handler registered for ``kind``."""
    return RESOURCE_KINDS.get(kind)


def register_resource_kind(handler: ResourceKind) -> None:
    """Register a :class:`ResourceKind` implementation in the global registry."""
    RESOURCE_KINDS.register(handler)
