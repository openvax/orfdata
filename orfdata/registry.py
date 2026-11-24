"""Shared registry base and build-receipt helpers for versioned ORFdata releases.

Both :class:`~orfdata.catalogs.CatalogRegistry` and
:class:`~orfdata.datasets.DatasetRegistry` manage JSON-defined pinned releases
that are either downloaded as published artifacts or built locally from pinned
sources. This module holds the four shared pieces: version resolution and
release validation, species/taxon filtering, published-artifact DataCache
registry construction, and build-receipt verification.
"""

import json
from dataclasses import replace
from importlib.resources import files
from pathlib import Path

import datacache

from .cache import content_digest, get_data_dir
from .reference import species_taxon


def load_registry(registry_cls, registry_path=None, *, cache_root=None, label=None):
    """Construct a :class:`ReleaseRegistry`, optionally reading JSON definitions from ``registry_path``.

    JSON syntax and schema errors in a user-supplied ``registry_path`` are
    reported as ``ValueError`` naming the file, while unexpected errors in the
    packaged registry keep their traceback.

    Parameters
    ----------
    registry_cls : type
        A :class:`ReleaseRegistry` subclass such as
        :class:`~orfdata.catalogs.CatalogRegistry` or
        :class:`~orfdata.datasets.DatasetRegistry`.
    registry_path : str or os.PathLike, optional
        Path to a JSON registry definition file. When ``None``, the packaged
        registry is loaded.
    cache_root : str or os.PathLike, optional
        Explicit cache root override forwarded to ``registry_cls``.
    label : str, optional
        Human-readable registry kind used in error messages (defaults to
        ``registry_cls.RESOURCE_KIND``).

    Returns
    -------
    ReleaseRegistry
        Initialized registry instance.

    Raises
    ------
    ValueError
        If ``registry_path`` contains invalid JSON or malformed definitions.
    """
    kind_label = label or getattr(registry_cls, "RESOURCE_KIND", "release")
    definitions = None
    if registry_path:
        try:
            definitions = json.loads(Path(registry_path).read_text())
        except ValueError as error:
            raise ValueError(f"Invalid {kind_label} registry {registry_path}: {error}") from error
    try:
        return registry_cls(definitions, cache_root=cache_root)
    except (KeyError, TypeError, AttributeError) as error:
        if not registry_path:  # The packaged registry: a bug, which keeps its traceback.
            raise
        raise ValueError(
            f"Invalid {kind_label} registry {registry_path}: {type(error).__name__}: {error}"
        ) from error


def resolve_species_taxon(species=None, taxon_id=None):
    """Resolve an optional species alias and optional NCBI taxon ID into one taxon ID.

    Parameters
    ----------
    species : str, optional
        Common or scientific species alias recognized by
        :func:`~orfdata.reference.species_taxon`.
    taxon_id : int, optional
        Explicit NCBI taxonomy ID.

    Returns
    -------
    int or None
        The resolved NCBI taxonomy ID, or ``None`` if neither argument was given.

    Raises
    ------
    ValueError
        If ``species`` is not a recognized alias, or if both ``species`` and
        ``taxon_id`` are supplied and disagree.
    """
    if species is not None:
        found = species_taxon(species)
        if found is None:
            raise ValueError("Unknown species alias; supply its explicit taxon_id")
        if taxon_id is not None and taxon_id != found:
            raise ValueError("Species and taxon_id disagree")
        return found
    return taxon_id


def build_store_path(build_area, resource_id, manifest, cache_root=None):
    """Return the local DataCache materialization directory for a built release.

    Parameters
    ----------
    build_area : str
        Subdirectory under :func:`~orfdata.cache.get_data_dir`, such as
        ``"catalog-imports"`` or ``"imports"``.
    resource_id : str
        Catalogue or dataset identifier.
    manifest : object
        Manifest object providing ``version`` and ``to_dict()``.
    cache_root : str or os.PathLike, optional
        Explicit cache root override.

    Returns
    -------
    pathlib.Path
        Content-addressed materialization directory for ``manifest``.
    """
    return (
        get_data_dir(cache_root)
        / build_area
        / resource_id
        / manifest.version
        / content_digest(manifest.to_dict())
    )


def build_published_registry(
    definitions,
    cache_area,
    cache_root=None,
    *,
    include_release=None,
    store_path=None,
):
    """Build a DataCache :class:`~datacache.VersionedDatasetRegistry` for published releases.

    Only releases with non-empty ``assets`` (and accepted by ``include_release``,
    when provided) are registered. If an entry's ``default_version`` is not
    among its published versions, the lexicographically smallest published
    version is used as DataCache's default.

    Parameters
    ----------
    definitions : dict
        Registry definitions mapping resource names to ``default_version`` and
        ``versions`` dictionaries.
    cache_area : str
        Cache subdirectory name under :func:`~orfdata.cache.get_data_dir`,
        such as ``"catalogs"`` or ``"datasets"``.
    cache_root : str or os.PathLike, optional
        Explicit cache root override used by the default ``store_path``.
    include_release : callable, optional
        Predicate ``include_release(release) -> bool`` filtering which
        asset-bearing releases belong in the DataCache registry.
    store_path : callable, optional
        Optional ``store_path(name, version) -> Path`` closure overriding the
        default ``<data_dir>/<cache_area>/<name>/<version>/<release_digest>``
        layout.

    Returns
    -------
    datacache.VersionedDatasetRegistry
        Registry covering the selected published releases.
    """
    published_defs = {}
    for name, entry in definitions.items():
        published = {
            version: release["assets"]
            for version, release in entry["versions"].items()
            if release["assets"] and (include_release is None or include_release(release))
        }
        if published:  # Every call names its version; the default only has to be a published one.
            default = entry["default_version"]
            published_defs[name] = dict(
                default_version=default if default in published else min(published),
                versions=published,
            )

    if store_path is None:
        def store_path(name, version):
            release = definitions[name]["versions"][version]
            return get_data_dir(cache_root) / cache_area / name / version / content_digest(release)

    return datacache.VersionedDatasetRegistry(published_defs, store_path=store_path)


def inspect_build_receipt(
    store_path,
    manifest,
    output_asset,
    *,
    verify_files=True,
    inventory_error="Built source/output inventory disagrees",
    integrity_error="Built source integrity pins disagree",
    transform_error="Built transform/manifest disagrees",
):
    """Inspect a local DataCache materialization receipt against a release manifest.

    Verifies that the recorded source and output inventory, each source's
    SHA-256, byte size, and ``verified`` flag, and the recorded transform match
    ``manifest``.

    Parameters
    ----------
    store_path : str or os.PathLike
        Directory passed to :func:`datacache.inspect_materialization`.
    manifest : object
        Release manifest providing ``sources``, ``schema``, and ``to_dict()``.
    output_asset : str
        Expected single output asset filename (e.g. ``"catalog.json.gz"`` or
        ``"dataset.sqlite"``).
    verify_files : bool, optional
        Whether DataCache should re-verify output file bytes against the receipt.
    inventory_error : str, optional
        Error message used when source or output keys disagree with ``manifest``.
    integrity_error : str, optional
        Error message used when a source's SHA-256, size, or verified flag
        disagrees with ``manifest``.
    transform_error : str, optional
        Error message used when the recorded transform disagrees with ``manifest``.

    Returns
    -------
    datacache.MaterializationInspection
        The DataCache inspection result, downgraded to ``status="invalid"`` with
        ``verified=False`` if receipt validation fails.
    """
    inspected = datacache.inspect_materialization(store_path, verify_files=verify_files)
    if inspected.status != "available":
        return inspected
    try:
        expected = {source.asset_id: source for source in manifest.sources}
        if set(inspected.sources) != set(expected) or set(inspected.files) != {output_asset}:
            raise ValueError(inventory_error)
        for asset_id, source in expected.items():
            recorded = inspected.sources[asset_id]
            if (
                recorded["sha256"] != source.sha256
                or recorded["size"] != source.size
                or recorded["verified"] is not True
            ):
                raise ValueError(integrity_error)
        if inspected.transform != dict(version=manifest.schema, options=manifest.to_dict()):
            raise ValueError(transform_error)
    except (KeyError, TypeError, ValueError) as exc:
        return replace(inspected, status="invalid", verified=False, error=exc)
    return inspected


class ReleaseRegistry:
    """Base registry for JSON-defined pinned releases backed by DataCache.

    Subclasses configure cache areas, asset names, manifest types, and
    resource-specific error messages via class attributes, and override
    :meth:`validate_entry`, :meth:`validate_release`,
    :meth:`include_in_data_registry`, :meth:`manifest_matches_taxon`, and
    :meth:`unpublished_error` for behavior unique to that resource kind.
    """

    RESOURCE_KIND = "release"
    PACKAGED_PATH = ""
    CACHE_AREA = ""
    BUILD_AREA = ""
    OUTPUT_ASSET = ""
    MANIFEST_TYPE = None
    MANIFEST_ID_ATTR = ""
    METADATA_KEYS = ("distribution", "assets")
    DEFAULT_VERSION_ERROR = "A release default must be a pinned version"
    IDENTITY_ERROR = "Registry and manifest identities disagree"
    ACCESS_TERMS_ERROR = "Artifact access terms are required"
    REDISTRIBUTABLE_ERROR = "Artifact redistribution must be explicit"
    UNKNOWN_NAME_HINT = ""
    BUILD_INVENTORY_ERROR = "Built source/output inventory disagrees"
    BUILD_INTEGRITY_ERROR = "Built source integrity pins disagree"
    BUILD_TRANSFORM_ERROR = "Built transform/manifest disagrees"

    def __init__(self, definitions=None, *, cache_root=None):
        self._packaged = definitions is None  # What 'orfdata data fetch' reads.
        if definitions is None:
            definitions = json.loads(files("orfdata").joinpath(self.PACKAGED_PATH).read_text())
        # Own a JSON copy so callers cannot mutate identities after validation.
        self._definitions = json.loads(json.dumps(definitions, allow_nan=False))
        self._manifests = {}
        self._cache_root = cache_root
        self._registry = None
        for name, entry in self._definitions.items():
            if entry["default_version"] not in entry["versions"]:
                raise ValueError(self.DEFAULT_VERSION_ERROR)
            self.validate_entry(name, entry)
            for version, release in entry["versions"].items():
                manifest = self.MANIFEST_TYPE.from_dict(release["manifest"])
                if (getattr(manifest, self.MANIFEST_ID_ATTR), manifest.version) != (name, version):
                    raise ValueError(self.IDENTITY_ERROR)
                self.validate_release(name, version, release, manifest)
                distribution = release["distribution"]
                if (
                    not isinstance(distribution.get("access_terms"), str)
                    or not distribution["access_terms"].strip()
                ):
                    raise ValueError(self.ACCESS_TERMS_ERROR)
                if type(distribution.get("redistributable")) is not bool:
                    raise ValueError(self.REDISTRIBUTABLE_ERROR)
                self._manifests[name, version] = manifest

    def validate_entry(self, name, entry):
        """Validate resource-level metadata beyond ``default_version`` and ``versions``."""

    def validate_release(self, name, version, release, manifest):
        """Validate release-level assets, archives, and buildability."""

    def include_in_data_registry(self, release):
        """Return whether an asset-bearing release is managed directly by DataCache."""
        return True

    def manifest_matches_taxon(self, manifest, taxon_id):
        """Return whether ``manifest`` matches ``taxon_id``."""
        return manifest.taxon_id == taxon_id

    def unpublished_error(self, name, version):
        """Return the error message raised when a build-only release is downloaded."""
        raise NotImplementedError

    def require_name(self, name):
        """Verify that ``name`` is registered, raising ``ValueError`` otherwise.

        Parameters
        ----------
        name : str
            Registered resource name.

        Raises
        ------
        ValueError
            If ``name`` is not present in this registry.
        """
        if name not in self._definitions:
            raise ValueError(f"Unknown {self.RESOURCE_KIND} {name!r}; {self.UNKNOWN_NAME_HINT}")

    def resolve_version(self, name, version=None):
        """Resolve ``version`` (or the pinned default) for a registered resource.

        Parameters
        ----------
        name : str
            Registered resource name.
        version : str, optional
            Explicit release version, or ``None`` for the pinned default.

        Returns
        -------
        str
            Resolved pinned version string.

        Raises
        ------
        ValueError
            If ``name`` or ``(name, version)`` is not registered.
        """
        self.require_name(name)
        version = self._definitions[name]["default_version"] if version is None else version
        if (name, version) not in self._manifests:
            raise ValueError(f"Unknown {self.RESOURCE_KIND} version {name!r} {version!r}")
        return version

    def list_manifests(self, *, species=None, taxon_id=None, predicate=None):
        """Return matching manifests in sorted ``(name, version)`` order without I/O.

        Parameters
        ----------
        species : str, optional
            Species alias filter.
        taxon_id : int, optional
            Explicit NCBI taxonomy ID filter.
        predicate : callable, optional
            Additional filter ``predicate(manifest) -> bool``.

        Returns
        -------
        tuple
            Matching manifest instances.
        """
        taxon_id = resolve_species_taxon(species=species, taxon_id=taxon_id)
        return tuple(
            manifest
            for _, manifest in sorted(self._manifests.items())
            if (taxon_id is None or self.manifest_matches_taxon(manifest, taxon_id))
            and (predicate is None or predicate(manifest))
        )

    def info(self, name, version=None):
        """Return immutable provenance and coverage for a pinned release."""
        return self._manifests[name, self.resolve_version(name, version)]

    def download_metadata(self, name, version=None):
        """Return a JSON copy of artifact access terms, locations, and integrity pins."""
        version = self.resolve_version(name, version)
        release = self._definitions[name]["versions"][version]
        return json.loads(
            json.dumps({key: release[key] for key in self.METADATA_KEYS if key in release})
        )

    def published(self, name, version=None):
        """Return whether the release has a published artifact to download; otherwise it is built locally."""
        version = self.resolve_version(name, version)
        return bool(self._definitions[name]["versions"][version]["assets"])

    def require_published(self, name, version=None):
        """Resolve ``version`` and verify that the release has a published artifact.

        Parameters
        ----------
        name : str
            Registered resource name.
        version : str, optional
            Explicit release version, or ``None`` for the pinned default.

        Returns
        -------
        str
            Resolved pinned version string.

        Raises
        ------
        ValueError
            If the release has no published artifact.
        """
        version = self.resolve_version(name, version)
        if not self.published(name, version):
            raise ValueError(self.unpublished_error(name, version))
        return version

    def release_store_path(self, name, version):
        """Return the content-addressed cache directory for a published release."""
        release = self._definitions[name]["versions"][version]
        return get_data_dir(self._cache_root) / self.CACHE_AREA / name / version / content_digest(release)

    def data_registry(self):
        """Return the lazily constructed DataCache registry for published releases."""
        if self._registry is None:
            self._registry = build_published_registry(
                self._definitions,
                self.CACHE_AREA,
                self._cache_root,
                include_release=self.include_in_data_registry,
                store_path=self.release_store_path,
            )
        return self._registry

    def build_store_path(self, manifest):
        """Return the local materialization directory for a built ``manifest``."""
        return build_store_path(
            self.BUILD_AREA,
            getattr(manifest, self.MANIFEST_ID_ATTR),
            manifest,
            self._cache_root,
        )

    def inspect_build(self, name, version=None, *, verify_files=True):
        """Inspect a local source build receipt offline against the registered manifest."""
        manifest = self.info(name, version)
        return inspect_build_receipt(
            self.build_store_path(manifest),
            manifest,
            self.OUTPUT_ASSET,
            verify_files=verify_files,
            inventory_error=self.BUILD_INVENTORY_ERROR,
            integrity_error=self.BUILD_INTEGRITY_ERROR,
            transform_error=self.BUILD_TRANSFORM_ERROR,
        )
