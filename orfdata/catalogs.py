"""Explicit optional catalogue downloads and builds using OpenVax DataCache."""

import gzip
import json

import datacache

from .cache import installed_path
from .catalog import ORFCatalog, load_catalog
from .catalog_manifest import CatalogManifest
from .downloads import install_hint
from .errors import NotInstalledError
from .registry import ReleaseRegistry, build_store_path

CATALOG_ASSET = "catalog.json.gz"
CATALOG_BUILD_AREA = "catalog-imports"
_ASSET = CATALOG_ASSET


def _open(path, manifest):
    catalog = load_catalog(path)
    manifest.validate_catalog(catalog)
    return catalog


def _copy_catalog(inputs, output, manifest):
    """Wrap pinned portable data with release provenance, retaining native facts."""
    asset = manifest.processing.get("catalog_asset")
    if asset not in inputs:
        raise ValueError("Copy importer requires a declared catalog_asset")
    catalog = load_catalog(inputs[asset])
    metadata = catalog.metadata
    if "catalog_manifest" in metadata:
        raise ValueError("Copy importer expects an unmanaged source catalogue; cannot replace release provenance")
    metadata["catalog_manifest"] = manifest.to_dict()
    catalog = ORFCatalog(catalog.records, metadata)
    manifest.validate_catalog(catalog)
    with open(output, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
        handle.write(json.dumps(catalog.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def _builder(manifest, *, for_build=False):
    """The packaged importer a release names; import paths in metadata are never executed."""
    importer = manifest.processing.get("importer")
    if importer == "catalog-copy-v1":
        return lambda inputs, output: _copy_catalog(inputs, output, manifest)
    if importer == "openprot-v1":
        from .sources.openprot import import_orfs

        return lambda inputs, output: import_orfs(inputs, output, manifest)
    if importer == "ncorf-mammals-v1":
        from .sources.ncorf_mammals import import_orfs

        return lambda inputs, output: import_orfs(inputs, output, manifest)
    if importer in ("riborf-v1", "riborf-v2"):
        from .sources.riborf import import_orfs, validate_transcript_manifest

        if for_build:
            validate_transcript_manifest(manifest)
        return lambda inputs, output: import_orfs(inputs, output, manifest)
    if importer == "riborf-native-v3":
        from .sources.riborf_native import import_orfs

        return lambda inputs, output: import_orfs(inputs, output, manifest)
    if importer in ("riborf-genome-v1", "riborf-genome-v2"):
        from .sources.riborf import import_genome_orfs

        return lambda inputs, output: import_genome_orfs(inputs, output, manifest)
    raise ValueError("No packaged importer for this catalogue release")


class CatalogRegistry(ReleaseRegistry):
    """Discover offline, explicitly acquire, then open the existing ORFCatalog.

    Each definition names a pinned default_version, manifests, distribution
    terms and either one checksummed catalog.json.gz asset or, for a source
    that orfdata does not redistribute, no assets: such a release is only
    built locally from its pinned sources. No method except download or build
    may acquire data. Bundled load_catalog() remains independent.
    """

    RESOURCE_KIND = "catalogue"
    PACKAGED_PATH = "data/catalogs.json"
    CACHE_AREA = "catalogs"
    BUILD_AREA = CATALOG_BUILD_AREA
    OUTPUT_ASSET = CATALOG_ASSET
    MANIFEST_TYPE = CatalogManifest
    MANIFEST_ID_ATTR = "catalog_id"
    METADATA_KEYS = ("distribution", "assets")
    DEFAULT_VERSION_ERROR = "A catalogue default must be a pinned version"
    IDENTITY_ERROR = "Registry and manifest catalogue identities disagree"
    ACCESS_TERMS_ERROR = "Catalogue artifact access terms are required"
    REDISTRIBUTABLE_ERROR = "Catalogue artifact redistribution must be explicit"
    UNKNOWN_NAME_HINT = "see 'orfdata catalogs list' or CatalogRegistry().list()"
    BUILD_INVENTORY_ERROR = "Built catalogue source/output inventory disagrees"
    BUILD_INTEGRITY_ERROR = "Built catalogue source integrity pins disagree"
    BUILD_TRANSFORM_ERROR = "Built catalogue transform/manifest disagrees"

    def validate_entry(self, name, entry):
        """Validate catalogue-specific ``priority`` and ``optional`` flags."""
        if type(entry.get("priority", 0)) is not int:
            raise ValueError(f"The {name} catalogue's priority must be an integer, lower merging first")
        if type(entry.get("optional", False)) is not bool:
            raise ValueError(f"The {name} catalogue's optional flag must be true or false")

    def validate_release(self, name, version, release, manifest):
        """Validate raw ``catalog.json.gz`` asset pins or buildable source importer."""
        if release["assets"]:
            if set(release["assets"]) != {_ASSET} or release["assets"][_ASSET].get("decompress", False):
                raise ValueError("A catalogue release requires raw catalog.json.gz bytes")
        else:
            _builder(manifest)  # A release with nothing to download must be buildable.

    def unpublished_error(self, name, version):
        """Return the build guidance error message for an unpublished catalogue release."""
        return (
            f"The {name} {version} ORF catalogue has no published artifact; build it from its "
            f"pinned sources with CatalogRegistry().build({name!r}, {version!r}) or "
            f"'orfdata catalogs build {name} --version {version}'"
        )

    def list(self, *, species=None, taxon_id=None):
        """Return matching versioned manifests without downloading data."""
        return self.list_manifests(species=species, taxon_id=taxon_id)

    def priority(self, name):
        """Where the catalogue comes when merged: lower first, so its IDs name shared ORFs (default 100)."""
        self.require_name(name)
        return self._definitions[name].get("priority", 100)

    def optional(self, name):
        """Whether the catalogue is used only when asked for (``catalogs="all"`` or its name), not by default."""
        self.require_name(name)
        return self._definitions[name].get("optional", False)

    def inspect(self, name, version=None, *, verify_files=True):
        """The installed artifact, or for a release that is only built, its local build."""
        if not self.published(name, version):
            return self.inspect_build(name, version, verify_files=verify_files)
        return self.data_registry().inspect(name, self.require_published(name, version), verify_files=verify_files)

    def download(self, name, version=None, *, force=False, source_paths=None, **options):
        """Explicitly install a pinned artifact; accept verified local overrides."""
        return self.data_registry().download(
            name, self.require_published(name, version), force=force, source_paths=source_paths, **options
        )

    def open(self, name, version=None, *, verify_files=True):
        """Validate installed bytes and content, then return an offline ORFCatalog.

        A release that is only built locally opens its build.
        """
        if not self.published(name, version):
            return self.open_build(name, version, verify_files=verify_files)
        manifest = self.info(name, version)
        version = self.require_published(name, version)

        def hint():
            download = install_hint(self, f"download({name!r}, {version!r})", name, "catalog", version)
            if self.inspect_build(name, version, verify_files=False).status == "available":
                return (
                    f"The {name} {version} ORF catalogue's published artifact is not downloaded, but a local "
                    f"build is installed: open it with open_build({name!r}, {version!r}) or "
                    f"'orfdata catalogs query {name} --built', or download the artifact with {download}"
                )
            return f"The {name} {version} ORF catalogue is not installed; download it with {download}"

        return _open(
            installed_path(self.data_registry(), name, version, hint, asset=_ASSET, verify_files=verify_files),
            manifest,
        )

    def build(self, name, version=None, *, force=False, source_paths=None, download_options=None):
        """Explicitly materialize an allowlisted catalogue source adapter."""
        manifest = self.info(name, version)
        return build_catalog(
            manifest,
            _builder(manifest, for_build=True),
            cache_root=self._cache_root,
            source_paths=source_paths,
            force=force,
            download_options=download_options,
        )

    def inspect_build(self, name, version=None, *, verify_files=True):
        """Inspect build receipts offline; their hashes are not released-artifact pins."""
        return super().inspect_build(name, version, verify_files=verify_files)

    def open_build(self, name, version=None, *, verify_files=True):
        """Open the local source build explicitly; never fall back to a download."""
        inspected = self.inspect_build(name, version, verify_files=verify_files)
        if inspected.status == "missing":
            version = self.resolve_version(name, version)
            raise NotInstalledError(
                f"The {name} {version} ORF catalogue has no local build; build it from its pinned sources with "
                + install_hint(
                    self,
                    f"build({name!r}, {version!r})",
                    name,
                    "catalog",
                    version,
                    build=self.published(name, version),
                )
            )
        if inspected.status != "available":
            raise datacache.FileValidationError(inspected.path, inspected.status) from inspected.error
        return _open(inspected.files[_ASSET].path, self.info(name, version))


def build_catalog(manifest, builder, *, cache_root=None, source_paths=None, force=False, download_options=None):
    """Materialize pinned sources with builder(inputs, output_path), atomically.

    The output is a portable ORF catalogue with its exact CatalogManifest in
    metadata['catalog_manifest']. Source and output hashes live in DataCache's
    receipt. Building does not register or publish a released artifact.
    """
    if not isinstance(manifest, CatalogManifest):
        raise TypeError("manifest must be a CatalogManifest")
    source_paths = {} if source_paths is None else dict(source_paths)
    if set(source_paths) - {s.asset_id for s in manifest.sources}:
        raise ValueError("Unknown local catalogue source asset override")
    sources = {}
    for asset in manifest.sources:
        source = dict(sha256=asset.sha256, size=asset.size)
        if asset.asset_id in source_paths:
            source["path"] = str(source_paths[asset.asset_id])
        else:
            source["url"] = asset.url
        sources[asset.asset_id] = source

    def build(inputs, outputs):
        builder(inputs, outputs[_ASSET])
        _open(outputs[_ASSET], manifest)

    paths = datacache.materialize(
        build_store_path(CATALOG_BUILD_AREA, manifest.catalog_id, manifest, cache_root),
        sources,
        transform=dict(version=manifest.schema, options=manifest.to_dict()),
        outputs={_ASSET: {}},
        builder=build,
        force=force,
        download_options=download_options,
    )
    return _open(paths[_ASSET], manifest)
