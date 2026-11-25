"""Explicit, checksummed dataset downloads using OpenVax DataCache conventions."""

import datacache

from . import artifact_archives
from .cache import content_digest, get_data_dir, installed_path
from .dataset import QuantitativeDataset
from .downloads import install_hint
from .errors import NotInstalledError
from .measurements import DatasetManifest
from .registry import ReleaseRegistry, build_store_path, inspect_build_receipt
from .source_store import SourceStore

DATASET_ASSET = "dataset.sqlite"
DATASET_BUILD_AREA = "imports"
DATASET_BUILD_INVENTORY_ERROR = "Built source/output inventory disagrees with the registered release"
DATASET_BUILD_INTEGRITY_ERROR = "Built source integrity pins disagree with the registered release"
DATASET_BUILD_TRANSFORM_ERROR = "Built transform/manifest disagrees with the registered release"


class DatasetRegistry(ReleaseRegistry):
    """Discover manifests offline; explicitly download before opening data.

    Definitions map names to a pinned default_version and versions. Each version
    contains a manifest, distribution access terms, and checksummed assets. A
    version with no assets is built locally from its pinned sources, and must
    name a packaged importer. No dataset is downloaded on construction,
    discovery, inspection or open.
    """

    RESOURCE_KIND = "dataset"
    PACKAGED_PATH = "data/datasets.json"
    CACHE_AREA = "datasets"
    BUILD_AREA = DATASET_BUILD_AREA
    OUTPUT_ASSET = DATASET_ASSET
    MANIFEST_TYPE = DatasetManifest
    MANIFEST_ID_ATTR = "dataset_id"
    METADATA_KEYS = ("distribution", "assets", "archive")
    DEFAULT_VERSION_ERROR = "A dataset default must be a pinned version"
    IDENTITY_ERROR = "Registry and manifest dataset identities disagree"
    ACCESS_TERMS_ERROR = "Dataset artifact access terms are required"
    REDISTRIBUTABLE_ERROR = "Dataset artifact redistribution must be declared explicitly"
    UNKNOWN_NAME_HINT = "use DatasetRegistry.list()"
    BUILD_INVENTORY_ERROR = DATASET_BUILD_INVENTORY_ERROR
    BUILD_INTEGRITY_ERROR = DATASET_BUILD_INTEGRITY_ERROR
    BUILD_TRANSFORM_ERROR = DATASET_BUILD_TRANSFORM_ERROR

    def validate_release(self, name, version, release, manifest):
        """Validate ``dataset.sqlite`` asset pins, optional archive config, or buildable importer."""
        if release["assets"]:
            if set(release["assets"]) != {self.OUTPUT_ASSET}:
                raise ValueError("A quantitative release must contain dataset.sqlite")
            if "archive" in release:
                artifact_archives.archive_options(
                    release["archive"], self.OUTPUT_ASSET, release["assets"][self.OUTPUT_ASSET]
                )
        else:
            if "archive" in release:
                raise ValueError("An archive requires its extracted dataset.sqlite integrity pin")
            from .sources import importer

            importer(manifest.processing.get("importer"))  # A release with nothing to download must be buildable.

    def include_in_data_registry(self, release):
        """Exclude multi-part archived releases, which are managed via :mod:`orfdata.artifact_archives`."""
        return "archive" not in release

    def manifest_matches_taxon(self, manifest, taxon_id):
        """Return whether any declared reference context of ``manifest`` matches ``taxon_id``."""
        return any(reference.taxon_id == taxon_id for reference in manifest.references)

    def unpublished_error(self, name, version):
        """Return the build guidance message for an unpublished dataset release."""
        release = self._definitions[name]["versions"][version]
        reason = "" if release["distribution"]["redistributable"] else ", because its source terms don't allow redistribution"
        return (
            f"The {name} {version} dataset has no published artifact{reason}; build it from its pinned "
            f"sources with " + install_hint(self, f"build({name!r}, {version!r})", name, "dataset", version)
        )

    def list(self, *, species=None, taxon_id=None, assay=None):
        """Return matching manifests for all registered versions, without I/O."""
        predicate = None if assay is None else (lambda manifest: assay in manifest.assay_kinds)
        return self.list_manifests(species=species, taxon_id=taxon_id, predicate=predicate)

    def search_references(self, name, version=None):
        """Discover declared native search FASTAs without downloading any files."""
        from .sources import search_references

        return search_references(self.info(name, version))

    def _search_reference(self, name, version, reference, asset_id):
        key = getattr(reference, "key", reference)
        choices = [
            spec
            for spec in self.search_references(name, version)
            if (key is None or spec.reference.key == key)
            and (asset_id is None or spec.source.asset_id == asset_id)
        ]
        if not choices:
            raise ValueError("No declared native search reference matches this dataset/context/asset")
        if len(choices) != 1:
            raise ValueError("Select an exact reference context and source asset")
        return choices[0]

    def inspect_reference(self, name, version=None, *, reference=None, asset_id=None, verify_files=True):
        """Inspect the separate reference cache offline; never repair or download."""
        from .search_reference import reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        return reference_registry(spec, self._cache_root).inspect(
            spec.dataset_id, spec.dataset_version, verify_files=verify_files
        )

    def download_reference(
        self,
        name,
        version=None,
        *,
        reference=None,
        asset_id=None,
        force=False,
        source_path=None,
        **options,
    ):
        """Explicitly acquire a complete pinned search FASTA without rebuilding abundance."""
        from .search_reference import reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        overrides = {spec.source.asset_id: source_path} if source_path is not None else None
        return reference_registry(spec, self._cache_root).download(
            spec.dataset_id, spec.dataset_version, force=force, source_paths=overrides, **options
        )

    def open_reference(self, name, version=None, *, reference=None, asset_id=None):
        """Open and validate an installed native search set, without network access."""
        from .search_reference import SearchReference, reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        asset, version = spec.source.asset_id, spec.dataset_version
        hint = (
            f"The {name} {version} search reference {asset} is not installed; download it with "
            + install_hint(
                self,
                f"download_reference({name!r}, {version!r}, asset_id={asset!r})",
                f"{name}:{asset}",
                "reference",
                version,
            )
        )
        path = installed_path(reference_registry(spec, self._cache_root), spec.dataset_id, version, hint)
        return SearchReference(path, spec)

    def _not_installed(self, name, version):
        return (
            f"The {name} {version} dataset is not installed; download it with "
            + install_hint(self, f"download({name!r}, {version!r})", name, "dataset", version)
        )

    def _not_built(self, name, version):
        version = self.resolve_version(name, version)
        return NotInstalledError(
            f"The {name} {version} dataset has no local build; build it from its pinned sources with "
            + install_hint(
                self,
                f"build({name!r}, {version!r})",
                name,
                "dataset",
                version,
                build=self.published(name, version),
            )
        )

    def inspect(self, name, version=None, *, verify_files=True):
        """Inspect local bytes offline; never download or repair a bad cache.

        A release that is only built locally reports its build.
        """
        if not self.published(name, version):
            return self.inspect_build(name, version, verify_files=verify_files)
        version = self.require_published(name, version)
        release = self._definitions[name]["versions"][version]
        if "archive" in release:
            return artifact_archives.inspect_archive(
                self.release_store_path(name, version),
                release["archive"],
                self.OUTPUT_ASSET,
                release["assets"][self.OUTPUT_ASSET],
                verify_files=verify_files,
            )
        return self.data_registry().inspect(name, version, verify_files=verify_files)

    def download(self, name, version=None, *, force=False, source_paths=None, **options):
        """Explicitly install verified files; forward timeout/retry/progress options.

        source_paths accepts user-downloaded files while retaining declared URLs
        and hashes. force installs a fresh bundle without deleting readers' files.
        For an archive, source_paths["dataset.sqlite"] is the ordered part paths.
        """
        version = self.require_published(name, version)
        release = self._definitions[name]["versions"][version]
        if "archive" in release:
            if source_paths is not None and set(source_paths) != {self.OUTPUT_ASSET}:
                raise ValueError("Archive overrides require dataset.sqlite's ordered part paths")

            def validate(path):
                if QuantitativeDataset(path).manifest.to_dict() != self.info(name, version).to_dict():
                    raise ValueError("Installed dataset manifest disagrees with registered release")

            path = artifact_archives.install_archive(
                self.release_store_path(name, version),
                release["archive"],
                self.OUTPUT_ASSET,
                release["assets"][self.OUTPUT_ASSET],
                source_paths=None if source_paths is None else source_paths[self.OUTPUT_ASSET],
                force=force,
                download_options=options,
                validator=validate,
            )
            return {self.OUTPUT_ASSET: path}
        asset_spec = release["assets"][self.OUTPUT_ASSET]
        dist = release.get("distribution", {})
        if (
            source_paths is None
            and options.get("resume")
            and asset_spec.get("decompress")
            and dist.get("download_sha256")
            and dist.get("download_size")
        ):
            if not force:
                inspected = self.data_registry().inspect(name, version, verify_files=True)
                if inspected.status == "available":
                    return {
                        self.OUTPUT_ASSET: str(
                            self.data_registry().local_path(
                                name, version, asset=self.OUTPUT_ASSET, verify_files=False
                            )
                        )
                    }
            gz_path = SourceStore(self._cache_root, seed_from_pyensembl=False).acquire(
                asset_spec["url"],
                dist["download_sha256"],
                dist["download_size"],
                download_options=options,
            )
            local_options = {k: v for k, v in options.items() if k != "resume"}
            return self.data_registry().download(
                name,
                version,
                force=force,
                source_paths={self.OUTPUT_ASSET: str(gz_path)},
                **local_options,
            )
        return self.data_registry().download(
            name, version, force=force, source_paths=source_paths, **options
        )

    def build(self, name, version=None, *, force=False, source_paths=None, download_options=None):
        """Rebuild a registered release from pinned source files, explicitly."""
        from .sources import importer

        manifest = self.info(name, version)
        transform = importer(manifest.processing.get("importer"))
        return build_dataset(
            manifest,
            lambda inputs, output: transform(inputs, output, manifest),
            cache_root=self._cache_root,
            source_paths=source_paths,
            force=force,
            download_options=download_options,
        )

    def open(self, name, version=None, *, verify_files=True):
        """Open an installed dataset after integrity and identity checks, offline.

        A release that is only built locally opens its build.
        """
        if not self.published(name, version):
            if self.inspect_build(name, version, verify_files=False).status == "missing":
                version = self.resolve_version(name, version)
                raise NotInstalledError(
                    self.unpublished_error(name, version).replace(
                        "has no published artifact", "hasn't been built, and has no published artifact"
                    )
                )
            return self.open_build(name, version, verify_files=verify_files)
        version = self.require_published(name, version)
        if "archive" in self._definitions[name]["versions"][version]:
            inspected = self.inspect(name, version, verify_files=verify_files)
            if inspected.status == "missing":
                raise NotInstalledError(self._not_installed(name, version))
            if inspected.status != "available":
                raise datacache.FileValidationError(inspected.path, inspected.status) from inspected.error
            path = inspected.files[self.OUTPUT_ASSET].path
        else:
            path = installed_path(
                self.data_registry(),
                name,
                version,
                self._not_installed(name, version),
                asset=self.OUTPUT_ASSET,
                verify_files=verify_files,
            )
        dataset = QuantitativeDataset(path)
        if dataset.manifest.to_dict() != self.info(name, version).to_dict():
            raise ValueError("Installed dataset manifest disagrees with registered release")
        return dataset

    def inspect_build(self, name, version=None, *, verify_files=True):
        """Inspect a source build offline, checking source pins and build identity.

        Output bytes are checked against the build receipt. ``verified`` remains
        false: the receipt is not an independently trusted published output hash.
        Original source acquisition locations, including local overrides, remain
        in ``sources``. Inspection never reads or reacquires those source files.
        """
        return super().inspect_build(name, version, verify_files=verify_files)

    def open_build(self, name, version=None, *, verify_files=True):
        """Open a local source materialization explicitly, without artifact fallback."""
        inspected = self.inspect_build(name, version, verify_files=verify_files)
        if inspected.status == "missing":
            raise self._not_built(name, version)
        if inspected.status != "available":
            raise datacache.FileValidationError(inspected.path, inspected.status) from inspected.error
        dataset = QuantitativeDataset(inspected.files[self.OUTPUT_ASSET].path)
        if dataset.manifest.to_dict() != self.info(name, version).to_dict():
            raise ValueError("Built dataset manifest disagrees with the registered release")
        return dataset


def build_dataset(
    manifest,
    builder,
    *,
    cache_root=None,
    source_paths=None,
    force=False,
    download_options=None,
):
    """Materialize a source-pinned import locally using DataCache.

    builder(source_paths, output_path) must create a QuantitativeDataset with
    the exact supplied manifest. Source files are read-only inputs. Failed builds
    are never published; hashes, processing version/options and output checksums
    are retained in DataCache's receipt. Source redistribution is not implied.
    """
    if not isinstance(manifest, DatasetManifest):
        raise TypeError("manifest must be a DatasetManifest")
    source_paths = {} if source_paths is None else dict(source_paths)
    if set(source_paths) - {s.asset_id for s in manifest.sources}:
        raise ValueError("Unknown local source asset override")
    if any(value is None for value in source_paths.values()):
        raise ValueError("Local source overrides must name a file or ordered archive parts")
    for asset in manifest.sources:
        if "archive" in asset.metadata:
            artifact_archives.archive_options(
                asset.metadata["archive"],
                asset.asset_id,
                dict(sha256=asset.sha256, size=asset.size),
            )
    destination = build_store_path(DATASET_BUILD_AREA, manifest.dataset_id, manifest, cache_root)
    if not force:
        cached = inspect_build_receipt(
            destination,
            manifest,
            DATASET_ASSET,
            inventory_error=DATASET_BUILD_INVENTORY_ERROR,
            integrity_error=DATASET_BUILD_INTEGRITY_ERROR,
            transform_error=DATASET_BUILD_TRANSFORM_ERROR,
        )
        if cached.status == "available":
            dataset = QuantitativeDataset(cached.files[DATASET_ASSET].path)
            if dataset.manifest.to_dict() != manifest.to_dict():
                raise ValueError("Built dataset manifest disagrees with registered release")
            return dataset
        if cached.status != "missing":
            raise datacache.FileValidationError(cached.path, cached.status) from cached.error
    store = SourceStore(cache_root)
    sources = {}
    for asset in manifest.sources:
        source = dict(sha256=asset.sha256, size=asset.size)
        supplied = source_paths.get(asset.asset_id)
        if "archive" in asset.metadata and (supplied is None or isinstance(supplied, (list, tuple))):
            # Explicit builds may acquire a distributed source component. Queries
            # and open/inspect paths never call this acquisition branch.
            path = artifact_archives.install_archive(
                get_data_dir(cache_root) / "components" / content_digest(asset.to_dict()),
                asset.metadata["archive"],
                asset.asset_id,
                source,
                source_paths=supplied,
                download_options=download_options,
            )
            source["path"] = str(path)
        elif supplied is not None:
            source["path"] = str(supplied)
        else:
            cached_source = store.acquire(
                asset.url,
                asset.sha256,
                asset.size,
                download_options=download_options,
            )
            source["path"] = str(cached_source)
        sources[asset.asset_id] = source

    def build(inputs, outputs):
        builder(inputs, outputs[DATASET_ASSET])
        dataset = QuantitativeDataset(outputs[DATASET_ASSET])
        if dataset.manifest.to_dict() != manifest.to_dict():
            raise ValueError("Importer produced a different dataset manifest")

    identity = manifest.to_dict()
    paths = datacache.materialize(
        destination,
        sources,
        transform=dict(version=manifest.schema, options=identity),
        outputs={DATASET_ASSET: {}},
        builder=build,
        force=force,
        download_options=download_options,
    )
    return QuantitativeDataset(paths[DATASET_ASSET])
