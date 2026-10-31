"""Explicit, checksummed dataset downloads using OpenVax DataCache conventions."""

import hashlib
import json
from dataclasses import replace
from importlib.resources import files
from pathlib import Path

from . import artifact_archives
from .dataset import QuantitativeDataset
from .downloads import install_hint
from .errors import NotInstalledError
from .measurements import DatasetManifest
from .reference import species_taxon


def _datacache():
    import datacache

    return datacache


def get_data_dir(override=None):
    """Resolve ORFDATA_CACHE, OPENVAX_DATA_CACHE, or the platform OpenVax cache.

    ORFdata owns only the orfdata/ child; lookup does not create directories.
    """
    return Path(_datacache().get_cache_root(
        "openvax", "ORFDATA_CACHE", "OPENVAX_DATA_CACHE", override=override)) / "orfdata"


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _build_store(manifest, cache_root):
    return (get_data_dir(cache_root) / "imports" / manifest.dataset_id /
            manifest.version / _digest(manifest.to_dict()))


def _installed_path(registry, name, version, hint, **options):
    """A DataCache registry's local path, or a NotInstalledError saying ``hint`` when nothing is installed.

    DataCache raises FileNotFoundError only for a missing bundle; damaged files
    raise its own validation error, which passes through.
    """
    try:
        return registry.local_path(name, version, **options)
    except FileNotFoundError as error:
        raise NotInstalledError(hint) from error


class DatasetRegistry:
    """Discover manifests offline; explicitly download before opening data.

    Definitions map names to a pinned default_version and versions. Each version
    contains a manifest, distribution access terms, and checksummed assets. A
    version with no assets is built locally from its pinned sources, and must
    name a packaged importer. No dataset is downloaded on construction,
    discovery, inspection or open.
    """

    def __init__(self, definitions=None, *, cache_root=None):
        self._packaged = definitions is None  # What 'orfdata data fetch' reads.
        if definitions is None:
            definitions = json.loads(files("orfdata").joinpath("data/datasets.json").read_text())
        # Own a JSON copy so callers cannot mutate identities after validation.
        self._definitions = json.loads(json.dumps(definitions, allow_nan=False))
        self._manifests = {}
        self._cache_root = cache_root
        self._registry = None
        for name, entry in self._definitions.items():
            if entry["default_version"] not in entry["versions"]:
                raise ValueError("A dataset default must be a pinned version")
            for version, release in entry["versions"].items():
                manifest = DatasetManifest.from_dict(release["manifest"])
                if (manifest.dataset_id, manifest.version) != (name, version):
                    raise ValueError("Registry and manifest dataset identities disagree")
                if release["assets"]:
                    if set(release["assets"]) != {"dataset.sqlite"}:
                        raise ValueError("A quantitative release must contain dataset.sqlite")
                    if "archive" in release:
                        artifact_archives.archive_options(release["archive"], "dataset.sqlite", release["assets"]["dataset.sqlite"])
                else:
                    if "archive" in release:
                        raise ValueError("An archive requires its extracted dataset.sqlite integrity pin")
                    from .sources import importer

                    importer(manifest.processing.get("importer"))  # A release with nothing to download must be buildable.
                distribution = release["distribution"]
                if not isinstance(distribution.get("access_terms"), str) or not distribution["access_terms"].strip():
                    raise ValueError("Dataset artifact access terms are required")
                if type(distribution.get("redistributable")) is not bool:
                    raise ValueError("Dataset artifact redistribution must be declared explicitly")
                self._manifests[name, version] = manifest

    def _version(self, name, version):
        if name not in self._definitions:
            raise ValueError(f"Unknown dataset {name!r}; use DatasetRegistry.list()")
        version = self._definitions[name]["default_version"] if version is None else version
        if (name, version) not in self._manifests:
            raise ValueError(f"Unknown dataset version {name!r} {version!r}")
        return version

    def list(self, *, species=None, taxon_id=None, assay=None):
        """Return matching manifests for all registered versions, without I/O."""
        if species is not None:
            found = species_taxon(species)
            if found is None:
                raise ValueError("Unknown species alias; supply its explicit taxon_id")
            if taxon_id is not None and taxon_id != found:
                raise ValueError("Species and taxon_id disagree")
            taxon_id = found
        return tuple(manifest for _, manifest in sorted(self._manifests.items())
                     if (taxon_id is None or any(r.taxon_id == taxon_id for r in manifest.references))
                     and (assay is None or assay in manifest.assay_kinds))

    def info(self, name, version=None):
        """Return immutable provenance and coverage for a pinned release."""
        return self._manifests[name, self._version(name, version)]

    def download_metadata(self, name, version=None):
        """Return a JSON copy of artifact access terms, locations and integrity pins."""
        version = self._version(name, version)
        release = self._definitions[name]["versions"][version]
        return json.loads(json.dumps({key: release[key] for key in ("distribution", "assets", "archive") if key in release}))

    def search_references(self, name, version=None):
        """Discover declared native search FASTAs without downloading any files."""
        from .sources import search_references

        return search_references(self.info(name, version))

    def _search_reference(self, name, version, reference, asset_id):
        key = getattr(reference, "key", reference)
        choices = [spec for spec in self.search_references(name, version)
                   if (key is None or spec.reference.key == key)
                   and (asset_id is None or spec.source.asset_id == asset_id)]
        if not choices:
            raise ValueError("No declared native search reference matches this dataset/context/asset")
        if len(choices) != 1:
            raise ValueError("Select an exact reference context and source asset")
        return choices[0]

    def inspect_reference(self, name, version=None, *, reference=None, asset_id=None, verify_files=True):
        """Inspect the separate reference cache offline; never repair or download."""
        from .search_reference import reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        return reference_registry(spec, self._cache_root).inspect(spec.dataset_id, spec.dataset_version,
                                                                 verify_files=verify_files)

    def download_reference(self, name, version=None, *, reference=None, asset_id=None,
                           force=False, source_path=None, **options):
        """Explicitly acquire a complete pinned search FASTA without rebuilding abundance."""
        from .search_reference import reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        overrides = {spec.source.asset_id: source_path} if source_path is not None else None
        return reference_registry(spec, self._cache_root).download(
            spec.dataset_id, spec.dataset_version, force=force, source_paths=overrides, **options)

    def open_reference(self, name, version=None, *, reference=None, asset_id=None):
        """Open and validate an installed native search set, without network access."""
        from .search_reference import SearchReference, reference_registry

        spec = self._search_reference(name, version, reference, asset_id)
        asset, version = spec.source.asset_id, spec.dataset_version
        hint = (f"The {name} {version} search reference {asset} is not installed; download it with "
                + install_hint(self, f"download_reference({name!r}, {version!r}, asset_id={asset!r})",
                               f"{name}:{asset}", "reference", version))
        path = _installed_path(reference_registry(spec, self._cache_root), spec.dataset_id, version, hint)
        return SearchReference(path, spec)

    def published(self, name, version=None):
        """Whether the release has a published artifact to download; otherwise it is built locally."""
        version = self._version(name, version)
        return bool(self._definitions[name]["versions"][version]["assets"])

    def _build_hint(self, name, version):
        release = self._definitions[name]["versions"][version]
        reason = "" if release["distribution"]["redistributable"] else ", because its source terms don't allow redistribution"
        return (f"The {name} {version} dataset has no published artifact{reason}; build it from its pinned "
                f"sources with " + install_hint(self, f"build({name!r}, {version!r})", name, "dataset", version))

    def _published(self, name, version):
        version = self._version(name, version)
        if not self.published(name, version):
            raise ValueError(self._build_hint(name, version))
        return version

    def _not_installed(self, name, version):
        return (f"The {name} {version} dataset is not installed; download it with "
                + install_hint(self, f"download({name!r}, {version!r})", name, "dataset", version))

    def _not_built(self, name, version):
        version = self._version(name, version)
        return NotInstalledError(
            f"The {name} {version} dataset has no local build; build it from its pinned sources with "
            + install_hint(self, f"build({name!r}, {version!r})", name, "dataset", version,
                           build=self.published(name, version)))

    def _data_registry(self):
        if self._registry is None:
            definitions = {}
            for name, entry in self._definitions.items():
                published = {v: r["assets"] for v, r in entry["versions"].items() if r["assets"] and "archive" not in r}
                if published:  # Every call names its version; the default only has to be a published one.
                    default = entry["default_version"]
                    definitions[name] = dict(default_version=default if default in published else min(published),
                                             versions=published)

            def store_path(name, version):
                release = self._definitions[name]["versions"][version]
                return get_data_dir(self._cache_root) / "datasets" / name / version / _digest(release)

            self._registry = _datacache().VersionedDatasetRegistry(definitions, store_path=store_path)
        return self._registry

    def _archive_store(self, name, version):
        release = self._definitions[name]["versions"][version]
        return get_data_dir(self._cache_root) / "datasets" / name / version / _digest(release)

    def inspect(self, name, version=None, *, verify_files=True):
        """Inspect local bytes offline; never download or repair a bad cache.

        A release that is only built locally reports its build.
        """
        if not self.published(name, version):
            return self.inspect_build(name, version, verify_files=verify_files)
        version = self._published(name, version)
        release = self._definitions[name]["versions"][version]
        if "archive" in release:
            return artifact_archives.inspect_archive(self._archive_store(name, version), release["archive"],
                                                      "dataset.sqlite", release["assets"]["dataset.sqlite"],
                                                      verify_files=verify_files)
        return self._data_registry().inspect(name, self._published(name, version), verify_files=verify_files)

    def download(self, name, version=None, *, force=False, source_paths=None, **options):
        """Explicitly install verified files; forward timeout/retry/progress options.

        source_paths accepts user-downloaded files while retaining declared URLs
        and hashes. force installs a fresh bundle without deleting readers' files.
        For an archive, source_paths["dataset.sqlite"] is the ordered part paths.
        """
        version = self._published(name, version)
        release = self._definitions[name]["versions"][version]
        if "archive" in release:
            if source_paths is not None and set(source_paths) != {"dataset.sqlite"}:
                raise ValueError("Archive overrides require dataset.sqlite's ordered part paths")

            def validate(path):
                if QuantitativeDataset(path).manifest.to_dict() != self.info(name, version).to_dict():
                    raise ValueError("Installed dataset manifest disagrees with registered release")

            path = artifact_archives.install_archive(
                self._archive_store(name, version), release["archive"], "dataset.sqlite", release["assets"]["dataset.sqlite"],
                source_paths=None if source_paths is None else source_paths["dataset.sqlite"], force=force,
                download_options=options, validator=validate)
            return {"dataset.sqlite": path}
        return self._data_registry().download(name, self._published(name, version), force=force,
                                             source_paths=source_paths, **options)

    def build(self, name, version=None, *, force=False, source_paths=None, download_options=None):
        """Rebuild a registered release from pinned source files, explicitly."""
        from .sources import importer

        manifest = self.info(name, version)
        transform = importer(manifest.processing.get("importer"))
        return build_dataset(manifest, lambda inputs, output: transform(inputs, output, manifest),
                             cache_root=self._cache_root, source_paths=source_paths, force=force,
                             download_options=download_options)

    def open(self, name, version=None, *, verify_files=True):
        """Open an installed dataset after integrity and identity checks, offline.

        A release that is only built locally opens its build.
        """
        if not self.published(name, version):
            if self.inspect_build(name, version, verify_files=False).status == "missing":
                version = self._version(name, version)
                raise NotInstalledError(self._build_hint(name, version).replace(
                    "has no published artifact", "hasn't been built, and has no published artifact"))
            return self.open_build(name, version, verify_files=verify_files)
        version = self._published(name, version)
        if "archive" in self._definitions[name]["versions"][version]:
            inspected = self.inspect(name, version, verify_files=verify_files)
            if inspected.status == "missing":
                raise NotInstalledError(self._not_installed(name, version))
            if inspected.status != "available":
                raise _datacache().FileValidationError(inspected.path, inspected.status) from inspected.error
            path = inspected.files["dataset.sqlite"].path
        else:
            path = _installed_path(self._data_registry(), name, version, self._not_installed(name, version),
                                   asset="dataset.sqlite", verify_files=verify_files)
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
        manifest = self.info(name, version)
        return _inspect_build(manifest, self._cache_root, verify_files=verify_files)

    def open_build(self, name, version=None, *, verify_files=True):
        """Open a local source materialization explicitly, without artifact fallback."""
        inspected = self.inspect_build(name, version, verify_files=verify_files)
        if inspected.status == "missing":
            raise self._not_built(name, version)
        if inspected.status != "available":
            raise _datacache().FileValidationError(inspected.path, inspected.status) from inspected.error
        dataset = QuantitativeDataset(inspected.files["dataset.sqlite"].path)
        if dataset.manifest.to_dict() != self.info(name, version).to_dict():
            raise ValueError("Built dataset manifest disagrees with the registered release")
        return dataset


def _inspect_build(manifest, cache_root, *, verify_files=True):
    inspected = _datacache().inspect_materialization(_build_store(manifest, cache_root),
                                                    verify_files=verify_files)
    if inspected.status != "available":
        return inspected
    try:
        expected = {source.asset_id: source for source in manifest.sources}
        if set(inspected.sources) != set(expected) or set(inspected.files) != {"dataset.sqlite"}:
            raise ValueError("Built source/output inventory disagrees with the registered release")
        for name, source in expected.items():
            recorded = inspected.sources[name]
            if (recorded["sha256"] != source.sha256 or recorded["size"] != source.size
                    or recorded["verified"] is not True):
                raise ValueError("Built source integrity pins disagree with the registered release")
        if inspected.transform != dict(version=manifest.schema, options=manifest.to_dict()):
            raise ValueError("Built transform/manifest disagrees with the registered release")
    except (KeyError, TypeError, ValueError) as exc:
        return replace(inspected, status="invalid", verified=False, error=exc)
    return inspected


def build_dataset(manifest, builder, *, cache_root=None, source_paths=None, force=False,
                  download_options=None):
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
            artifact_archives.archive_options(asset.metadata["archive"], asset.asset_id,
                                              dict(sha256=asset.sha256, size=asset.size))
    if not force and any("archive" in asset.metadata for asset in manifest.sources):
        cached = _inspect_build(manifest, cache_root)
        if cached.status == "available":
            dataset = QuantitativeDataset(cached.files["dataset.sqlite"].path)
            if dataset.manifest.to_dict() != manifest.to_dict():
                raise ValueError("Built dataset manifest disagrees with registered release")
            return dataset
        if cached.status != "missing":
            raise _datacache().FileValidationError(cached.path, cached.status) from cached.error
    sources = {}
    for asset in manifest.sources:
        source = dict(sha256=asset.sha256, size=asset.size)
        supplied = source_paths.get(asset.asset_id)
        if "archive" in asset.metadata and (supplied is None or isinstance(supplied, (list, tuple))):
            # Explicit builds may acquire a distributed source component. Queries
            # and open/inspect paths never call this acquisition branch.
            path = artifact_archives.install_archive(
                get_data_dir(cache_root) / "components" / _digest(asset.to_dict()),
                asset.metadata["archive"], asset.asset_id, source, source_paths=supplied,
                download_options=download_options)
            source["path"] = str(path)
        elif supplied is not None:
            source["path"] = str(supplied)
        else:
            source["url"] = asset.url
        sources[asset.asset_id] = source

    def build(inputs, outputs):
        builder(inputs, outputs["dataset.sqlite"])
        dataset = QuantitativeDataset(outputs["dataset.sqlite"])
        if dataset.manifest.to_dict() != manifest.to_dict():
            raise ValueError("Importer produced a different dataset manifest")

    identity = manifest.to_dict()
    destination = _build_store(manifest, cache_root)
    paths = _datacache().materialize(
        destination, sources, transform=dict(version=manifest.schema, options=identity),
        outputs={"dataset.sqlite": {}}, builder=build, force=force,
        download_options=download_options)
    return QuantitativeDataset(paths["dataset.sqlite"])
