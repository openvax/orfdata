"""Explicit, checksummed dataset downloads using OpenVax DataCache conventions."""

import hashlib
import json
from importlib.resources import files
from pathlib import Path

from .dataset import QuantitativeDataset
from .measurements import DatasetManifest
from .reference import species_taxon


def _datacache():
    try:
        import datacache
    except ImportError as exc:
        raise ImportError("Install orfdata[data] to manage quantitative dataset files") from exc
    return datacache


def get_data_dir(override=None):
    """Resolve ORFDATA_CACHE, OPENVAX_DATA_CACHE, or the platform OpenVax cache.

    ORFdata owns only the orfdata/ child; lookup does not create directories.
    """
    return Path(_datacache().get_cache_root(
        "openvax", "ORFDATA_CACHE", "OPENVAX_DATA_CACHE", override=override)) / "orfdata"


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class DatasetRegistry:
    """Discover manifests offline; explicitly download before opening data.

    Definitions map names to a pinned default_version and versions. Each version
    contains a manifest, distribution access terms, and checksummed assets. No
    dataset is downloaded on construction, discovery, inspection or open.
    """

    def __init__(self, definitions=None, *, cache_root=None):
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
                if set(release["assets"]) != {"dataset.sqlite"}:
                    raise ValueError("A quantitative release must contain dataset.sqlite")
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
        return json.loads(json.dumps({key: release[key] for key in ("distribution", "assets")}))

    def _data_registry(self):
        if self._registry is None:
            definitions = {name: dict(default_version=entry["default_version"], versions={
                version: release["assets"] for version, release in entry["versions"].items()})
                for name, entry in self._definitions.items()}

            def store_path(name, version):
                release = self._definitions[name]["versions"][version]
                return get_data_dir(self._cache_root) / "datasets" / name / version / _digest(release)

            self._registry = _datacache().VersionedDatasetRegistry(definitions, store_path=store_path)
        return self._registry

    def inspect(self, name, version=None, *, verify_files=True):
        """Inspect local bytes offline; never download or repair a bad cache."""
        version = self._version(name, version)
        return self._data_registry().inspect(name, version, verify_files=verify_files)

    def download(self, name, version=None, *, force=False, source_paths=None, **options):
        """Explicitly install verified files; forward timeout/retry/progress options.

        source_paths accepts user-downloaded files while retaining declared URLs
        and hashes. force installs a fresh bundle without deleting readers' files.
        """
        version = self._version(name, version)
        return self._data_registry().download(name, version, force=force,
                                             source_paths=source_paths, **options)

    def open(self, name, version=None, *, verify_files=True):
        """Open an installed dataset after integrity and identity checks, offline."""
        version = self._version(name, version)
        path = self._data_registry().local_path(name, version, asset="dataset.sqlite",
                                                verify_files=verify_files)
        dataset = QuantitativeDataset(path)
        if dataset.manifest.to_dict() != self.info(name, version).to_dict():
            raise ValueError("Installed dataset manifest disagrees with registered release")
        return dataset


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
    sources = {}
    for asset in manifest.sources:
        source = dict(sha256=asset.sha256, size=asset.size)
        if asset.asset_id in source_paths:
            source["path"] = str(source_paths[asset.asset_id])
        else:
            source["url"] = asset.url
        sources[asset.asset_id] = source

    def build(inputs, outputs):
        builder(inputs, outputs["dataset.sqlite"])
        dataset = QuantitativeDataset(outputs["dataset.sqlite"])
        if dataset.manifest.to_dict() != manifest.to_dict():
            raise ValueError("Importer produced a different dataset manifest")

    identity = manifest.to_dict()
    destination = (get_data_dir(cache_root) / "imports" / manifest.dataset_id /
                   manifest.version / _digest(identity))
    paths = _datacache().materialize(
        destination, sources, transform=dict(version=manifest.schema, options=identity),
        outputs={"dataset.sqlite": {}}, builder=build, force=force,
        download_options=download_options)
    return QuantitativeDataset(paths["dataset.sqlite"])
