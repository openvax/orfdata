"""Explicit optional catalogue downloads and builds using OpenVax DataCache."""

import gzip
import json
from dataclasses import replace
from importlib.resources import files

from .catalog import ORFCatalog, load_catalog
from .catalog_manifest import CatalogManifest
from .datasets import _datacache, _digest, get_data_dir
from .reference import species_taxon

_ASSET = "catalog.json.gz"


def _build_store(manifest, cache_root):
    return (get_data_dir(cache_root) / "catalog-imports" / manifest.catalog_id /
            manifest.version / _digest(manifest.to_dict()))


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


class CatalogRegistry:
    """Discover offline, explicitly acquire, then open the existing ORFCatalog.

    Each definition names a pinned default_version, manifests, distribution
    terms and one checksummed catalog.json.gz asset. No method except download
    or build may acquire data. Bundled load_catalog() remains independent.
    """

    def __init__(self, definitions=None, *, cache_root=None):
        if definitions is None:
            definitions = json.loads(files("orfdata").joinpath("data/catalogs.json").read_text())
        self._definitions = json.loads(json.dumps(definitions, allow_nan=False))
        self._cache_root = cache_root
        self._manifests = {}
        self._registry = None
        for name, entry in self._definitions.items():
            if entry["default_version"] not in entry["versions"]:
                raise ValueError("A catalogue default must be a pinned version")
            for version, release in entry["versions"].items():
                manifest = CatalogManifest.from_dict(release["manifest"])
                if (manifest.catalog_id, manifest.version) != (name, version):
                    raise ValueError("Registry and manifest catalogue identities disagree")
                if set(release["assets"]) != {_ASSET} or release["assets"][_ASSET].get("decompress", False):
                    raise ValueError("A catalogue release requires raw catalog.json.gz bytes")
                distribution = release["distribution"]
                if not isinstance(distribution.get("access_terms"), str) or not distribution["access_terms"].strip():
                    raise ValueError("Catalogue artifact access terms are required")
                if type(distribution.get("redistributable")) is not bool:
                    raise ValueError("Catalogue artifact redistribution must be explicit")
                self._manifests[name, version] = manifest

    def _version(self, name, version):
        if name not in self._definitions:
            raise ValueError(f"Unknown catalogue {name!r}; see 'orfdata catalogs list' or CatalogRegistry().list()")
        version = self._definitions[name]["default_version"] if version is None else version
        if (name, version) not in self._manifests:
            raise ValueError(f"Unknown catalogue version {name!r} {version!r}")
        return version

    def list(self, *, species=None, taxon_id=None):
        """Return matching versioned manifests without downloading data."""
        if species is not None:
            found = species_taxon(species)
            if found is None:
                raise ValueError("Unknown species alias; supply its explicit taxon_id")
            if taxon_id is not None and taxon_id != found:
                raise ValueError("Species and taxon_id disagree")
            taxon_id = found
        return tuple(m for _, m in sorted(self._manifests.items()) if taxon_id is None or m.taxon_id == taxon_id)

    def info(self, name, version=None):
        return self._manifests[name, self._version(name, version)]

    def download_metadata(self, name, version=None):
        release = self._definitions[name]["versions"][self._version(name, version)]
        return json.loads(json.dumps({k: release[k] for k in ("distribution", "assets")}))

    def _data_registry(self):
        if self._registry is None:
            definitions = {name: dict(default_version=entry["default_version"], versions={
                version: release["assets"] for version, release in entry["versions"].items()})
                for name, entry in self._definitions.items()}

            def store_path(name, version):
                release = self._definitions[name]["versions"][version]
                return get_data_dir(self._cache_root) / "catalogs" / name / version / _digest(release)

            self._registry = _datacache().VersionedDatasetRegistry(definitions, store_path=store_path)
        return self._registry

    def inspect(self, name, version=None, *, verify_files=True):
        return self._data_registry().inspect(name, self._version(name, version), verify_files=verify_files)

    def download(self, name, version=None, *, force=False, source_paths=None, **options):
        """Explicitly install a pinned artifact; accept verified local overrides."""
        return self._data_registry().download(name, self._version(name, version), force=force,
                                             source_paths=source_paths, **options)

    def open(self, name, version=None, *, verify_files=True):
        """Validate installed bytes and content, then return an offline ORFCatalog."""
        manifest = self.info(name, version)
        path = self._data_registry().local_path(name, manifest.version, asset=_ASSET, verify_files=verify_files)
        return _open(path, manifest)

    def build(self, name, version=None, *, force=False, source_paths=None, download_options=None):
        """Explicitly materialize an allowlisted catalogue source adapter."""
        manifest = self.info(name, version)
        if manifest.processing.get("importer") != "catalog-copy-v1":
            raise ValueError("No packaged importer for this catalogue release")
        return build_catalog(manifest, lambda inputs, output: _copy_catalog(inputs, output, manifest),
                             cache_root=self._cache_root, source_paths=source_paths, force=force,
                             download_options=download_options)

    def inspect_build(self, name, version=None, *, verify_files=True):
        """Inspect build receipts offline; their hashes are not released-artifact pins."""
        manifest = self.info(name, version)
        inspected = _datacache().inspect_materialization(_build_store(manifest, self._cache_root),
                                                        verify_files=verify_files)
        if inspected.status != "available":
            return inspected
        try:
            expected = {s.asset_id: s for s in manifest.sources}
            if set(inspected.sources) != set(expected) or set(inspected.files) != {_ASSET}:
                raise ValueError("Built catalogue source/output inventory disagrees")
            for name, source in expected.items():
                recorded = inspected.sources[name]
                if (recorded["sha256"] != source.sha256 or recorded["size"] != source.size
                        or recorded["verified"] is not True):
                    raise ValueError("Built catalogue source integrity pins disagree")
            if inspected.transform != dict(version=manifest.schema, options=manifest.to_dict()):
                raise ValueError("Built catalogue transform/manifest disagrees")
        except (KeyError, TypeError, ValueError) as exc:
            return replace(inspected, status="invalid", verified=False, error=exc)
        return inspected

    def open_build(self, name, version=None, *, verify_files=True):
        """Open the local source build explicitly; never fall back to a download."""
        inspected = self.inspect_build(name, version, verify_files=verify_files)
        if inspected.status == "missing":
            raise FileNotFoundError(inspected.path)
        if inspected.status != "available":
            raise _datacache().FileValidationError(inspected.path, inspected.status) from inspected.error
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

    paths = _datacache().materialize(
        _build_store(manifest, cache_root), sources,
        transform=dict(version=manifest.schema, options=manifest.to_dict()),
        outputs={_ASSET: {}}, builder=build, force=force, download_options=download_options)
    return _open(paths[_ASSET], manifest)
