"""Verbatim source documents a viral collection is built from, pinned so they can be installed.

A viral collection downloads a GenBank or PaVE document per genome, in about 35
requests. These are those documents as fetched on a known day, in one pinned
archive, so a collection can be rebuilt, or checked against its published counts,
without the network and without upstream revisions changing the result. As with
chain files, :func:`install_records` (or ``orfdata data fetch``) downloads and
checks the archive once; :func:`read_records` then reads it offline.
"""

import json
import tarfile
from dataclasses import dataclass
from importlib.resources import files

from .cache import content_digest, get_data_dir, pinned_registry
from .errors import NotInstalledError


@dataclass(frozen=True)
class PinnedRecords:
    """An archive of source documents pinned by URL, size and SHA-256."""

    name: str
    species: str
    title: str
    url: str
    sha256: str
    size: int
    fetched: str
    orfdata_version: str
    provider: str

    @property
    def filename(self):
        return self.url.rsplit("/", 1)[-1]


def pinned_records():
    """The pinned archives, without downloading any."""
    pins = json.loads(files("orfdata").joinpath("data/source_records.json").read_text())
    return tuple(PinnedRecords(name, **entry) for name, entry in sorted(pins.items()))


def _pinned(name):
    found = {pin.name: pin for pin in pinned_records()}
    if name not in found:
        raise ValueError(f"No pinned source records {name!r}; pinned: {', '.join(found) or 'none'}")
    return found[name]


def records_registry(pin, cache_root=None):
    """Return the DataCache registry that installs a pinned archive.

    The release is named ``pin.name``, version ``"1"``, and its store is
    ``source-records/<name>/<digest of the pin>`` in the cache, so a new pin never
    reuses an old file. Nothing is downloaded until the registry is asked to.
    """
    assets = {pin.filename: dict(url=pin.url, sha256=pin.sha256, size=pin.size)}
    return pinned_registry(pin.name, "1", assets,
                           get_data_dir(cache_root) / "source-records" / pin.name
                           / content_digest(assets))


def install_records(name, *, cache_root=None, force=False, source_path=None, download_options=None):
    """Download the pinned archive and return its path; its SHA-256 and size must match the pin.

    ``source_path`` installs a copy you already have instead, checked the same way.
    """
    pin = _pinned(name)
    registry = records_registry(pin, cache_root)
    registry.download(pin.name, "1", force=force,
                      source_paths=None if source_path is None else {pin.filename: str(source_path)},
                      **(download_options or {}))
    return registry.local_path(pin.name, "1", asset=pin.filename)


def read_records(name, *, cache_root=None):
    """The installed archive's documents as ``{file name: bytes}``, still gzipped as published.

    File names read ``<abbreviation>__<source>__<identifier>.<suffix>.gz``. Nothing
    is downloaded: an archive that is not installed raises NotInstalledError.
    """
    pin = _pinned(name)
    registry = records_registry(pin, cache_root)
    if registry.inspect(pin.name, "1", verify_files=True).status != "available":
        raise NotInstalledError(
            f"{name!r} is not installed or fails its checksum; install it with: "
            f"orfdata data fetch {name} --kind records")
    path = registry.local_path(pin.name, "1", asset=pin.filename)
    with tarfile.open(path, "r:gz") as archive:
        return {member.name.removeprefix("./"): archive.extractfile(member).read()
                for member in archive if member.isfile()}
