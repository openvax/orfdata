"""Where ORFdata keeps what it downloads and builds, and the DataCache conventions it follows there.

Everything ORFdata installs lives under :func:`get_data_dir`: the ``orfdata/``
directory of ``ORFDATA_CACHE``, else ``OPENVAX_DATA_CACHE``, else the
platform's OpenVax cache. Each kind of resource has its own area there, such as
``catalogs/``, ``datasets/``, ``chains/``, ``genomes/`` and ``references/``, and
a store's directory is named by the :func:`content_digest` of the pins it holds,
so a changed pin never reuses an old store.

DataCache downloads, verifies, publishes and inspects the files. This module is
the one place ORFdata's own conventions around it live: the cache directory, how
stores are named, a registry for one pinned release, and the default download
timeout.
"""

import hashlib
import json
from pathlib import Path

import datacache

from .errors import NotInstalledError

# Seconds any ORFdata download may stall before it is retried, then fails; the CLI's --timeout default too.
DOWNLOAD_TIMEOUT = 60


def get_data_dir(override=None):
    """Return the directory ORFdata owns in the cache, without creating it.

    Parameters
    ----------
    override : str or os.PathLike, optional
        A cache root to use instead of ``ORFDATA_CACHE``, ``OPENVAX_DATA_CACHE``
        or the platform's OpenVax cache; what ``cache_root=`` and
        ``--cache-root`` pass.

    Returns
    -------
    pathlib.Path
        The ``orfdata`` directory of the selected cache root.
    """
    return Path(datacache.get_cache_root(
        "openvax", "ORFDATA_CACHE", "OPENVAX_DATA_CACHE", override=override)) / "orfdata"


def content_digest(value):
    """Return the SHA-256 of ``value`` written as canonical JSON, the name of a store holding it.

    Parameters
    ----------
    value : object
        JSON-serializable pins or manifest; the order of mapping keys doesn't
        change the digest.

    Returns
    -------
    str
        64 lowercase hexadecimal digits.

    Raises
    ------
    TypeError
        If ``value`` can't be written as JSON.
    """
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def pinned_registry(name, version, assets, store):
    """Return a DataCache registry of one pinned release, installed in ``store``.

    Nothing is read, downloaded or created until one of the registry's methods
    is called.

    Parameters
    ----------
    name : str
        The name the registry knows the release by.
    version : str
        The release's only version, which is also its default.
    assets : dict
        DataCache asset definitions keyed by file name, each with ``url``,
        ``sha256`` and ``size``, and optionally ``decompress``.
    store : str or os.PathLike
        The directory DataCache installs the release in.

    Returns
    -------
    datacache.VersionedDatasetRegistry
        Its ``download``, ``inspect`` and ``local_path`` take ``name`` and
        ``version``.
    """
    store = Path(store)
    return datacache.VersionedDatasetRegistry(
        {name: dict(default_version=version, versions={version: assets})},
        store_path=lambda name, version: store)


def fetch_pinned_file(url, sha256, size, directory, *, download_options=None):
    """Download one file pinned by SHA-256 and size into ``directory``, and return its path.

    The file is installed through DataCache in ``directory/<sha256>/``, so files
    with different pins never share a store. An installed file that still
    matches its pin is reused.

    Parameters
    ----------
    url : str
        Where to download the file; its last path segment is the file name.
    sha256 : str
        The file's expected SHA-256, as lowercase hexadecimal.
    size : int
        The file's expected size in bytes.
    directory : str or os.PathLike
        The directory that holds one store per pinned file.
    download_options : dict, optional
        DataCache download options, such as ``timeout`` or ``resume``.

    Returns
    -------
    str
        The installed file's path.

    Raises
    ------
    datacache.FileValidationError
        If the downloaded bytes don't match the pin.
    """
    name = url.rsplit("/", 1)[-1]
    registry = pinned_registry("source", "1", {name: dict(url=url, sha256=sha256, size=size)},
                               Path(directory) / sha256)
    registry.download("source", "1", **(download_options or {}))
    return str(registry.local_path("source", "1", asset=name))


def installed_path(registry, name, version, hint, **options):
    """Return an installed release's path from a DataCache registry, without downloading anything.

    Parameters
    ----------
    registry : datacache.VersionedDatasetRegistry
        The registry the release was installed through.
    name : str
        The release's name in ``registry``.
    version : str
        The release's version.
    hint : str or callable
        The message of the :class:`~orfdata.NotInstalledError` raised when the
        release isn't installed, or a function returning it, called only then.
    **options
        Passed on to ``registry.local_path``, such as ``asset`` and
        ``verify_files``.

    Returns
    -------
    pathlib.Path
        The ``asset`` file when one is named or the release has only one,
        otherwise the release's directory.

    Raises
    ------
    NotInstalledError
        If nothing is installed for the release.
    datacache.FileValidationError
        If installed files fail verification; DataCache's own error passes
        through.
    """
    try:
        return registry.local_path(name, version, **options)
    except FileNotFoundError as error:
        raise NotInstalledError(hint() if callable(hint) else hint) from error
