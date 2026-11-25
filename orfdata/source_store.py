"""Content-addressed SHA-256 cache for verified build inputs and resumable source transfers.

Build inputs shared across catalogues, datasets and ``orfdata setup`` runs —
such as Ensembl GTF/cDNA files or multi-species supplementary tables — are
stored once under ``<orfdata-cache>/sources/<sha256>/`` through DataCache.

Every lookup re-verifies the file's SHA-256 digest (and byte size when pinned)
before returning a path; bytes are never trusted by filename alone. When a
pinned source is not yet in ORFdata's store, :class:`SourceStore` can
optionally seed a verified hit from an existing local PyEnsembl cache without
touching the network.

Cached build sources can be inspected and removed explicitly with
:func:`prune_source_store` or ``orfdata data prune --sources --yes``.
"""

from __future__ import annotations

import os
import re
import shutil
import stat
from os import PathLike
from pathlib import Path
from typing import Mapping

from .cache import file_sha256, get_data_dir, pinned_registry

SOURCE_STORE_AREA = "sources"
_SHA256_RE = re.compile(r"[0-9a-f]{64}")


def _validate_sha256(sha256: str) -> str:
    if not isinstance(sha256, str) or _SHA256_RE.fullmatch(sha256) is None:
        raise ValueError(f"Source store keys must be lowercase 64-character SHA-256 digests, got {sha256!r}")
    return sha256


def _validate_size(size: int | None) -> int | None:
    if size is not None and (type(size) is not int or size < 0):
        raise ValueError(f"Expected source size must be a nonnegative integer or None, got {size!r}")
    return size


def source_store_dir(cache_root: str | PathLike[str] | None = None) -> Path:
    """Return the content-addressed source cache directory under ORFdata's cache root.

    Does not create the directory.

    Parameters
    ----------
    cache_root : str or os.PathLike, optional
        Explicit cache root override instead of ``ORFDATA_CACHE``,
        ``OPENVAX_DATA_CACHE`` or the platform default.

    Returns
    -------
    pathlib.Path
        ``<orfdata-cache>/sources``.
    """
    return get_data_dir(cache_root) / SOURCE_STORE_AREA


def pyensembl_cache_dirs() -> tuple[Path, ...]:
    """Return existing local PyEnsembl cache directories, without creating or downloading anything.

    Checks ``PYENSEMBL_CACHE_DIR`` and DataCache's resolved ``pyensembl`` cache
    directory. Only directories that already exist on disk are returned.

    Returns
    -------
    tuple of pathlib.Path
        Deduplicated existing PyEnsembl cache directories.
    """
    import datacache

    candidates: list[Path] = []
    env_dir = os.environ.get("PYENSEMBL_CACHE_DIR")
    if env_dir:
        candidates.append(Path(os.path.expanduser(env_dir)))
    try:
        resolved = datacache.resolve_cache_dir(
            "pyensembl",
            "PYENSEMBL_CACHE_DIR",
            append_package_to_env=False,
        )
        candidates.append(Path(resolved.path))
    except Exception:
        pass
    seen: set[Path] = set()
    existing: list[Path] = []
    for path in candidates:
        try:
            key = path.resolve()
        except OSError:
            continue
        if key not in seen and path.is_dir():
            seen.add(key)
            existing.append(path)
    return tuple(existing)


def _verify_candidate(path: Path, sha256: str, size: int | None) -> bool:
    """Return True iff ``path`` is a regular file matching ``size`` (when given) and ``sha256``."""
    try:
        info = path.stat()
    except OSError:
        return False
    if not stat.S_ISREG(info.st_mode):
        return False
    if size is not None and info.st_size != size:
        return False
    try:
        return file_sha256(path) == sha256
    except OSError:
        return False


def _dir_usage(path: Path) -> tuple[int, int]:
    """Return ``(total_bytes, file_count)`` for files under ``path`` without following symlinks."""
    total_bytes = 0
    file_count = 0
    for root, _dirs, filenames in os.walk(path, followlinks=False):
        for name in filenames:
            try:
                info = os.lstat(os.path.join(root, name))
            except OSError:
                continue
            total_bytes += info.st_size
            file_count += 1
    return total_bytes, file_count


class SourceStore:
    """Content-addressed SHA-256 store for checksum-pinned source files.

    Each pinned source lives in ``<orfdata-cache>/sources/<sha256>/`` managed
    through DataCache. Hits are always re-verified by SHA-256 before use, and
    downloads honour DataCache options such as ``timeout`` and ``resume``.

    Parameters
    ----------
    cache_root : str or os.PathLike, optional
        Explicit ORFdata cache root override.
    seed_from_pyensembl : bool, optional
        When ``True`` (the default), search existing local PyEnsembl cache
        directories for a byte-verified match before downloading.
    """

    def __init__(
        self,
        cache_root: str | PathLike[str] | None = None,
        *,
        seed_from_pyensembl: bool = True,
    ):
        self._cache_root = cache_root
        self._seed_from_pyensembl = bool(seed_from_pyensembl)

    @property
    def root(self) -> Path:
        """Directory holding per-SHA-256 source stores."""
        return source_store_dir(self._cache_root)

    def _asset_registry(self, url: str, sha256: str, size: int | None, filename: str):
        asset: dict[str, object] = {"url": url, "sha256": sha256}
        if size is not None:
            asset["size"] = size
        return pinned_registry("source", "1", {filename: asset}, self.root / sha256)

    def lookup(
        self,
        sha256: str,
        *,
        size: int | None = None,
        filename: str | None = None,
    ) -> Path | None:
        """Return a re-verified local path for ``sha256`` if present, else ``None``.

        Never trusts bytes by filename or cached receipt alone: every candidate
        file is hashed with SHA-256 before being returned.

        Parameters
        ----------
        sha256 : str
            Expected lowercase 64-character SHA-256 hex digest.
        size : int, optional
            Expected file size in bytes.
        filename : str, optional
            Preferred filename to check first inside the content-addressed store.

        Returns
        -------
        pathlib.Path or None
            Verified path on disk, or ``None`` if absent or damaged.
        """
        sha256 = _validate_sha256(sha256)
        size = _validate_size(size)
        store_dir = self.root / sha256
        if not store_dir.is_dir():
            return None

        # Prefer DataCache's published current bundle pointer when present.
        candidates: list[Path] = []
        current_link = store_dir / "current"
        if current_link.is_dir():
            if filename:
                preferred = current_link / filename
                if preferred.is_file():
                    candidates.append(preferred)
            for child in sorted(current_link.iterdir(), key=lambda p: p.name):
                if child.name.startswith(".datacache-") or child in candidates:
                    continue
                if child.is_file():
                    candidates.append(child)

        # Also inspect bundle generations if no current symlink exists.
        bundles_dir = store_dir / "bundles"
        if not candidates and bundles_dir.is_dir():
            generations = sorted(
                (p for p in bundles_dir.iterdir() if p.is_dir()),
                key=lambda p: p.name,
                reverse=True,
            )
            for gen in generations:
                for child in sorted(gen.iterdir(), key=lambda p: p.name):
                    if child.name.startswith(".datacache-") or not child.is_file():
                        continue
                    candidates.append(child)

        for candidate in candidates:
            if _verify_candidate(candidate, sha256, size):
                return candidate
        return None

    def find_pyensembl_seed(
        self,
        sha256: str,
        *,
        size: int | None = None,
        filename: str | None = None,
    ) -> Path | None:
        """Search existing PyEnsembl cache directories for a file matching ``sha256``.

        Only regular files matching ``size`` (when provided) are hashed, and no
        file is ever accepted on filename alone.

        Parameters
        ----------
        sha256 : str
            Expected lowercase 64-character SHA-256 hex digest.
        size : int, optional
            Expected byte size used to skip non-matching files quickly.
        filename : str, optional
            Filename to inspect before other size-matching files.

        Returns
        -------
        pathlib.Path or None
            Verified PyEnsembl cache file, or ``None`` if no match exists.
        """
        sha256 = _validate_sha256(sha256)
        size = _validate_size(size)
        preferred: list[Path] = []
        others: list[Path] = []
        for directory in pyensembl_cache_dirs():
            for walk_root, _dirs, filenames in os.walk(directory, followlinks=False):
                for name in sorted(filenames):
                    if name.startswith(".") or name.endswith(".db"):
                        continue
                    candidate = Path(walk_root) / name
                    try:
                        info = candidate.stat()
                    except OSError:
                        continue
                    if not stat.S_ISREG(info.st_mode):
                        continue
                    if size is not None and info.st_size != size:
                        continue
                    if filename is not None and name == filename:
                        preferred.append(candidate)
                    elif size is not None:
                        others.append(candidate)
        for candidate in (*preferred, *others):
            if _verify_candidate(candidate, sha256, size):
                return candidate
        return None

    def acquire(
        self,
        url: str,
        sha256: str,
        size: int | None = None,
        *,
        filename: str | None = None,
        download_options: Mapping[str, object] | None = None,
        force: bool = False,
    ) -> Path:
        """Return a verified local path for a pinned source, caching it by ``sha256``.

        Re-verifies any existing store hit before reusing it. On a miss, checks
        local PyEnsembl caches (when ``seed_from_pyensembl=True``) and otherwise
        downloads through DataCache with ``download_options`` (including
        ``timeout`` and ``resume``).

        Parameters
        ----------
        url : str
            Remote source URL.
        sha256 : str
            Expected lowercase 64-character SHA-256 hex digest.
        size : int, optional
            Expected byte size.
        filename : str, optional
            Asset filename inside the store (defaults to the URL's basename).
        download_options : mapping, optional
            DataCache download options such as ``timeout``, ``show_progress``,
            or ``resume``.
        force : bool, optional
            When ``True``, bypass an existing valid hit and re-acquire.

        Returns
        -------
        pathlib.Path
            Verified local path inside ``<orfdata-cache>/sources/<sha256>/``.
        """
        sha256 = _validate_sha256(sha256)
        size = _validate_size(size)
        name = filename or url.rsplit("/", 1)[-1] or "source"
        if not force:
            hit = self.lookup(sha256, size=size, filename=name)
            if hit is not None:
                return hit

        registry = self._asset_registry(url, sha256, size, name)
        store_dir = self.root / sha256
        # Replace if force was requested or a damaged/stale store already occupies this digest directory.
        replace_existing = force or store_dir.exists()
        opts = dict(download_options or {})

        if self._seed_from_pyensembl and not force:
            seed = self.find_pyensembl_seed(sha256, size=size, filename=name)
            if seed is not None:
                local_opts = {k: v for k, v in opts.items() if k != "resume"}
                registry.download(
                    "source",
                    "1",
                    force=replace_existing,
                    source_paths={name: str(seed)},
                    **local_opts,
                )
                return Path(registry.local_path("source", "1", asset=name, verify_files=True))

        registry.download("source", "1", force=replace_existing, **opts)
        return Path(registry.local_path("source", "1", asset=name, verify_files=True))


def prune_source_store(
    cache_root: str | PathLike[str] | None = None,
    *,
    apply: bool = False,
) -> dict[str, object]:
    """Inspect or remove cached source files under ``<orfdata-cache>/sources/``.

    Parameters
    ----------
    cache_root : str or os.PathLike, optional
        Explicit ORFdata cache root override.
    apply : bool, optional
        When ``False`` (the default), report what would be removed without
        deleting anything. When ``True``, delete cached source directories.

    Returns
    -------
    dict
        Summary with ``path``, ``dry_run``, ``removed_count``,
        ``reclaimable_bytes``, and ``entries``.
    """
    root = source_store_dir(cache_root)
    entries: list[dict[str, object]] = []
    reclaimable_bytes = 0
    if root.is_dir():
        for child in sorted(root.iterdir(), key=lambda p: p.name):
            if not child.is_dir():
                continue
            size_bytes, file_count = _dir_usage(child)
            reclaimable_bytes += size_bytes
            entries.append(
                {
                    "sha256": child.name,
                    "path": str(child),
                    "size_bytes": size_bytes,
                    "file_count": file_count,
                }
            )
            if apply:
                shutil.rmtree(child)
    return {
        "path": str(root),
        "dry_run": not apply,
        "removed_count": len(entries),
        "reclaimable_bytes": reclaimable_bytes,
        "entries": entries,
    }
