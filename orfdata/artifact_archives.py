"""Pinned single-file tar archives, including ordered distribution parts.

DataCache owns acquisition, extraction, transactions and cache publication.
ORFdata additionally pins the extracted file independently of cache receipts.
"""

import re
import stat
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path


def _pin(value):
    if not isinstance(value, Mapping) or set(value) != {"sha256", "size"}:
        raise ValueError("An archive member requires exactly its SHA-256 and size")
    if not isinstance(value["sha256"], str) or not re.fullmatch("[0-9a-f]{64}", value["sha256"]):
        raise ValueError("An archive requires a lowercase SHA-256 digest")
    if type(value["size"]) is not int or value["size"] < 0:
        raise ValueError("An archive size must be a nonnegative integer")


def archive_options(archive, member, pin):
    """Validate declarative metadata without filesystem or network access."""
    if not isinstance(member, str) or not member or member in (".", "..", ".datacache-manifest.json") or any(
            c in member for c in "/\\:"):
        raise ValueError("An archive member must be a single regular filename")
    _pin(pin)
    if not isinstance(archive, Mapping) or set(archive) != {"format", "parts", "sha256", "size"}:
        raise ValueError("An archive requires format, ordered parts, SHA-256 and size")
    if archive["format"] != "tar.gz":
        raise ValueError("Unsupported dataset archive format")
    _pin({key: archive[key] for key in ("sha256", "size")})
    parts = archive["parts"]
    if not isinstance(parts, (list, tuple)) or not 1 <= len(parts) <= 1000:
        raise ValueError("An archive requires 1–1000 ordered parts")
    sources = []
    for part in parts:
        if not isinstance(part, Mapping) or set(part) != {"url", "sha256", "size"}:
            raise ValueError("Every archive part requires URL, SHA-256 and size")
        _pin({key: part[key] for key in ("sha256", "size")})
        if not isinstance(part["url"], str) or not part["url"].strip() or part["size"] == 0:
            raise ValueError("Archive parts require nonempty URLs and bytes")
        sources.append(dict(part))
    if len({p["url"] for p in sources}) != len(sources) or sum(p["size"] for p in sources) != archive["size"]:
        raise ValueError("Archive parts must have distinct URLs and sum to the assembled size")
    return dict(sources=sources, expected_sha256=archive["sha256"], expected_size=archive["size"])


def _member(tree, member, pin, *, verify_files):
    """Check one independently pinned file, never a hash supplied by the cache."""
    import datacache

    tree = Path(tree)
    if {p.name for p in tree.iterdir()} - {member, ".datacache-manifest.json"}:
        raise ValueError("Dataset archive contains unexpected files")
    path = tree / member
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError("Dataset archive member must be a regular file without links")
    inspected = datacache.inspect_file(path, expected_sha256=pin["sha256"] if verify_files else None,
                                      expected_size=pin["size"])
    if inspected.status != "available":
        raise inspected.error
    return inspected


def install_archive(destination, archive, member, pin, *, source_paths=None, force=False,
                    download_options=None, validator=None):
    """Install ordered parts atomically, validating the member before publication.

    Local overrides are an ordered sequence, one path per declared part; they
    retain the original URLs and checksums. No custom downloader or tar reader.
    """
    import datacache

    from .datasets import _resumable

    options = archive_options(archive, member, pin)
    if source_paths is not None:
        if not isinstance(source_paths, (list, tuple)) or len(source_paths) != len(options["sources"]):
            raise ValueError("Supply one local path for every ordered archive part")
        options["sources"] = [dict(source, path=Path(path)) for source, path in zip(options["sources"], source_paths)]

    def validate(tree):
        checked = _member(tree, member, pin, verify_files=True)
        if validator is not None:
            return validator(Path(checked.path))

    tree = datacache.install_archive(destination, **options, force=force, verified=True,
                                    max_members=1, max_extracted_size=pin["size"],
                                    download_options=_resumable(download_options), validator=validate)
    return Path(tree) / member


def inspect_archive(destination, archive, member, pin, *, verify_files=True):
    """Check archive identity and the independently pinned member entirely offline."""
    import datacache

    if type(verify_files) is not bool:
        raise ValueError("verify_files must be a boolean")
    options = archive_options(archive, member, pin)
    # Validate the managed archive identity, then hash its only member against
    # the independent registry pin. This avoids hashing a large file twice and
    # rejects an attacker changing both payload and its cached receipt hash.
    inspected = datacache.inspect_archive(destination, **options, verify_files=False)
    if inspected.status != "available":
        return inspected
    try:
        checked = _member(inspected.bundle, member, pin, verify_files=verify_files)
    except PermissionError as error:
        return replace(inspected, status="inaccessible", verified=False, error=error)
    except (OSError, ValueError) as error:
        return replace(inspected, status="invalid", verified=False, error=error)
    return replace(inspected, files={member: checked}, verified=checked.verified)
