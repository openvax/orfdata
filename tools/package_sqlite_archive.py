"""Create bounded, deterministic split tar.gz archives of pinned local files."""

import argparse
import gzip
import hashlib
import json
import os
import shutil
import tarfile
import time
from pathlib import Path

from orfdata.artifact_archives import archive_options
from orfdata.sources.pfizer_native import file_sha256

CODE_SHA256 = file_sha256(__file__)


class _PartsWriter:
    def __init__(self, directory, prefix, part_bytes, max_bytes):
        self.directory, self.prefix = directory, prefix
        self.part_bytes, self.max_bytes = part_bytes, max_bytes
        self.count, self.digest = 0, hashlib.sha256()
        self.parts, self.owned, self.handle = [], [], None

    def _close_part(self):
        if self.handle is not None:
            self.handle.flush()
            os.fsync(self.handle.fileno())
            self.handle.close()
            self.handle = None
            self.parts.append(dict(name=self.owned[-1].name, size=self.part_count, sha256=self.part_digest.hexdigest()))

    def write(self, block):
        if self.count + len(block) > self.max_bytes:
            raise ValueError("Compressed archive exceeds its explicit byte limit")
        view = memoryview(block)
        while view:
            if self.handle is None:
                path = self.directory / f"{self.prefix}.part{len(self.parts):04d}"
                self.handle = path.open("xb")
                self.owned.append(path)
                self.part_count, self.part_digest = 0, hashlib.sha256()
            chunk = view[:self.part_bytes - self.part_count]
            if self.handle.write(chunk) != len(chunk):
                raise OSError("Short archive write")
            self.count += len(chunk)
            self.part_count += len(chunk)
            self.digest.update(chunk)
            self.part_digest.update(chunk)
            view = view[len(chunk):]
            if self.part_count == self.part_bytes:
                self._close_part()
        return len(block)

    def close(self):
        self._close_part()

    def discard(self):
        if self.handle is not None:
            self.handle.close()
            self.handle = None
        for path in self.owned:
            path.unlink()


class _PartsReader:
    def __init__(self, paths):
        self.paths, self.handle = iter(paths), None

    def read(self, size):
        if size < 0:
            raise ValueError("Archive verification requires bounded reads")
        result = bytearray()
        while len(result) < size:
            if self.handle is None:
                path = next(self.paths, None)
                if path is None:
                    break
                self.handle = path.open("rb")
            block = self.handle.read(size - len(result))
            if block:
                result.extend(block)
            else:
                self.handle.close()
                self.handle = None
        return bytes(result)

    def close(self):
        if self.handle is not None:
            self.handle.close()
            self.handle = None


def archive_verified(source, output, *, member, url_prefix, part_bytes, max_bytes, expected_sha256=None):
    """Write a new split archive and verify every restored byte without a spool.

    This primitive establishes byte preservation, not biological acceptance.
    The consumer-artifact packager separately requires full-source verification.
    """
    if any(type(v) is not int or v < 1 for v in (part_bytes, max_bytes)):
        raise ValueError("Archive byte limits must be positive integers")
    if part_bytes >= 2**31 or max_bytes > 1000 * part_bytes:
        raise ValueError("Use parts smaller than 2 GiB and at most 1000 parts")
    if not isinstance(url_prefix, str) or not url_prefix.startswith("https://") or not url_prefix.endswith("/"):
        raise ValueError("Supply the intended HTTPS release URL prefix ending in /")
    if not isinstance(member, str) or not member or Path(member).name != member or any(c in member for c in "\\:"):
        raise ValueError("Archive member must be a single filename")
    if member in (".", "..", ".datacache-manifest.json"):
        raise ValueError("Reserved archive member name")
    source, output = Path(source), Path(output)
    if output.exists():
        raise FileExistsError(output)
    source_size, source_hash = source.stat().st_size, file_sha256(source)
    if expected_sha256 is not None and source_hash != expected_sha256:
        raise ValueError("Source differs from the expected SHA-256")
    if shutil.disk_usage(output.parent).free < max_bytes + 1024**2:
        raise ValueError("Insufficient free disk for the bounded compressed archive")
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    writer = _PartsWriter(output, member + ".tar.gz", part_bytes, max_bytes)
    try:
        with gzip.GzipFile(filename="", fileobj=writer, mode="wb", mtime=0, compresslevel=6) as compressed:
            # GNU numeric headers support files above 8 GiB without extra members.
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.GNU_FORMAT) as archive:
                info = tarfile.TarInfo(member)
                info.size, info.mode, info.mtime = source_size, 0o644, 0
                with source.open("rb") as handle:
                    archive.addfile(info, handle)
        writer.close()
        restored_hash, restored_size = hashlib.sha256(), 0
        reader = _PartsReader(writer.owned)
        try:
            with gzip.GzipFile(fileobj=reader, mode="rb") as restored:
                with tarfile.open(fileobj=restored, mode="r|") as archive:
                    info = archive.next()
                    if info is None or info.name != member or not info.isfile() or info.size != source_size:
                        raise ValueError("Archive member differs from its pinned source")
                    with archive.extractfile(info) as handle:
                        while block := handle.read(1024**2):
                            restored_size += len(block)
                            if restored_size > source_size:
                                raise ValueError("Archive expanded beyond its pinned source size")
                            restored_hash.update(block)
                    if archive.next() is not None:
                        raise ValueError("Archive contains an unexpected member")
                # Finish gzip CRC/trailer checking, with bounded tar padding.
                tail_size = 0
                while block := restored.read(1024**2):
                    tail_size += len(block)
                    if tail_size > tarfile.RECORDSIZE or any(block):
                        raise ValueError("Unexpected bytes after the archive member")
        finally:
            reader.close()
        if restored_size != source_size or restored_hash.hexdigest() != source_hash or file_sha256(source) != source_hash:
            raise ValueError("Archive does not reproduce the unchanged source bytes")
        for path, part in zip(writer.owned, writer.parts):
            if path.stat().st_size != part["size"] or file_sha256(path) != part["sha256"]:
                raise ValueError("Archive part changed during verification")
        archive = dict(format="tar.gz", sha256=writer.digest.hexdigest(), size=writer.count,
                       parts=[dict(url=url_prefix + p["name"], size=p["size"], sha256=p["sha256"]) for p in writer.parts])
        archive_options(archive, member, dict(size=source_size, sha256=source_hash))
        return dict(archive=archive, installed_size=source_size, installed_sha256=source_hash,
                    elapsed_seconds=time.monotonic() - started, part_bytes=part_bytes, max_bytes=max_bytes,
                    compression=dict(format="tar.gz", level=6, mtime=0, tar_format="GNU"),
                    tool_sha256=CODE_SHA256)
    except BaseException:
        writer.discard()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path, help="New directory for parts and receipts")
    parser.add_argument("--member", required=True, help="Extracted filename, such as evidence.sqlite")
    parser.add_argument("--url-prefix", required=True, help="Intended HTTPS release URL prefix ending in /")
    parser.add_argument("--part-mib", required=True, type=int)
    parser.add_argument("--max-compressed-mib", required=True, type=int)
    parser.add_argument("--expected-sha256", help="Require this independently recorded source checksum before writing")
    args = parser.parse_args()
    result = archive_verified(args.source, args.output, member=args.member, url_prefix=args.url_prefix,
                              part_bytes=args.part_mib * 1024**2, max_bytes=args.max_compressed_mib * 1024**2,
                              expected_sha256=args.expected_sha256)
    for name, value in (("archive.json", result["archive"]), ("archive-receipt.json", result)):
        with (args.output / name).open("x") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
