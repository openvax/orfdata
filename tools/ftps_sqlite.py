"""Opt-in, bounded read-only access to a deposited SQLite database over FTPS.

This extraction tool is not imported by the normal ORFdata APIs. APSW delegates
all SQLite page, index and overflow decoding to SQLite itself. Acquired chunk
hashes establish integrity only for those bytes, not for the full remote file.
"""

import ftplib
import hashlib
import json
import resource
import ssl
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path


def sha256(data):
    return hashlib.sha256(data).hexdigest()


class SourceChanged(ValueError):
    """The remote asset no longer matches the pinned identity."""


class LimitExceeded(ValueError):
    """An explicit extraction resource limit would be exceeded."""


@dataclass(frozen=True)
class Identity:
    host: str
    path: str
    size: int
    mdtm: str


class FTPSRanges:
    """Anonymous explicit FTPS; each bounded stream is closed at its limit.

    Counters measure decrypted data payload returned by recv, including failed
    attempts and verification reads. They exclude TLS/control/network overhead
    and bytes a server sends ahead of the bounded client read.
    """

    def __init__(self, identity, *, transfer_limit, request_limit=10000,
                 seconds_limit=3600, retries=2, timeout=30):
        self.identity = identity
        self.transfer_limit = transfer_limit
        self.request_limit = request_limit
        self.seconds_limit = seconds_limit
        self.retries = retries
        self.timeout = timeout
        self.payload_bytes = 0
        self.requests = 0
        self.failures = []
        self.started = time.monotonic()

    def _connect(self):
        ftp = ftplib.FTP_TLS(context=ssl.create_default_context(), timeout=self.timeout)
        try:
            ftp.connect(self.identity.host, 21)
            ftp.login("anonymous", "anonymous@")
            ftp.prot_p()
            ftp.voidcmd("TYPE I")
        except BaseException:
            ftp.close()
            raise
        return ftp

    def check(self, ftp=None):
        own = ftp is None
        ftp = ftp or self._connect()
        try:
            found = Identity(self.identity.host, self.identity.path,
                             ftp.size(self.identity.path), ftp.sendcmd("MDTM " + self.identity.path))
            if found != self.identity:
                raise SourceChanged(f"Remote identity changed: {found!r}")
        finally:
            if own:
                ftp.close()

    def _check_budget(self, amount):
        if time.monotonic() - self.started > self.seconds_limit:
            raise LimitExceeded("Elapsed-time limit exceeded")
        if self.requests >= self.request_limit:
            raise LimitExceeded("Range request limit exceeded")
        if self.payload_bytes + amount > self.transfer_limit:
            raise LimitExceeded("Transfer payload limit exceeded")

    def read(self, offset, length):
        if offset < 0 or length <= 0 or offset + length > self.identity.size:
            raise ValueError("Range outside remote file")
        for attempt in range(self.retries + 1):
            self._check_budget(length)
            self.requests += 1
            ftp, stream = None, None
            try:
                ftp = self._connect()
                self.check(ftp)
                reply = ftp.sendcmd("REST " + str(offset))
                if not reply.startswith("350"):
                    raise ValueError("Server did not accept REST")
                stream = ftp.transfercmd("RETR " + self.identity.path)
                chunks, received = [], 0
                while received < length:
                    if time.monotonic() - self.started > self.seconds_limit:
                        raise LimitExceeded("Elapsed-time limit exceeded during range read")
                    chunk = stream.recv(min(65536, length - received))
                    if not chunk:
                        raise EOFError(f"Short range: {received} of {length} bytes")
                    self.payload_bytes += len(chunk)
                    received += len(chunk)
                    chunks.append(chunk)
                # Do not drain RETR: the remote file can be many gigabytes.
                stream.close()
                stream = None
                ftp.close()
                ftp = None
                self.check()
                return b"".join(chunks)
            except (OSError, EOFError, ftplib.error_temp) as error:
                self.failures.append(dict(offset=offset, length=length, attempt=attempt + 1,
                                          error=type(error).__name__))
                if attempt == self.retries:
                    raise
            finally:
                if stream is not None:
                    stream.close()
                if ftp is not None:
                    ftp.close()
        raise AssertionError("Unreachable retry state")


class ChunkCache:
    """Persistent checksummed chunks, with verified shifted overlap per chunk."""

    def __init__(self, source, directory, *, chunk_size=1024 * 1024,
                 cache_limit=64 * 1024 * 1024):
        if chunk_size < 4096 or cache_limit < chunk_size:
            raise ValueError("Cache must hold at least one SQLite-sized chunk")
        self.source = source
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.chunk_size = chunk_size
        self.cache_limit = cache_limit
        self.hits = 0
        self.read_calls = 0
        self.identity = asdict(source.identity)
        self.manifest_path = self.directory / "chunks.json"
        self.chunks = {}
        self.source.check()
        if self.manifest_path.exists():
            manifest = json.loads(self.manifest_path.read_text())
            if manifest["identity"] != self.identity or manifest["chunk_size"] != chunk_size:
                raise SourceChanged("Cache belongs to another source identity or chunk size")
            self.chunks = manifest["chunks"]
        expected = {"chunks.json", "receipt.json"} | {str(k) + ".bin" for k in self.chunks}
        if any(p.name not in expected for p in self.directory.iterdir()):
            raise ValueError("Untracked cache files; inspect the interrupted extraction before reuse")
        self.cached_bytes = sum(item["length"] for item in self.chunks.values())
        if self.cached_bytes > cache_limit:
            raise LimitExceeded("Existing cache exceeds cache limit")
        # Identity includes a stable content anchor in addition to SIZE/MDTM.
        anchor = source.read(0, min(4096, source.identity.size))
        self.header_sha256 = sha256(anchor)
        if self.manifest_path.exists() and manifest["header_sha256"] != self.header_sha256:
            raise SourceChanged("Remote header changed despite matching metadata")
        self._save()

    def _save(self):
        manifest = dict(identity=self.identity, header_sha256=self.header_sha256,
                        chunk_size=self.chunk_size, chunks=self.chunks)
        temporary = self.manifest_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        temporary.replace(self.manifest_path)

    def _chunk(self, offset):
        key = str(offset)
        path = self.directory / (key + ".bin")
        if key in self.chunks:
            data = path.read_bytes()
            entry = self.chunks[key]
            if len(data) != entry["length"] or sha256(data) != entry["sha256"]:
                raise ValueError("Cached chunk checksum/length mismatch")
            self.hits += 1
            return data
        length = min(self.chunk_size, self.source.identity.size - offset)
        # Reserve space before the network operation; do not evict and re-fetch
        # indefinitely under an apparently small disk budget.
        if self.cached_bytes + length > self.cache_limit:
            raise LimitExceeded("Cache limit exceeded")
        data = self.source.read(offset, length)
        shift = 17
        if length > shift:
            overlap = self.source.read(offset + shift, length - shift)
            if data[shift:] != overlap:
                raise ValueError("Shifted range overlap mismatch: REST offset not reliable")
        if offset == 0 and sha256(data[:min(4096, length)]) != self.header_sha256:
            raise SourceChanged("Remote header changed during chunk acquisition")
        path.write_bytes(data)
        self.chunks[key] = dict(length=length, sha256=sha256(data),
                                verification="independent REST +17 full overlap")
        self.cached_bytes += length
        self._save()
        return data

    def read(self, offset, length):
        self.read_calls += 1
        if offset < 0 or length <= 0 or offset + length > self.source.identity.size:
            raise ValueError("SQLite requested bytes outside the deposited file")
        output = []
        end = offset + length
        while offset < end:
            base = offset // self.chunk_size * self.chunk_size
            chunk = self._chunk(base)
            stop = min(end - base, len(chunk))
            output.append(chunk[offset - base:stop])
            offset = base + stop
        return b"".join(output)

    def receipt(self):
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return dict(identity=self.identity, header_sha256=self.header_sha256,
                    integrity_scope="metadata/header pin plus verified acquired chunks; no full-file checksum",
                    payload_bytes=self.source.payload_bytes, range_requests=self.source.requests,
                    failures=self.source.failures, cached_payload_bytes=self.cached_bytes,
                    disk_bytes=sum(p.stat().st_size for p in self.directory.iterdir() if p.is_file()),
                    peak_process_rss_bytes=peak if sys.platform == "darwin" else peak * 1024,
                    elapsed_seconds=time.monotonic() - self.source.started,
                    cache_hits=self.hits, sqlite_reads=self.read_calls,
                    transfer_limit=self.source.transfer_limit, cache_limit=self.cache_limit,
                    request_limit=self.source.request_limit, seconds_limit=self.source.seconds_limit,
                    accounting="recv data payload, including overlap/retries; excludes network/TLS overhead and server read-ahead")


def connect(cache, *, memory_limit=128 * 1024 * 1024):
    """Return a read-only APSW connection and its live VFS (keep both alive)."""
    import apsw

    class File:
        def xRead(self, amount, offset):
            return cache.read(offset, amount)

        def xFileSize(self):
            return cache.source.identity.size

        def xClose(self):
            pass

        def xLock(self, level):
            pass

        def xUnlock(self, level):
            pass

        def xCheckReservedLock(self):
            return False

        def xFileControl(self, op, ptr):
            return False

        def xDeviceCharacteristics(self):
            return apsw.SQLITE_IOCAP_IMMUTABLE

        def xSectorSize(self):
            return 4096

        def xWrite(self, data, offset):
            raise apsw.ReadOnlyError("Deposited file is read-only")

        def xTruncate(self, size):
            raise apsw.ReadOnlyError("Deposited file is read-only")

    class VFS(apsw.VFS):
        def xOpen(self, name, flags):
            if not flags[0] & apsw.SQLITE_OPEN_MAIN_DB or flags[0] & apsw.SQLITE_OPEN_READWRITE:
                raise apsw.ReadOnlyError("Only the read-only main database is available")
            flags[1] = apsw.SQLITE_OPEN_READONLY
            return File()

        def xFullPathname(self, name):
            return name

        def xAccess(self, name, flags):
            return False

        def xDelete(self, name, syncdir):
            raise apsw.ReadOnlyError("Deposited file is read-only")

    name = f"orfdata-ftps-{id(cache)}"
    vfs = VFS(name, "")
    database = apsw.Connection("deposited.sqlite", flags=apsw.SQLITE_OPEN_READONLY, vfs=name)
    database.pragma("query_only", True)
    database.pragma("trusted_schema", False)
    database.pragma("temp_store", "memory")
    database.pragma("cache_size", -8192)
    database.pragma("hard_heap_limit", memory_limit)
    database.limit(apsw.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
    return database, vfs
