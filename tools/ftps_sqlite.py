"""Opt-in, bounded read-only access to a deposited SQLite database over FTPS.

This extraction tool is not imported by the normal ORFdata APIs. APSW delegates
all SQLite page, index and overflow decoding to SQLite itself. Acquired chunk
hashes establish integrity only for those bytes, not for the full remote file.
"""

import fcntl
import ftplib
import hashlib
import importlib.metadata
import json
import os
import re
import resource
import ssl
import sys
import time
from bisect import bisect_right, insort
from dataclasses import asdict, dataclass
from pathlib import Path

from tools.extraction_budget import CODE_SHA256 as BUDGET_CODE_SHA256
from tools.extraction_budget import LimitExceeded, ScratchBudget, TransferLedger

CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


class SourceChanged(ValueError):
    """The remote asset no longer matches the pinned identity."""


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
                 seconds_limit=3600, retries=2, timeout=30, reuse_control=False):
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
        self.ledger = None
        self.reuse_control = reuse_control
        self._control = None
        self.control_connections = 0
        self.bounded_close_responses = 0

    def _connect(self):
        self.control_connections += 1
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

    def _get_control(self):
        if self._control is None:
            self._control = self._connect()
        return self._control

    def _finish_bounded_transfer(self, ftp):
        """Consume the observed server response and prove command synchronization.

        MassIVE returns a final 550 network-write error after a deliberately
        bounded data socket is closed. Only that observed response is accepted,
        only after the requested payload length was received. A following 200
        NOOP and fresh SIZE/MDTM must succeed before the channel is reused.
        Chunk overlap checks independently establish returned-offset integrity.
        """
        try:
            ftp.voidresp()
        except ftplib.error_perm as error:
            message = str(error)
            if not (self.identity.host == "massive-ftp.ucsd.edu" and
                    message.startswith("550 Issue during transfer: network error:") and
                    "write:" in message and ("connection reset by peer" in message or "broken pipe" in message)):
                raise
            self.bounded_close_responses += 1
        except ftplib.error_temp as error:
            if not str(error).startswith("426"):
                raise
            ftp.voidresp()  # A 426 may be followed by its final 226 response.
        if not ftp.voidcmd("NOOP").startswith("200"):
            raise ValueError("FTPS control responses are not synchronized")

    def close(self):
        if self._control is not None:
            self._control.close()
            self._control = None

    def check(self, ftp=None):
        own = ftp is None
        ftp = ftp or (self._get_control() if self.reuse_control else self._connect())
        try:
            found = Identity(self.identity.host, self.identity.path,
                             ftp.size(self.identity.path), ftp.sendcmd("MDTM " + self.identity.path))
            if found != self.identity:
                raise SourceChanged(f"Remote identity changed: {found!r}")
        finally:
            if own and not self.reuse_control:
                ftp.close()

    def _check_budget(self, amount):
        if time.monotonic() - self.started > self.seconds_limit:
            raise LimitExceeded("Elapsed-time limit exceeded")
        usage = self.ledger.summary() if self.ledger else dict(
            range_requests=self.requests, charged_payload_bytes=self.payload_bytes)
        if usage["range_requests"] >= self.request_limit:
            raise LimitExceeded("Range request limit exceeded")
        if usage["charged_payload_bytes"] + amount > self.transfer_limit:
            raise LimitExceeded("Transfer payload limit exceeded")

    def read(self, offset, length):
        if offset < 0 or length <= 0 or offset + length > self.identity.size:
            raise ValueError("Range outside remote file")
        for attempt in range(self.retries + 1):
            self._check_budget(length)
            reservation = self.ledger.reserve(offset, length, byte_limit=self.transfer_limit,
                                               request_limit=self.request_limit) if self.ledger else None
            self.requests += 1
            ftp, stream = None, None
            initial_bytes = self.payload_bytes
            status = "failed"
            try:
                ftp = self._get_control() if self.reuse_control else self._connect()
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
                if self.reuse_control:
                    self._finish_bounded_transfer(ftp)
                    self.check(ftp)
                else:
                    ftp.close()
                    ftp = None
                    self.check()
                status = "complete"
                return b"".join(chunks)
            except (OSError, EOFError, ftplib.error_temp) as error:
                self.failures.append(dict(offset=offset, length=length, attempt=attempt + 1,
                                          error=type(error).__name__))
                if attempt == self.retries:
                    raise
            finally:
                try:
                    if stream is not None:
                        stream.close()
                    if ftp is not None and (not self.reuse_control or status != "complete"):
                        ftp.close()
                        if self._control is ftp:
                            self._control = None
                finally:
                    if reservation is not None:
                        self.ledger.settle(reservation, self.payload_bytes - initial_bytes, status=status)
        raise AssertionError("Unreachable retry state")


class ChunkCache:
    """Persistent checksummed chunks, with verified shifted overlap per chunk."""

    def __init__(self, source, directory, *, chunk_size=1024 * 1024,
                 cache_limit=64 * 1024 * 1024, prior_receipts=(), scratch=None, read_ahead=None):
        if not 4096 <= chunk_size <= 16 * 1024 * 1024 or chunk_size & (chunk_size - 1):
            raise ValueError("Chunk size must be a power of two from 4 KiB to 16 MiB")
        if cache_limit < chunk_size:
            raise ValueError("Cache must hold at least one SQLite-sized chunk")
        self.source = source
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.directory = self.directory.resolve()
        self.scratch = scratch or ScratchBudget(self.directory, cache_limit + 32 * 1024 * 1024)
        self.scratch.require(4096)
        self.chunk_size = chunk_size
        self.read_ahead = chunk_size if read_ahead is None else read_ahead
        if not chunk_size <= self.read_ahead <= 16 * 1024 * 1024 or self.read_ahead % chunk_size:
            raise ValueError("Read ahead must be a chunk multiple up to 16 MiB")
        self.cache_limit = cache_limit
        self.hits = 0
        self.read_calls = 0
        self.identity = asdict(source.identity)
        self.manifest_path = self.directory / "chunks.json"
        self.ledger = None
        self.lock = (self.directory / ".lock").open("a+b")
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.lock.close()
            self.lock = None
            raise ValueError("Cache is already in use by another reader") from None
        try:
            self._initialize(prior_receipts)
        except BaseException:
            self.close()
            raise

    def _initialize(self, prior_receipts):
        ledger_path = self.directory / "transfers.sqlite"
        if self.manifest_path.exists() and not ledger_path.exists() and not prior_receipts:
            raise ValueError("Legacy cache requires its prior transfer receipts before reuse")
        self.chunks = {}
        self.pending = None
        self.source.check()
        if self.manifest_path.exists():
            if self.manifest_path.stat().st_size > 16 * 1024 * 1024:
                raise LimitExceeded("Cache manifest exceeds metadata size limit")
            manifest = json.loads(self.manifest_path.read_text())
            previous_size = manifest["chunk_size"]
            if (manifest["identity"] != self.identity or previous_size < self.chunk_size
                    or previous_size % self.chunk_size):
                raise SourceChanged("Cache belongs to another source identity or chunk size")
            self.chunks = manifest["chunks"]
            self.pending = manifest.get("pending")
        for key, entry in self.chunks.items():
            self._validate_entry(key, entry)
        previous_end = 0
        for key in sorted(self.chunks, key=int):
            if int(key) < previous_end:
                raise ValueError("Overlapping cached chunks")
            previous_end = int(key) + self.chunks[key]["length"]
        if self.pending is not None:
            self._validate_entry(self.pending["key"], self.pending["entry"])
            begin = int(self.pending["key"])
            end = begin + self.pending["entry"]["length"]
            if any(begin < int(key) + entry["length"] and int(key) < end
                   for key, entry in self.chunks.items()):
                raise ValueError("Pending cache entry overlaps a committed chunk")
        self.scratch.reserve_sqlite(ledger_path, 4 * 1024 * 1024)
        self.ledger = TransferLedger(ledger_path, self.identity, prior_receipts=prior_receipts)
        self.source.ledger = self.ledger
        expected = {"chunks.json", "chunks.json.tmp", "chunks.tmp", "receipt.json", "transfers.sqlite", ".lock"}
        expected.update(str(k) + ".bin" for k in self.chunks)
        if self.pending is not None:
            expected.update({self.pending["key"] + ".bin", self.pending["key"] + ".part"})
        if any(p.name not in expected for p in self.directory.iterdir()):
            raise ValueError("Untracked cache files; inspect the interrupted extraction before reuse")
        if any(p.is_symlink() for p in self.directory.iterdir()):
            raise ValueError("Cache files must not be symbolic links")
        self.header_sha256 = manifest["header_sha256"] if self.manifest_path.exists() else None
        self._recover_pending()
        self.chunk_offsets = sorted(map(int, self.chunks))
        self.cached_bytes = sum(item["length"] for item in self.chunks.values())
        if self.cached_bytes > self.cache_limit:
            raise LimitExceeded("Existing cache exceeds cache limit")
        # Identity includes a stable content anchor in addition to SIZE/MDTM.
        anchor = self.source.read(0, min(4096, self.source.identity.size))
        self.header_sha256 = sha256(anchor)
        if self.manifest_path.exists() and manifest["header_sha256"] != self.header_sha256:
            raise SourceChanged("Remote header changed despite matching metadata")
        self._save()

    def _validate_entry(self, key, entry):
        if not isinstance(key, str) or not re.fullmatch(r"0|[1-9][0-9]*", key):
            raise ValueError("Invalid cache chunk offset")
        offset = int(key)
        if offset >= self.source.identity.size or offset % self.chunk_size:
            raise ValueError("Invalid cache chunk alignment or bounds")
        length = entry.get("length")
        if (type(length) is not int or not 0 < length <= 16 * 1024 * 1024
                or offset + length > self.source.identity.size
                or (length % self.chunk_size and offset + length != self.source.identity.size)
                or not isinstance(entry.get("sha256"), str)
                or not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"])):
            raise ValueError("Invalid cache chunk length or digest")

    def _recover_pending(self):
        if self.pending is None:
            return
        key, entry = self.pending["key"], self.pending["entry"]
        path = self.directory / (key + ".bin")
        temporary = self.directory / (key + ".part")
        candidate = path if path.exists() else temporary
        if candidate.exists():
            data = candidate.read_bytes() if candidate.stat().st_size == entry["length"] else b""
            if len(data) == entry["length"] and sha256(data) == entry["sha256"]:
                candidate.replace(path)
                self.chunks[key] = entry
            elif candidate == path:
                raise ValueError("Interrupted cache commit has a corrupt final chunk")
        temporary.unlink(missing_ok=True)
        self.pending = None
        self._save()

    def _save(self):
        manifest = dict(identity=self.identity, header_sha256=self.header_sha256,
                        chunk_size=self.chunk_size, chunks=self.chunks, pending=self.pending)
        data = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode()
        if len(data) > 16 * 1024 * 1024:
            raise LimitExceeded("Cache manifest exceeds metadata size limit")
        self.scratch.atomic_write(self.manifest_path, data)

    def _chunk(self, offset):
        key = str(offset)
        path = self.directory / (key + ".bin")
        if key in self.chunks:
            entry = self.chunks[key]
            data = path.read_bytes() if path.stat().st_size == entry["length"] else b""
            if len(data) != entry["length"] or sha256(data) != entry["sha256"]:
                raise ValueError("Cached chunk checksum/length mismatch")
            self.hits += 1
            return data
        length = min(self.read_ahead, self.source.identity.size - offset)
        following = bisect_right(self.chunk_offsets, offset)
        if following < len(self.chunk_offsets):
            length = min(length, self.chunk_offsets[following] - offset)
        # Reserve space before the network operation; do not evict and re-fetch
        # indefinitely under an apparently small disk budget.
        if self.cached_bytes + length > self.cache_limit:
            raise LimitExceeded("Cache limit exceeded")
        # Leave space for the chunk plus both manifest generations before any
        # acquisition. Actual metadata size is checked again at each write.
        metadata_size = self.manifest_path.stat().st_size if self.manifest_path.exists() else 0
        self.scratch.require(length + 2 * metadata_size + 16384)
        data = self.source.read(offset, length)
        shift = 17
        if length > shift:
            overlap = self.source.read(offset + shift, length - shift)
            if data[shift:] != overlap:
                raise ValueError("Shifted range overlap mismatch: REST offset not reliable")
        if offset == 0 and sha256(data[:min(4096, length)]) != self.header_sha256:
            raise SourceChanged("Remote header changed during chunk acquisition")
        entry = dict(length=length, sha256=sha256(data),
                     verification="independent REST +17 full overlap")
        self.pending = dict(key=key, entry=entry)
        self._save()
        temporary = path.with_suffix(".part")
        with temporary.open("wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        self.scratch.usage()
        temporary.replace(path)
        self.chunks[key] = entry
        insort(self.chunk_offsets, offset)
        self.pending = None
        self.cached_bytes += length
        self._save()
        return data

    def read(self, offset, length):
        self.read_calls += 1
        if time.monotonic() - self.source.started > self.source.seconds_limit:
            raise LimitExceeded("Elapsed-time limit exceeded during SQLite read")
        if offset < 0 or length <= 0 or offset + length > self.source.identity.size:
            raise ValueError("SQLite requested bytes outside the deposited file")
        output = []
        end = offset + length
        while offset < end:
            base = offset // self.chunk_size * self.chunk_size
            index = bisect_right(self.chunk_offsets, offset) - 1
            if index >= 0:
                cached_start = self.chunk_offsets[index]
                if offset < cached_start + self.chunks[str(cached_start)]["length"]:
                    base = cached_start
            chunk = self._chunk(base)
            stop = min(end - base, len(chunk))
            output.append(chunk[offset - base:stop])
            offset = base + stop
        return b"".join(output)

    def receipt(self):
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        versions = dict(python=sys.version, platform=sys.platform,
                        ftps_sqlite_sha256=CODE_SHA256, extraction_budget_sha256=BUDGET_CODE_SHA256)
        for name in ("apsw", "orfdata"):
            try:
                versions[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                versions[name] = None
        if "apsw" in sys.modules:
            versions["sqlite"] = sys.modules["apsw"].sqlitelibversion()
        return dict(identity=self.identity, header_sha256=self.header_sha256,
                    software=versions,
                    integrity_scope="metadata/header pin plus verified acquired chunks; no full-file checksum",
                    payload_bytes=self.source.payload_bytes, range_requests=self.source.requests,
                    cumulative_transfer=self.ledger.summary(),
                    scratch=self.scratch.receipt(),
                    failures=self.source.failures, cached_payload_bytes=self.cached_bytes,
                    disk_bytes=sum(p.stat().st_size for p in self.directory.iterdir() if p.is_file()),
                    peak_process_rss_bytes=peak if sys.platform == "darwin" else peak * 1024,
                    elapsed_seconds=time.monotonic() - self.source.started,
                    cache_hits=self.hits, sqlite_reads=self.read_calls,
                    transfer_limit=self.source.transfer_limit, cache_limit=self.cache_limit,
                    request_limit=self.source.request_limit, seconds_limit=self.source.seconds_limit,
                    chunk_size=self.chunk_size, read_ahead=self.read_ahead,
                    reuse_control=self.source.reuse_control,
                    control_connections=self.source.control_connections,
                    bounded_close_responses=self.source.bounded_close_responses,
                    accounting="recv data payload, including overlap/retries; excludes network/TLS overhead and server read-ahead")

    def close(self):
        self.source.close()
        if self.ledger is not None:
            self.ledger.close()
            if self.source.ledger is self.ledger:
                self.source.ledger = None
            self.ledger = None
        if self.lock is not None:
            self.lock.close()
            self.lock = None


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
