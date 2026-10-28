"""Durable accounting for bounded acquisition, including interrupted processes."""

import hashlib
import json
import os
import sqlite3
from pathlib import Path

CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


class LimitExceeded(ValueError):
    """An explicit extraction resource limit would be exceeded."""


class ScratchBudget:
    """Bound owned scratch/output files, including atomic-write duplicates.

    SQLite files have fixed reservations for the database and rollback journal;
    their connections must separately set max_page_count. Other writes reserve
    their complete temporary file before writing. One extractor owns this tree.
    Both logical bytes and allocated file blocks count toward the limit.
    """

    def __init__(self, directory, limit):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        if type(limit) is not int or limit <= 0:
            raise ValueError("Scratch limit must be a positive byte count")
        self.limit = limit
        self.reservations = {}
        self.peak_bytes = 0
        self.require(0)

    def _path(self, path):
        path = Path(path).absolute()
        if not path.is_relative_to(self.directory) or path.is_symlink():
            raise ValueError("Output path is outside the owned scratch directory")
        if path.parent.resolve() != path.parent:
            raise ValueError("Scratch output parents must not contain symbolic links")
        return path

    def reserve_sqlite(self, path, maximum):
        path = self._path(path)
        # DELETE journals store at most one original copy of every page, plus
        # per-page and sector headers. Three database capacities cover both.
        self.reservations[path] = 3 * maximum
        self.require(0)

    def usage(self):
        total, unreserved = 0, 0
        for path in self.directory.rglob("*"):
            if path.is_symlink():
                raise ValueError("Scratch tree must not contain symbolic links")
            if not path.is_file():
                continue
            stat = path.stat()
            size = max(stat.st_size, stat.st_blocks * 512)
            total += size
            if not any(path == name or path == Path(str(name) + "-journal")
                       for name in self.reservations):
                unreserved += size
        self.peak_bytes = max(self.peak_bytes, total)
        return total, unreserved + sum(self.reservations.values())

    def require(self, extra):
        actual, reserved = self.usage()
        if max(actual, reserved) + extra > self.limit:
            raise LimitExceeded("Scratch/output disk limit exceeded")

    def atomic_write(self, path, data):
        path = self._path(path)
        temporary = path.with_name(path.name + ".tmp")
        self._path(temporary)
        unit = max(4096, path.parent.stat().st_blksize)
        self.require(((len(data) + unit - 1) // unit) * unit)
        with temporary.open("wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        self.usage()  # Observe the duplicate before replacing its predecessor.
        temporary.replace(path)
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)

    def receipt(self):
        actual, reserved = self.usage()
        return dict(limit_bytes=self.limit, current_bytes=actual, peak_observed_bytes=self.peak_bytes,
                    committed_capacity_bytes=reserved,
                    accounting="max(logical size, allocated file blocks), including atomic temporaries; "
                               "SQLite database/journal capacity reserved in advance")


class TransferLedger:
    """Reserve each request before opening it; settle actual received payload.

    Unsettled requests retain their full reservation after a crash. They must
    not silently restore the budget: the process cannot prove how many bytes
    reached its receive buffer. All completed attempts, including retries, are
    retained. Transactions serialize reservations across connections.
    """

    def __init__(self, path, identity, *, prior_receipts=(), disk_limit=4 * 1024 * 1024):
        self.path = Path(path)
        self.connection = sqlite3.connect(self.path, timeout=30)
        try:
            self._initialize(identity, prior_receipts, disk_limit)
        except BaseException:
            self.close()
            raise

    def _initialize(self, identity, prior_receipts, disk_limit):
        self.connection.execute("PRAGMA page_size=4096")
        pages = disk_limit // 4096
        actual = self.connection.execute(f"PRAGMA max_page_count={pages}").fetchone()[0]
        if actual > pages:
            raise LimitExceeded("Existing transfer ledger exceeds its disk limit")
        self.connection.execute("PRAGMA synchronous=FULL")
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 1024 * 1024)
        self.connection.executescript("""
            CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS usage (
                id INTEGER PRIMARY KEY CHECK(id=1), charged INTEGER NOT NULL,
                received INTEGER NOT NULL, requests INTEGER NOT NULL);
            CREATE TABLE IF NOT EXISTS requests (
                id INTEGER PRIMARY KEY, offset INTEGER NOT NULL, length INTEGER NOT NULL,
                received INTEGER, status TEXT NOT NULL DEFAULT 'reserved');
        """)
        encoded = json.dumps(identity, sort_keys=True)
        row = self.connection.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
        if row is not None and row[0] != encoded:
            raise ValueError("Transfer ledger belongs to another source identity")
        if row is None:
            total_bytes, total_requests, provenance = 0, 0, []
            for receipt_path in prior_receipts:
                if Path(receipt_path).stat().st_size > 4 * 1024 * 1024:
                    raise LimitExceeded("Prior receipt exceeds metadata size limit")
                data = Path(receipt_path).read_bytes()
                receipt = json.loads(data)
                if receipt.get("identity") != identity:
                    raise ValueError("Prior transfer receipt belongs to another source identity")
                size, count = receipt["payload_bytes"], receipt["range_requests"]
                if type(size) is not int or size < 0 or type(count) is not int or count < 0:
                    raise ValueError("Invalid prior transfer counters")
                digest = hashlib.sha256(data).hexdigest()
                if digest in {p["sha256"] for p in provenance}:
                    raise ValueError("Duplicate prior transfer receipt")
                provenance.append(dict(sha256=digest, payload_bytes=size, range_requests=count))
                total_bytes += size
                total_requests += count
            with self.connection:
                self.connection.execute("INSERT INTO metadata VALUES ('identity',?)", (encoded,))
                self.connection.execute("INSERT INTO metadata VALUES ('prior_receipts',?)",
                                        (json.dumps(provenance, sort_keys=True),))
                self.connection.execute("INSERT INTO usage VALUES (1,?,?,?)",
                                        (total_bytes, total_bytes, total_requests))

    def reserve(self, offset, length, *, byte_limit, request_limit):
        if type(offset) is not int or offset < 0 or type(length) is not int or length <= 0:
            raise ValueError("Invalid transfer reservation")
        self.connection.execute("BEGIN IMMEDIATE")
        try:
            charged, _, count = self.connection.execute("SELECT charged,received,requests FROM usage").fetchone()
            if charged + length > byte_limit:
                raise LimitExceeded("Cumulative transfer payload limit exceeded")
            if count >= request_limit:
                raise LimitExceeded("Cumulative range request limit exceeded")
            cursor = self.connection.execute("INSERT INTO requests(offset,length) VALUES (?,?)", (offset, length))
            self.connection.execute("UPDATE usage SET charged=charged+?, requests=requests+1", (length,))
            self.connection.commit()
            return cursor.lastrowid
        except BaseException:
            self.connection.rollback()
            raise

    def settle(self, request, received, *, status):
        self.connection.execute("BEGIN IMMEDIATE")
        try:
            row = self.connection.execute("SELECT length,received FROM requests WHERE id=?", (request,)).fetchone()
            if row is None or row[1] is not None:
                raise ValueError("Unknown or already settled transfer request")
            if type(received) is not int or not 0 <= received <= row[0] or status == "reserved":
                raise ValueError("Invalid transfer settlement")
            self.connection.execute("UPDATE requests SET received=?,status=? WHERE id=?", (received, status, request))
            self.connection.execute("UPDATE usage SET charged=charged-?,received=received+?",
                                    (row[0] - received, received))
            self.connection.commit()
        except BaseException:
            self.connection.rollback()
            raise

    def summary(self):
        charged, received, count = self.connection.execute("SELECT charged,received,requests FROM usage").fetchone()
        unsettled = self.connection.execute("SELECT count(*) FROM requests WHERE received IS NULL").fetchone()[0]
        return dict(charged_payload_bytes=charged, verified_received_payload_bytes=received,
                    range_requests=count, unsettled_requests=unsettled,
                    uncertain_reserved_bytes=charged - received)

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
