"""Bounded source access: independent SQLite traversal and hostile transports."""

import json
import sqlite3
from dataclasses import asdict

import pytest

from tools.ftps_sqlite import (
    ChunkCache,
    FTPSRanges,
    Identity,
    LimitExceeded,
    SourceChanged,
    connect,
)


class LocalSource(FTPSRanges):
    def __init__(self, data, **kwargs):
        super().__init__(Identity("fixture", "/db", len(data), "213 pinned"), **kwargs)
        self.data = data
        self.changed = False
        self.ignore_offset = False

    def check(self, ftp=None):
        if self.changed:
            raise SourceChanged("Fixture identity changed")

    def _connect(self):
        source = self

        class Stream:
            def __init__(self, offset):
                self.offset = offset

            def recv(self, length):
                data = source.data[self.offset:self.offset + length]
                self.offset += len(data)
                return data

            def close(self):
                pass

        class FTP:
            def sendcmd(self, command):
                self.offset = 0 if source.ignore_offset else int(command.split()[1])
                return "350 Restart accepted"

            def transfercmd(self, command):
                return Stream(self.offset)

            def close(self):
                pass

        return FTP()


@pytest.fixture
def sqlite_bytes(tmp_path):
    path = tmp_path / "local.sqlite"
    connection = sqlite3.connect(path)
    connection.execute("pragma page_size=4096")
    connection.execute("create table peptides(id integer primary key, sequence text, payload blob)")
    rows = [(i, f"PEPTIDE{i:05d}", bytes([i % 251]) * (20000 if i % 37 == 0 else 100))
            for i in range(800)]
    connection.executemany("insert into peptides values (?, ?, ?)", rows)
    connection.execute("create index sequences on peptides(sequence)")
    connection.commit()
    connection.close()
    return path.read_bytes(), rows


def test_sqlite_index_interior_overflow_and_cross_chunk_reads(sqlite_bytes, tmp_path):
    pytest.importorskip("apsw")
    data, rows = sqlite_bytes
    source = LocalSource(data, transfer_limit=4 * len(data))
    cache = ChunkCache(source, tmp_path / "cache", chunk_size=8192, cache_limit=2 * len(data))
    db, vfs = connect(cache)
    try:
        plan = list(db.execute("explain query plan select id,payload from peptides where sequence=?",
                               ("PEPTIDE00370",)))
        assert "INDEX sequences" in str(plan)
        assert list(db.execute("select id,payload from peptides where sequence=?", ("PEPTIDE00370",))) == [
            (370, rows[370][2])]
        assert list(db.execute("select * from peptides order by id")) == rows
        assert cache.read(4090, 8220) == data[4090:12310]
        with pytest.raises(Exception, match="readonly|read-only"):
            db.execute("delete from peptides")
    finally:
        db.close()
        vfs.unregister()
    assert cache.cached_bytes <= len(data)
    assert source.payload_bytes < 3 * len(data)
    cache.close()


def test_ignored_offset_is_rejected_before_cache_commit(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    source = LocalSource(data, transfer_limit=100000)
    cache = ChunkCache(source, tmp_path / "cache", chunk_size=8192)
    source.ignore_offset = True
    with pytest.raises(ValueError, match="overlap mismatch"):
        cache.read(8192, 4096)
    assert not cache.chunks and not list((tmp_path / "cache").glob("*.bin"))
    cache.close()


def test_source_identity_header_and_cache_corruption(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    directory = tmp_path / "cache"
    source = LocalSource(data, transfer_limit=100000)
    cache = ChunkCache(source, directory, chunk_size=8192)
    assert cache.read(0, 100) == data[:100]
    (directory / "0.bin").write_bytes(b"wrong")
    with pytest.raises(ValueError, match="checksum/length"):
        cache.read(0, 100)
    cache.close()
    source.changed = True
    with pytest.raises(SourceChanged):
        ChunkCache(source, directory, chunk_size=8192)
    source.changed = False
    source.data = b"changed" + data[7:]
    with pytest.raises(SourceChanged, match="header changed"):
        ChunkCache(source, directory, chunk_size=8192)


def test_budgets_bounds_and_resumption(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    directory = tmp_path / "cache"
    source = LocalSource(data, transfer_limit=21000)
    cache = ChunkCache(source, directory, chunk_size=8192, cache_limit=8192)
    assert cache.read(0, 100) == data[:100]
    with pytest.raises(LimitExceeded, match="Cache limit"):
        cache.read(8192, 10)
    with pytest.raises(ValueError, match="outside"):
        cache.read(len(data) - 1, 2)
    used = cache.receipt()["cumulative_transfer"]["charged_payload_bytes"]
    assert used == source.payload_bytes
    cache.close()
    with pytest.raises(LimitExceeded, match="Transfer"):
        ChunkCache(LocalSource(data, transfer_limit=used + 4095), directory, chunk_size=8192)
    resumed = ChunkCache(LocalSource(data, transfer_limit=used + 4096), directory,
                         chunk_size=8192, cache_limit=8192)
    assert resumed.read(0, 100) == data[:100]
    assert resumed.source.payload_bytes == 4096
    assert resumed.receipt()["cumulative_transfer"]["charged_payload_bytes"] == used + 4096
    resumed.close()
    with pytest.raises(LimitExceeded, match="Transfer"):
        ChunkCache(LocalSource(data, transfer_limit=4095), tmp_path / "small")


def test_short_read_retries_count_all_payload_and_recheck_identity(monkeypatch, tmp_path):
    from tools.extraction_budget import TransferLedger

    source = FTPSRanges(Identity("fixture", "/db", 10, "213 pinned"), transfer_limit=30)
    source.ledger = TransferLedger(tmp_path / "transfers.sqlite", asdict(source.identity))
    calls = []

    class Stream:
        def __init__(self, short):
            self.data = b"12" if short else b"1234567890"

        def recv(self, length):
            data, self.data = self.data[:length], self.data[length:]
            return data

        def close(self):
            pass

    class FTP:
        def size(self, path):
            return 10

        def sendcmd(self, command):
            return "213 pinned" if command.startswith("MDTM") else "350 Restart accepted"

        def transfercmd(self, command):
            calls.append(command)
            return Stream(short=len(calls) == 1)

        def close(self):
            pass

    monkeypatch.setattr(source, "_connect", FTP)
    assert source.read(0, 10) == b"1234567890"
    assert source.payload_bytes == 12 and source.requests == 2 and len(source.failures) == 1
    assert source.ledger.summary()["verified_received_payload_bytes"] == 12
    assert source.ledger.summary()["charged_payload_bytes"] == 12
    assert list(source.ledger.connection.execute("SELECT received,status FROM requests ORDER BY id")) == [
        (2, "failed"), (10, "complete")]
    source.transfer_limit = 15
    with pytest.raises(LimitExceeded, match="Transfer"):
        source.read(0, 10)
    assert len(calls) == 2
    source.ledger.close()


def test_remote_metadata_change_is_not_retried(monkeypatch):
    source = FTPSRanges(Identity("fixture", "/db", 10, "213 pinned"), transfer_limit=100)

    class FTP:
        def size(self, path):
            return 11

        def sendcmd(self, command):
            return "213 changed"

        def close(self):
            pass

    monkeypatch.setattr(source, "_connect", FTP)
    with pytest.raises(SourceChanged):
        source.read(0, 10)
    assert source.requests == 1 and source.payload_bytes == 0 and not source.failures


def test_cache_lock_and_invalid_manifest_offsets(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    cache = ChunkCache(LocalSource(data, transfer_limit=100000), tmp_path / "cache", chunk_size=8192)
    with pytest.raises(ValueError, match="already in use"):
        ChunkCache(LocalSource(data, transfer_limit=100000), cache.directory, chunk_size=8192)
    cache.read(0, 100)
    cache.close()
    manifest = json.loads(cache.manifest_path.read_text())
    manifest["chunks"]["../outside"] = manifest["chunks"].pop("0")
    cache.manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="chunk offset"):
        ChunkCache(LocalSource(data, transfer_limit=100000), cache.directory, chunk_size=8192)


@pytest.mark.parametrize("stage", ["before_write", "partial_write", "after_write", "after_rename"])
def test_interrupted_chunk_commits_resume_without_losing_transfer_charge(sqlite_bytes, tmp_path, stage):
    from tools.ftps_sqlite import sha256

    data, _ = sqlite_bytes
    cache = ChunkCache(LocalSource(data, transfer_limit=100000), tmp_path / "cache", chunk_size=8192)
    chunk = cache.source.read(8192, 8192)
    overlap = cache.source.read(8209, 8175)
    assert chunk[17:] == overlap
    cache.pending = dict(key="8192", entry=dict(length=8192, sha256=sha256(chunk)))
    cache._save()
    if stage != "before_write":
        (cache.directory / "8192.part").write_bytes(chunk[:100] if stage == "partial_write" else chunk)
    if stage == "after_rename":
        (cache.directory / "8192.part").replace(cache.directory / "8192.bin")
    charged = cache.receipt()["cumulative_transfer"]["charged_payload_bytes"]
    cache.close()
    resumed = ChunkCache(LocalSource(data, transfer_limit=100000), cache.directory, chunk_size=8192)
    assert resumed.read(8192, 8192) == data[8192:16384]
    assert resumed.source.payload_bytes == (4096 if stage in ("after_write", "after_rename") else 20463)
    assert resumed.receipt()["cumulative_transfer"]["charged_payload_bytes"] == charged + resumed.source.payload_bytes
    assert resumed.pending is None and not list(resumed.directory.glob("*.part"))
    resumed.close()


def test_smaller_index_reads_reuse_larger_verified_chunks(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    cache = ChunkCache(LocalSource(data, transfer_limit=100000), tmp_path / "cache", chunk_size=16384)
    assert cache.read(0, 100) == data[:100]
    cache.close()
    resumed = ChunkCache(LocalSource(data, transfer_limit=100000), cache.directory, chunk_size=4096)
    assert resumed.read(8192, 10) == data[8192:8202]
    assert resumed.source.payload_bytes == 4096  # Only the header pin is read again.
    assert resumed.read(16380, 20) == data[16380:16400]
    assert resumed.source.payload_bytes == 4096 + 4096 + 4079
    assert resumed.cached_bytes == 16384 + 4096
    resumed.close()


def test_cache_scratch_limit_prevents_network_acquisition(sqlite_bytes, tmp_path):
    from tools.extraction_budget import ScratchBudget

    data, _ = sqlite_bytes
    scratch = ScratchBudget(tmp_path / "cache", 13 * 1024 * 1024)
    source = LocalSource(data, transfer_limit=100000)
    cache = ChunkCache(source, scratch.directory, chunk_size=8192, scratch=scratch)
    scratch.limit = scratch.receipt()["committed_capacity_bytes"] + 4096
    with pytest.raises(LimitExceeded, match="Scratch"):
        cache.read(0, 100)
    assert source.payload_bytes == 4096
    cache.close()
