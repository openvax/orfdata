"""Bounded source access: independent SQLite traversal and hostile transports."""

import sqlite3

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

    def read(self, offset, length):
        self._check_budget(length)
        self.check()
        self.requests += 1
        self.payload_bytes += length
        if self.ignore_offset:
            offset = 0
        return self.data[offset:offset + length]


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


def test_ignored_offset_is_rejected_before_cache_commit(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    source = LocalSource(data, transfer_limit=100000)
    cache = ChunkCache(source, tmp_path / "cache", chunk_size=8192)
    source.ignore_offset = True
    with pytest.raises(ValueError, match="overlap mismatch"):
        cache.read(8192, 4096)
    assert not cache.chunks and not list((tmp_path / "cache").glob("*.bin"))


def test_source_identity_header_and_cache_corruption(sqlite_bytes, tmp_path):
    data, _ = sqlite_bytes
    directory = tmp_path / "cache"
    source = LocalSource(data, transfer_limit=100000)
    cache = ChunkCache(source, directory, chunk_size=8192)
    assert cache.read(0, 100) == data[:100]
    (directory / "0.bin").write_bytes(b"wrong")
    with pytest.raises(ValueError, match="checksum/length"):
        cache.read(0, 100)
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
    resumed = ChunkCache(LocalSource(data, transfer_limit=4096), directory,
                         chunk_size=8192, cache_limit=8192)
    assert resumed.read(0, 100) == data[:100]
    assert resumed.source.payload_bytes == 4096
    with pytest.raises(LimitExceeded, match="Transfer"):
        ChunkCache(LocalSource(data, transfer_limit=4095), tmp_path / "small")


def test_short_read_retries_count_all_payload_and_recheck_identity(monkeypatch):
    source = FTPSRanges(Identity("fixture", "/db", 10, "213 pinned"), transfer_limit=30)
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
    source.transfer_limit = 15
    with pytest.raises(LimitExceeded, match="Transfer"):
        source.read(0, 10)
    assert len(calls) == 2


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
