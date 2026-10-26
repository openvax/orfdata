"""Resumable selection of native SQLite facts, without spectra or a file copy."""

import hashlib
import json
import math
import re
import sqlite3
import struct
from pathlib import Path

from tools.extraction_budget import LimitExceeded

CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def quoted(identifier):
    if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", identifier):
        raise ValueError("Invalid native table/column identifier")
    return '"' + identifier + '"'


def canonical_row(row):
    """Unambiguous typed bytes; integer IDs and floating point bits stay exact."""
    result = []
    for value in row:
        if value is None:
            tag, payload = b"n", b""
        elif type(value) is int:
            tag, payload = b"i", struct.pack(">q", value)
        elif type(value) is float:
            if math.isnan(value):
                raise ValueError("Scalar NaN cannot round-trip through SQLite; retain its source BLOB")
            tag, payload = b"f", struct.pack(">d", value)
        elif isinstance(value, str):
            tag, payload = b"s", value.encode("utf-8")
        elif isinstance(value, bytes):
            tag, payload = b"b", value
        else:
            raise TypeError("Unsupported native SQLite value")
        result.extend((tag, struct.pack(">Q", len(payload)), payload))
    return b"".join(result)


def source_schema(database):
    fields = ("type", "name", "table", "rootpage", "sql")
    rows = [dict(zip(fields, row)) for row in database.execute(
        "SELECT type,name,tbl_name,rootpage,sql FROM sqlite_master ORDER BY type,name")]
    data = (json.dumps(rows, sort_keys=True, indent=2) + "\n").encode()
    return rows, hashlib.sha256(data).hexdigest()


class NativeExport:
    """Checkpoint each output batch atomically with its source cursor and digest.

    Scans retain every selected source row, including negative/noncontiguous
    row IDs. Referenced acquisitions use bounded indexed key lookups, driven by
    previously exported native link tables. Missing links remain explicit.
    The output is staging evidence; completion of this selection does not assert
    normalization, reference mapping or biological acceptance of the dataset.
    """

    def __init__(self, source, output, plan, identity, scratch, *, output_limit,
                 batch_rows=500, batch_bytes=16 * 1024 * 1024, lookup_keys=32):
        if min(batch_rows, batch_bytes, lookup_keys) < 1:
            raise ValueError("Export batch limits must be positive")
        self.source, self.plan, self.scratch = source, plan, scratch
        self.path = Path(output).resolve()
        self.batch_rows, self.batch_bytes, self.lookup_keys = batch_rows, batch_bytes, lookup_keys
        self.tables = {item["table"]: item for item in plan["tables"]}
        if len(self.tables) != len(plan["tables"]):
            raise ValueError("Duplicate table in native selection")
        _, schema_hash = source_schema(source)
        if plan["source_schema_sha256"] != schema_hash:
            raise ValueError("Source schema differs from the inspected extraction plan")
        for item in plan["tables"]:
            table = item["table"]
            columns = item["columns"]
            if not columns or len(columns) != len(set(columns)):
                raise ValueError("Native selection needs unique columns")
            available = {row[1] for row in source.execute(f"PRAGMA table_info({quoted(table)})")}
            if not set(columns) <= available:
                raise ValueError("Selected native column is absent: " + table)
            for column in columns:
                quoted(column)
            for selector in item.get("select_from", []):
                if (selector["table"] not in self.tables or not set(selector["columns"]) <=
                        set(self.tables[selector["table"]]["columns"])):
                    raise ValueError("Acquisition selector requires an exported native link")
                if len(selector["columns"]) != len(item["keys"]):
                    raise ValueError("Acquisition selector key arity differs")
            if not set(item.get("keys", [])) <= set(columns):
                raise ValueError("Selected acquisition must retain its lookup keys")
        scratch.reserve_sqlite(self.path, output_limit)
        self.database = sqlite3.connect(self.path)
        try:
            self._initialize(identity, output_limit)
        except BaseException:
            self.close()
            raise

    def _initialize(self, identity, output_limit):
        db = self.database
        db.execute("PRAGMA page_size=4096")
        if db.execute(f"PRAGMA max_page_count={output_limit // 4096}").fetchone()[0] > output_limit // 4096:
            raise LimitExceeded("Existing native export exceeds output limit")
        db.execute("PRAGMA synchronous=FULL")
        db.execute("PRAGMA journal_mode=DELETE")
        db.execute("PRAGMA temp_store=MEMORY")
        db.execute("PRAGMA cache_size=-8192")
        db.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
        db.executescript("""
            CREATE TABLE IF NOT EXISTS _export_metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS _export_progress(
                name TEXT PRIMARY KEY, cursor TEXT, rows INTEGER NOT NULL,
                selected_keys INTEGER NOT NULL, missing_keys INTEGER NOT NULL,
                digest TEXT NOT NULL, complete INTEGER NOT NULL);
            CREATE TABLE IF NOT EXISTS _export_missing(
                name TEXT NOT NULL, source_key TEXT NOT NULL,
                reason TEXT NOT NULL, PRIMARY KEY(name,source_key));
        """)
        configuration = json.dumps(dict(plan=self.plan, identity=identity), sort_keys=True, separators=(",", ":"))
        previous = db.execute("SELECT value FROM _export_metadata WHERE key='configuration'").fetchone()
        if previous is not None and previous[0] != configuration:
            raise ValueError("Export checkpoint belongs to another source or selection")
        with db:
            db.execute("INSERT OR IGNORE INTO _export_metadata VALUES ('configuration',?)", (configuration,))
            for item in self.plan["tables"]:
                columns = ",".join(quoted(c) for c in item["columns"])
                # Typeless native columns preserve SQLite storage classes and int64 IDs.
                db.execute(f"CREATE TABLE IF NOT EXISTS {quoted(item['table'])}"
                           f"(_source_rowid INTEGER PRIMARY KEY,{columns})")
                db.execute("INSERT OR IGNORE INTO _export_progress VALUES (?,NULL,0,0,0,?,0)",
                           (item["table"], "0" * 64))

    def progress(self, table):
        row = self.database.execute("SELECT cursor,rows,selected_keys,missing_keys,digest,complete "
                                    "FROM _export_progress WHERE name=?", (table,)).fetchone()
        result = dict(zip(("cursor", "rows", "selected_keys", "missing_keys", "digest", "complete"), row))
        result["cursor"] = json.loads(result["cursor"]) if result["cursor"] is not None else None
        return result

    def _commit_batch(self, item, rows, cursor, *, selected=0, missing=(), complete=False):
        name = item["table"]
        previous = self.progress(name)
        digest = bytes.fromhex(previous["digest"])
        for row in rows:
            digest = hashlib.sha256(digest + canonical_row(row)).digest()
        columns = "_source_rowid," + ",".join(quoted(c) for c in item["columns"])
        placeholders = ",".join("?" for _ in range(len(item["columns"]) + 1))
        with self.database:
            self.database.executemany(f"INSERT INTO {quoted(name)}({columns}) VALUES ({placeholders})", rows)
            self.database.executemany("INSERT INTO _export_missing VALUES (?,?,?)", (
                (name, json.dumps(key, separators=(",", ":")), "No native row for the referenced key")
                for key in missing))
            self.database.execute("UPDATE _export_progress SET cursor=?,rows=?,selected_keys=?,"
                                  "missing_keys=?,digest=?,complete=? WHERE name=?",
                                  (json.dumps(cursor), previous["rows"] + len(rows),
                                   previous["selected_keys"] + selected, previous["missing_keys"] + len(missing),
                                   digest.hex(), int(complete), name))
        self.scratch.usage()

    def _scan(self, item):
        state = self.progress(item["table"])
        columns = ",".join(quoted(c) for c in item["columns"])
        where, parameters = ("", ()) if state["cursor"] is None else ("WHERE rowid>?", (state["cursor"],))
        cursor = self.source.execute(f"SELECT rowid,{columns} FROM {quoted(item['table'])} "
                                     f"{where} ORDER BY rowid LIMIT ?", (*parameters, self.batch_rows))
        rows, size = [], 0
        for row in cursor:
            row = tuple(row)
            encoded = canonical_row(row)
            if len(encoded) > self.batch_bytes:
                raise LimitExceeded("Native row exceeds the export batch byte limit")
            if size + len(encoded) > self.batch_bytes:
                break
            rows.append(row)
            size += len(encoded)
        self._commit_batch(item, rows, rows[-1][0] if rows else state["cursor"], complete=not rows)

    def _queue(self, item):
        name = item["table"]
        queue = '"_selection_' + name + '"'
        keys = ["k" + str(i) for i in range(len(item["keys"]))]
        if self.database.execute("SELECT 1 FROM sqlite_master WHERE name=?", ("_selection_" + name,)).fetchone():
            return queue, keys
        for selector in item["select_from"]:
            if not self.progress(selector["table"])["complete"]:
                raise ValueError("Complete native link export before acquisition selection: " + selector["table"])
        definitions = ",".join(f"{k} INTEGER NOT NULL CHECK(typeof({k})='integer')" for k in keys)
        with self.database:
            self.database.execute(f"CREATE TABLE {queue}({definitions},PRIMARY KEY({','.join(keys)})) WITHOUT ROWID")
            for selector in item["select_from"]:
                columns = ",".join(quoted(c) for c in selector["columns"])
                self.database.execute(f"INSERT INTO {queue} SELECT {columns} FROM {quoted(selector['table'])} "
                                      f"WHERE TRUE ON CONFLICT({','.join(keys)}) DO NOTHING")
        return queue, keys

    def _lookup(self, item):
        queue, queue_columns = self._queue(item)
        state = self.progress(item["table"])
        columns = ",".join(queue_columns)
        parameters = () if state["cursor"] is None else tuple(state["cursor"])
        where = "" if not parameters else f"WHERE ({columns}) > ({','.join('?' for _ in parameters)})"
        keys = list(self.database.execute(f"SELECT {columns} FROM {queue} {where} ORDER BY {columns} LIMIT ?",
                                          (*parameters, self.lookup_keys)))
        if not keys:
            self._commit_batch(item, [], state["cursor"], complete=True)
            return
        native_keys = ",".join(quoted(k) for k in item["keys"])
        group = "(" + ",".join("?" for _ in item["keys"]) + ")"
        values = ",".join(group for _ in keys)
        selected_columns = ",".join(quoted(c) for c in item["columns"])
        sql = (f"SELECT rowid,{selected_columns} FROM {quoted(item['table'])} "
               f"WHERE ({native_keys}) IN (VALUES {values}) ORDER BY {native_keys},rowid")
        bindings = tuple(value for key in keys for value in key)
        query_plan = [r[-1] for r in self.source.execute("EXPLAIN QUERY PLAN " + sql, bindings)]
        if not any("SEARCH " in detail and "INDEX" in detail for detail in query_plan):
            raise ValueError("Native acquisition lookup must use an index: " + repr(query_plan))
        rows, size, found = [], 0, set()
        positions = [item["columns"].index(k) + 1 for k in item["keys"]]
        for row in self.source.execute(sql, bindings):
            row = tuple(row)
            size += len(canonical_row(row))
            if size > self.batch_bytes:
                raise LimitExceeded("Native lookup batch exceeds its byte limit; reduce lookup_keys")
            rows.append(row)
            found.add(tuple(row[i] for i in positions))
        missing = sorted(set(keys) - found)
        self._commit_batch(item, rows, keys[-1], selected=len(keys), missing=missing)

    def copy(self, table, *, max_batches=None, on_batch=None):
        """Copy to completion or a caller's explicit batch checkpoint."""
        item = self.tables[table]
        batches = 0
        while not self.progress(table)["complete"]:
            if max_batches is not None and batches >= max_batches:
                break
            self._lookup(item) if item.get("select_from") else self._scan(item)
            batches += 1
            if on_batch is not None:
                on_batch(table, self.progress(table))
        return self.progress(table)

    def summary(self):
        tables = {name: self.progress(name) for name in self.tables}
        logical = json.dumps({name: {k: row[k] for k in ("rows", "selected_keys", "missing_keys", "digest", "complete")}
                              for name, row in tables.items()}, sort_keys=True, separators=(",", ":"))
        return dict(selection_complete=all(row["complete"] for row in tables.values()), tables=tables,
                    logical_sha256=hashlib.sha256(logical.encode()).hexdigest(),
                    scope="native selected facts; normalized dataset acceptance remains separate")

    def close(self):
        if self.database is not None:
            self.database.close()
            self.database = None
