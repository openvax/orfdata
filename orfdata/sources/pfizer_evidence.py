"""Compact source records and normalized numeric vectors, with no acquisition.

This is the evidence component of the canine dataset, not a vendor SQLite file.
Every selected scalar, opaque byte field and numeric vector can be reconstructed
exactly. Native IDs and confidence/modification definitions remain separate from
peptide/protein attribution and the common measurement projection.
"""

import json
import math
import re
import sqlite3
import struct
from contextlib import closing
from pathlib import Path

from ..numeric_vector import NumericVector

ENCODING = "orfdata.pfizer-evidence.v2"
SUPPORTED_ENCODINGS = ("orfdata.pfizer-evidence.v1", ENCODING)
SCHEMA = """
CREATE TABLE evidence_metadata(name TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE evidence_table(i INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL,
    role TEXT NOT NULL, columns_json TEXT NOT NULL);
CREATE TABLE evidence_record(i INTEGER PRIMARY KEY, table_i INTEGER NOT NULL REFERENCES evidence_table(i),
    source_rowid INTEGER NOT NULL, source_key TEXT NOT NULL, sequence TEXT, attributes TEXT NOT NULL,
    UNIQUE(table_i,source_rowid));
CREATE INDEX evidence_record_key ON evidence_record(table_i,source_key);
CREATE INDEX evidence_record_sequence ON evidence_record(table_i,sequence);
CREATE TABLE evidence_binary(record_i INTEGER NOT NULL REFERENCES evidence_record(i), column_name TEXT NOT NULL,
    value_kind TEXT NOT NULL CHECK(value_kind IN ('bytes','float64')), data BLOB NOT NULL,
    PRIMARY KEY(record_i,column_name)) WITHOUT ROWID;
CREATE TABLE evidence_vector(record_i INTEGER NOT NULL REFERENCES evidence_record(i), column_name TEXT NOT NULL,
    value_type TEXT NOT NULL, element_count INTEGER NOT NULL, values_data BLOB NOT NULL,
    presence_data BLOB NOT NULL, source_null INTEGER NOT NULL CHECK(source_null IN (0,1)),
    PRIMARY KEY(record_i,column_name)) WITHOUT ROWID;
"""


def record_keys(table):
    """Source-declared entity IDs; link rows retain their separate source rowid."""
    if table in ("TargetProteins", "DecoyProteins"):
        return ("UniqueSequenceID",)
    if table in ("TargetProteinGroups", "DecoyProteinGroups"):
        return ("ProteinGroupID",)
    if table in ("TargetPeptideGroups", "DecoyPeptideGroups"):
        return ("PeptideGroupID",)
    if table in ("TargetPsms", "DecoyPsms"):
        return ("WorkflowID", "PeptideID")
    if table in ("MSnSpectrumInfo", "QuanSpectrumInfo"):
        return ("WorkflowID", "SpectrumID")
    if table == "ReporterQuanResults":
        return ("WorkflowID", "QuanResultID", "QuanChannelID")
    if table == "FoundModifications":
        return ("ModificationID",)
    if table == "ProteinIdentificationGroups":
        return ("WorkflowID", "ProteinIdentificationGroupID")
    return ("_source_rowid",)


def source_key(values):
    # Decimal strings with explicit types prevent int64 IDs from losing precision
    # when an exported key is read by a runtime whose JSON numbers are doubles.
    if any(type(v) is not int for v in values):
        raise TypeError("Native source key components must be integers")
    return json.dumps([str(v) for v in values], separators=(",", ":"))


def encode_scalars(attributes, columns):
    """Store scalar presence and values in the declared source-column order.

    A present null is different from an omitted binary/vector field. The mask
    is hexadecimal text, so it never depends on JSON integer precision.
    """
    if not set(attributes) <= set(columns):
        raise ValueError("Scalar fields differ from the declared source columns")
    mask = sum(1 << index for index, name in enumerate(columns) if name in attributes)
    values = [attributes[name] for name in columns if name in attributes]
    return json.dumps([format(mask, "x"), values], separators=(",", ":"), allow_nan=False)


def decode_scalars(payload, columns):
    """Read v2 ordered scalars and the original v1 JSON-object representation."""
    encoded = json.loads(payload)
    if isinstance(encoded, dict):
        if not set(encoded) <= set(columns):
            raise ValueError("Scalar fields differ from the declared source columns")
        result = encoded
    else:
        if (not isinstance(encoded, list) or len(encoded) != 2 or not isinstance(encoded[0], str)
                or re.fullmatch(r"0|[1-9a-f][0-9a-f]*", encoded[0]) is None
                or not isinstance(encoded[1], list)):
            raise ValueError("Invalid compact scalar presence/value encoding")
        mask, values = int(encoded[0], 16), encoded[1]
        if mask.bit_length() > len(columns) or mask.bit_count() != len(values):
            raise ValueError("Compact scalar presence/value count differs from declared columns")
        names = (name for index, name in enumerate(columns) if mask & (1 << index))
        result = dict(zip(names, values, strict=True))
    if any(type(value) not in (str, int, float, bool, type(None))
           or isinstance(value, float) and not math.isfinite(value) for value in result.values()):
        raise ValueError("Scalar payload contains a non-scalar or nonfinite value")
    return result


def column_layouts(connection):
    return {name: tuple(json.loads(columns)) for name, columns in connection.execute(
        "SELECT name,columns_json FROM evidence_table")}


class EvidenceWriter:
    """Insert into a caller-owned SQLite transaction and explicit disk budget.

    It does not commit, choose a source selection, certify completion or open
    files. The extraction/build command owns those lifecycle and quota checks.
    """

    def __init__(self, connection, plan, decoder):
        self.connection, self.decoder = connection, decoder
        connection.executescript(SCHEMA)
        self.tables = {}
        for index, item in enumerate(plan["tables"], 1):
            self.tables[item["table"]] = index, tuple(item["columns"])
            connection.execute("INSERT INTO evidence_table VALUES (?,?,?,?)",
                               (index, item["table"], item["role"], json.dumps(item["columns"])))
        connection.execute("INSERT INTO evidence_metadata VALUES ('encoding',?)", (json.dumps(ENCODING),))
        connection.execute("INSERT INTO evidence_metadata VALUES ('state','\"building\"')")

    def insert(self, table, row):
        index, columns = self.tables[table]
        if set(row) != {*columns, "_source_rowid"}:
            raise ValueError("Evidence record columns differ from the declared source selection")
        attributes, binaries, vectors = {}, [], []
        for name in columns:
            value = row[name]
            data_type = self.decoder.tables.get(table)
            definition = self.decoder.columns.get((data_type["DataTypeID"], name)) if data_type else None
            if definition is not None and definition["SpecialValueType"] == "DataDistribution":
                _, distribution, _ = self.decoder.distribution(table, name)
                kind = self.decoder.numeric_types[distribution["CustomDataType"]]
                if kind in ("Int", "Int64", "Double"):
                    vector = self.decoder.vector(table, name, value)
                    vectors.append((name, vector.value_type, vector.count, vector.values, vector.present, int(vector.source_null)))
                    continue
            if isinstance(value, bytes):
                binaries.append((name, "bytes", value))
            elif isinstance(value, float) and not math.isfinite(value):
                binaries.append((name, "float64", struct.pack("<d", value)))
            else:
                attributes[name] = value
        key = source_key(tuple(row[k] for k in record_keys(table)))
        cursor = self.connection.execute(
            "INSERT INTO evidence_record(table_i,source_rowid,source_key,sequence,attributes) VALUES (?,?,?,?,?)",
            (index, row["_source_rowid"], key, row.get("Sequence"),
             encode_scalars(attributes, columns)))
        identifier = cursor.lastrowid
        self.connection.executemany("INSERT INTO evidence_binary VALUES (?,?,?,?)",
                                    ((identifier, *record) for record in binaries))
        self.connection.executemany("INSERT INTO evidence_vector VALUES (?,?,?,?,?,?,?)",
                                    ((identifier, *record) for record in vectors))
        return identifier


class EvidenceStore:
    """Read a completed local evidence component without downloads or vendor tools."""

    def __init__(self, path):
        self.path = Path(path).resolve()
        with closing(self.connect()) as db:
            self.metadata = {key: json.loads(value) for key, value in db.execute("SELECT name,value FROM evidence_metadata")}
            if self.metadata.get("encoding") not in SUPPORTED_ENCODINGS or self.metadata.get("state") != "complete":
                raise ValueError("Evidence component is incomplete or uses an unsupported encoding")
            self.tables = {name: (index, tuple(json.loads(columns))) for index, name, columns in db.execute(
                "SELECT i,name,columns_json FROM evidence_table ORDER BY i")}

    def connect(self):
        db = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True)
        db.execute("PRAGMA query_only=ON")
        db.execute("PRAGMA cache_size=-8192")
        return db

    @staticmethod
    def _vector(row):
        kind, count, values, present, source_null = row
        return NumericVector(kind, count, values, present, bool(source_null))

    @classmethod
    def restore(cls, db, record, *, columns=None):
        """Restore exact native fields; numeric vectors retain original source bytes."""
        from .pfizer_vendor import DistributionDecoder

        identifier, source_rowid, attributes = record
        if columns is None:
            layout = db.execute("SELECT t.columns_json FROM evidence_record r "
                                "JOIN evidence_table t ON t.i=r.table_i WHERE r.i=?", (identifier,)).fetchone()
            if layout is None:
                raise ValueError("Stored evidence record has no declared column layout")
            columns = tuple(json.loads(layout[0]))
        result = decode_scalars(attributes, columns)
        result["_source_rowid"] = source_rowid
        for column, kind, data in db.execute("SELECT column_name,value_kind,data FROM evidence_binary WHERE record_i=?",
                                             (identifier,)):
            if column in result:
                raise ValueError("Stored evidence field has conflicting scalar/binary/vector payloads")
            if kind == "float64":
                if len(data) != 8:
                    raise ValueError("Invalid nonfinite scalar encoding")
                result[column] = struct.unpack("<d", data)[0]
            else:
                result[column] = data
        for row in db.execute("SELECT column_name,value_type,element_count,values_data,presence_data,source_null "
                              "FROM evidence_vector WHERE record_i=?", (identifier,)):
            if row[0] in result:
                raise ValueError("Stored evidence field has conflicting scalar/binary/vector payloads")
            result[row[0]] = DistributionDecoder.native_bytes(cls._vector(row[1:]))
        if set(result) != {*columns, "_source_rowid"}:
            raise ValueError("Stored evidence fields differ from the declared source selection")
        return result

    def records(self, table, *, key=None, sequence=None):
        """Stream complete native records in source row order, with exact ID/sequence filters."""
        table_id, columns = self.tables[table]
        where, bindings = ["table_i=?"], [table_id]
        if key is not None:
            where.append("source_key=?")
            bindings.append(source_key(tuple(key)))
        if sequence is not None:
            where.append("sequence=?")
            bindings.append(sequence)
        sql = "SELECT i,source_rowid,attributes FROM evidence_record WHERE " + " AND ".join(where) + " ORDER BY source_rowid"
        with closing(self.connect()) as db:
            for record in db.execute(sql, bindings):
                yield self.restore(db, record, columns=columns)

    def vector(self, table, key, column):
        """Get an independently typed numeric vector for one original record."""
        with closing(self.connect()) as db:
            rows = db.execute(
                "SELECT v.value_type,v.element_count,v.values_data,v.presence_data,v.source_null "
                "FROM evidence_record r JOIN evidence_vector v ON r.i=v.record_i "
                "WHERE r.table_i=? AND r.source_key=? AND v.column_name=?",
                (self.tables[table][0], source_key(tuple(key)), column)).fetchall()
            if len(rows) != 1:
                raise ValueError("Expected one native numeric vector for the requested record/column")
            return self._vector(rows[0])
