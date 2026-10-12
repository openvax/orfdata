"""Streaming, reference-partitioned quantitative datasets stored as SQLite."""

import hashlib
import json
import os
import sqlite3
from contextlib import closing
from pathlib import Path

from .measurements import (
    MEASUREMENT_STATUSES,
    QUANTITIES,
    Assay,
    DatasetManifest,
    Measurement,
    Sample,
    TargetMapping,
)
from .reference import ReferenceContext, ScopedTarget

_SCHEMA = """
CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE reference (key TEXT PRIMARY KEY, record TEXT NOT NULL);
CREATE TABLE target (
    reference_key TEXT NOT NULL REFERENCES reference(key),
    key TEXT NOT NULL, kind TEXT NOT NULL, namespace TEXT NOT NULL,
    identifier TEXT NOT NULL, PRIMARY KEY(reference_key, key));
CREATE TABLE sample (
    reference_key TEXT NOT NULL REFERENCES reference(key),
    id TEXT NOT NULL, record TEXT NOT NULL, PRIMARY KEY(reference_key, id));
CREATE TABLE assay (
    reference_key TEXT NOT NULL, id TEXT NOT NULL, sample_id TEXT NOT NULL,
    record TEXT NOT NULL, PRIMARY KEY(reference_key, id),
    UNIQUE(reference_key, id, sample_id),
    FOREIGN KEY(reference_key, sample_id) REFERENCES sample(reference_key, id));
CREATE TABLE mapping (
    reference_key TEXT NOT NULL REFERENCES reference(key),
    id TEXT NOT NULL, record TEXT NOT NULL, PRIMARY KEY(reference_key, id));
CREATE TABLE measurement (
    reference_key TEXT NOT NULL, id TEXT NOT NULL, target_key TEXT NOT NULL,
    sample_id TEXT NOT NULL, assay_id TEXT NOT NULL, quantity TEXT NOT NULL,
    status TEXT NOT NULL, record TEXT NOT NULL, PRIMARY KEY(reference_key, id),
    FOREIGN KEY(reference_key, target_key) REFERENCES target(reference_key, key),
    FOREIGN KEY(reference_key, assay_id, sample_id)
        REFERENCES assay(reference_key, id, sample_id));
CREATE TABLE input_measurement (
    reference_key TEXT NOT NULL, output_id TEXT NOT NULL, input_id TEXT NOT NULL,
    PRIMARY KEY(reference_key, output_id, input_id),
    FOREIGN KEY(reference_key, output_id) REFERENCES measurement(reference_key, id),
    FOREIGN KEY(reference_key, input_id) REFERENCES measurement(reference_key, id));
CREATE INDEX target_lookup ON target(reference_key, namespace, kind, identifier);
CREATE INDEX observation_lookup ON measurement(reference_key, target_key, sample_id, assay_id);
CREATE INDEX assay_lookup ON measurement(reference_key, assay_id, quantity, status);
"""


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _target_key(target):
    return hashlib.sha256(_json(target.key).encode()).hexdigest()


class QuantitativeDataset:
    """A read-only release. Constructing or querying it never downloads data.

    Local IDs are unique only within one dataset release AND reference context.
    Iterators preserve that context even when IDs collide between species.
    Registry-managed files are checked against pinned hashes before opening.
    """

    def __init__(self, path):
        self.path = Path(path).resolve()
        with closing(self._connect()) as connection:
            if connection.execute("PRAGMA user_version").fetchone()[0] != 1:
                raise ValueError("Unsupported quantitative database schema")
            row = connection.execute("SELECT value FROM metadata WHERE name='manifest'").fetchone()
            if row is None:
                raise ValueError("Missing quantitative manifest")
            self.manifest = DatasetManifest.from_dict(json.loads(row[0]))
        self._references = {r.key: r for r in self.manifest.references}

    def _connect(self):
        # mode=ro prevents SQLite's default creation of an empty missing database.
        connection = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True)
        connection.execute("PRAGMA query_only=ON")
        return connection

    @classmethod
    def create(cls, path, manifest, *, samples, assays, measurements, mappings=()):
        """Write a new dataset from iterators, validating scope and foreign keys.

        Never overwrite an existing file. An interrupted/invalid build removes
        only the new file it created. For atomic managed publication, call this
        inside DataCache's materialize builder. Derived observations must follow
        their input observations, preventing cyclic or dangling provenance.
        """
        if not isinstance(manifest, DatasetManifest):
            raise TypeError("manifest must be a DatasetManifest")
        path = Path(path)
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        os.close(descriptor)
        connection = None
        try:
            connection = sqlite3.connect(path)
            connection.execute("PRAGMA foreign_keys=ON")
            connection.executescript(_SCHEMA)
            connection.execute("PRAGMA user_version=1")
            connection.execute("INSERT INTO metadata VALUES ('manifest', ?)",
                               (_json(manifest.to_dict()),))
            for reference in manifest.references:
                connection.execute("INSERT INTO reference VALUES (?, ?)",
                                   (reference.key, _json(reference.to_dict())))
            for sample in samples:
                if not isinstance(sample, Sample):
                    raise TypeError("samples must contain Samples")
                connection.execute("INSERT INTO sample VALUES (?, ?, ?)",
                                   (sample.reference.key, sample.sample_id, _json(sample.to_dict())))
            for assay in assays:
                if not isinstance(assay, Assay):
                    raise TypeError("assays must contain Assays")
                if assay.kind not in manifest.assay_kinds:
                    raise ValueError("Assay kind absent from manifest coverage")
                connection.execute("INSERT INTO assay VALUES (?, ?, ?, ?)",
                                   (assay.reference.key, assay.assay_id, assay.sample_id,
                                    _json(assay.to_dict())))
            for mapping in mappings:
                if not isinstance(mapping, TargetMapping):
                    raise TypeError("mappings must contain TargetMappings")
                connection.execute("INSERT INTO mapping VALUES (?, ?, ?)",
                                   (mapping.reference.key, mapping.mapping_id,
                                    _json(mapping.to_dict())))
            source_ids = {s.asset_id for s in manifest.sources}
            for observation in measurements:
                if not isinstance(observation, Measurement):
                    raise TypeError("measurements must contain Measurements")
                if (observation.dataset_id, observation.dataset_version) != (manifest.dataset_id, manifest.version):
                    raise ValueError("Measurement belongs to a different dataset release")
                source = observation.provenance.get("source_asset")
                if source is not None and source not in source_ids:
                    raise ValueError("Measurement provenance names an unknown source asset")
                if source is None and not observation.input_measurements:
                    raise ValueError("Measurement requires a source_asset or input measurements")
                reference_key = observation.target.reference.key
                for identifier in observation.input_measurements:
                    if connection.execute(
                        "SELECT 1 FROM measurement WHERE reference_key=? AND id=?",
                        (reference_key, identifier),
                    ).fetchone() is None:
                        raise ValueError("Derived inputs must precede outputs in the same reference context")
                target = observation.target
                target_key = _target_key(target)
                connection.execute("INSERT OR IGNORE INTO target VALUES (?, ?, ?, ?, ?)",
                                   (reference_key, target_key, target.kind,
                                    target.namespace, target.identifier))
                record = observation.to_dict()
                record.pop("target")
                connection.execute("INSERT INTO measurement VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                                   (reference_key, observation.measurement_id, target_key,
                                    observation.sample_id, observation.assay_id,
                                    observation.quantity, observation.status, _json(record)))
                connection.executemany("INSERT INTO input_measurement VALUES (?, ?, ?)",
                                       ((reference_key, observation.measurement_id, identifier)
                                        for identifier in observation.input_measurements))
            counts = {table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
                      for table in ("sample", "assay", "mapping", "target", "measurement")}
            connection.execute("INSERT INTO metadata VALUES ('counts', ?)", (_json(counts),))
            connection.commit()
            connection.close()
            connection = None
            return cls(path)
        except BaseException:
            if connection is not None:
                connection.close()
            path.unlink()
            raise

    def _scope(self, reference):
        if reference is None:
            return None
        key = reference.key if isinstance(reference, ReferenceContext) else reference
        if key not in self._references:
            raise ValueError("Reference context is not present in this dataset")
        return key

    def _records(self, table, cls, *, reference=None):
        scope = self._scope(reference)
        query = f"SELECT record FROM {table}"
        params = ()
        if scope is not None:
            query += " WHERE reference_key=?"
            params = (scope,)
        query += " ORDER BY reference_key, id"
        with closing(self._connect()) as connection:
            for row in connection.execute(query, params):
                yield cls.from_dict(json.loads(row[0]))

    def samples(self, *, reference=None):
        return self._records("sample", Sample, reference=reference)

    def assays(self, *, reference=None):
        return self._records("assay", Assay, reference=reference)

    def mappings(self, *, reference=None):
        return self._records("mapping", TargetMapping, reference=reference)

    def measurements(self, *, target=None, reference=None, sample_id=None,
                     assay_id=None, quantity=None, status=None):
        """Iterate native observations; a target match is exact, never orthology.

        Unqualified IDs across multiple references are rejected. To browse all
        species together, omit ID filters; each result still carries its scope.
        """
        if quantity is not None and quantity not in QUANTITIES:
            raise ValueError("Unknown measurement quantity")
        if status is not None and status not in MEASUREMENT_STATUSES:
            raise ValueError("Unknown measurement status")
        scope = self._scope(reference)
        if target is not None:
            if not isinstance(target, ScopedTarget):
                raise TypeError("target must be a ScopedTarget")
            target_scope = self._scope(target.reference)
            if scope is not None and scope != target_scope:
                raise ValueError("Target and requested reference disagree")
            scope = target_scope
        if scope is None and len(self._references) > 1 and (sample_id or assay_id):
            raise ValueError("An ID query across multiple references requires an explicit reference")
        filters = dict(reference_key=scope, target_key=_target_key(target) if target else None,
                       sample_id=sample_id, assay_id=assay_id, quantity=quantity, status=status)
        selected = [(key, value) for key, value in filters.items() if value is not None]
        query = ("SELECT m.record, m.reference_key, t.kind, t.namespace, t.identifier "
                 "FROM measurement m JOIN target t ON m.reference_key=t.reference_key "
                 "AND m.target_key=t.key")
        if selected:
            query += " WHERE " + " AND ".join(f"m.{key}=?" for key, _ in selected)
        query += " ORDER BY m.reference_key, m.id"
        with closing(self._connect()) as connection:
            for row in connection.execute(query, [value for _, value in selected]):
                values = json.loads(row[0])
                values["target"] = ScopedTarget(self._references[row[1]], row[2], row[4], row[3])
                yield Measurement(**values)

    @property
    def counts(self):
        with closing(self._connect()) as connection:
            result = json.loads(connection.execute(
                "SELECT value FROM metadata WHERE name='counts'").fetchone()[0])
        return result

    def export_parquet(self, path, **query):
        """Export query results in bounded batches; requires orfdata[columnar]."""
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise ImportError("Install orfdata[columnar] to export Parquet") from exc
        schema = pa.schema([(name, pa.string()) for name in
                            ("dataset_id", "version", "reference_key", "measurement_id",
                             "sample_id", "assay_id", "target_kind", "target_namespace",
                             "target_identifier", "quantity", "unit", "status", "origin",
                             "scale", "measurement_json")] + [("value_float64", pa.float64())])
        with pq.ParquetWriter(path, schema) as writer:
            batch = []
            for value in self.measurements(**query):
                batch.append(dict(dataset_id=self.manifest.dataset_id, version=self.manifest.version,
                                  reference_key=value.target.reference.key,
                                  measurement_id=value.measurement_id,
                                  sample_id=value.sample_id, assay_id=value.assay_id,
                                  target_kind=value.target.kind,
                                  target_namespace=value.target.namespace,
                                  target_identifier=value.target.identifier, quantity=value.quantity,
                                  unit=value.unit, status=value.status, origin=value.origin,
                                  scale=value.scale,
                                  value_float64=float(value.value) if value.value is not None else None,
                                  measurement_json=_json(value.to_dict())))
                if len(batch) == 4096:
                    writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                    batch = []
            if batch:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
