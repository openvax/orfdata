"""Streaming, reference-partitioned quantitative datasets stored as SQLite."""

import hashlib
import json
import os
import sqlite3
from contextlib import closing
from dataclasses import MISSING, fields
from functools import lru_cache
from pathlib import Path

from .measurements import (
    MEASUREMENT_STATUSES,
    QUANTITIES,
    Assay,
    DatasetManifest,
    Measurement,
    MeasurementReference,
    Sample,
    TargetMapping,
    _plain,
)
from .pairing import (
    PAIRING_BASES,
    TEMPORAL_RELATIONS,
    AssayPair,
    MatchedMeasurements,
    TargetRelation,
    _check_one_to_one,
)
from .reference import ReferenceContext, ScopedTarget
from .reference_mapping import MappedMeasurement, ReferenceMapping

_SCHEMA = """
CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE reference (i INTEGER PRIMARY KEY, key TEXT UNIQUE NOT NULL,
    record TEXT NOT NULL, UNIQUE(key, i));
CREATE TABLE target (
    i INTEGER PRIMARY KEY, reference_i INTEGER NOT NULL,
    reference_key TEXT NOT NULL,
    key TEXT NOT NULL, kind TEXT NOT NULL, namespace TEXT NOT NULL,
    identifier TEXT NOT NULL, UNIQUE(reference_key, key), UNIQUE(reference_i, i),
    FOREIGN KEY(reference_key, reference_i) REFERENCES reference(key, i));
CREATE TABLE sample (
    reference_key TEXT NOT NULL REFERENCES reference(key),
    id TEXT NOT NULL, record TEXT NOT NULL, PRIMARY KEY(reference_key, id));
CREATE TABLE assay (
    i INTEGER PRIMARY KEY, reference_i INTEGER NOT NULL,
    reference_key TEXT NOT NULL, id TEXT NOT NULL, sample_id TEXT NOT NULL,
    record TEXT NOT NULL, UNIQUE(reference_key, id), UNIQUE(reference_i, i),
    UNIQUE(reference_key, id, sample_id),
    FOREIGN KEY(reference_key, reference_i) REFERENCES reference(key, i),
    FOREIGN KEY(reference_key, sample_id) REFERENCES sample(reference_key, id));
CREATE TABLE mapping (
    reference_key TEXT NOT NULL REFERENCES reference(key),
    id TEXT NOT NULL, record TEXT NOT NULL, PRIMARY KEY(reference_key, id));
CREATE TABLE observation (
    i INTEGER PRIMARY KEY, reference_i INTEGER NOT NULL,
    id TEXT NOT NULL, target_i INTEGER NOT NULL, assay_i INTEGER NOT NULL,
    quantity TEXT NOT NULL, status TEXT NOT NULL, record TEXT NOT NULL,
    UNIQUE(reference_i, id),
    FOREIGN KEY(reference_i, target_i) REFERENCES target(reference_i, i),
    FOREIGN KEY(reference_i, assay_i) REFERENCES assay(reference_i, i));
CREATE TABLE observation_input (
    output_i INTEGER NOT NULL REFERENCES observation(i),
    input_i INTEGER NOT NULL REFERENCES observation(i), PRIMARY KEY(output_i, input_i));
CREATE VIEW measurement AS
    SELECT r.key reference_key, m.id, t.key target_key, a.sample_id,
           a.id assay_id, m.quantity, m.status, m.record
    FROM observation m JOIN reference r ON r.i=m.reference_i
    JOIN target t ON t.i=m.target_i JOIN assay a ON a.i=m.assay_i;
CREATE VIEW input_measurement AS
    SELECT r.key reference_key, o.id output_id, s.key input_reference_key, i.id input_id
    FROM observation_input p JOIN observation o ON o.i=p.output_i
    JOIN observation i ON i.i=p.input_i JOIN reference r ON r.i=o.reference_i
    JOIN reference s ON s.i=i.reference_i;
CREATE INDEX target_lookup ON target(reference_key, namespace, kind, identifier);
CREATE INDEX observation_lookup ON observation(reference_i, target_i, assay_i);
CREATE INDEX assay_lookup ON observation(reference_i, assay_i, quantity, status);
CREATE TABLE assay_pair (
    reference_key TEXT NOT NULL REFERENCES reference(key), id TEXT NOT NULL,
    right_reference_key TEXT NOT NULL REFERENCES reference(key),
    left_assay_id TEXT, right_assay_id TEXT, record TEXT NOT NULL,
    PRIMARY KEY(reference_key, id),
    FOREIGN KEY(reference_key, left_assay_id) REFERENCES assay(reference_key, id),
    FOREIGN KEY(right_reference_key, right_assay_id) REFERENCES assay(reference_key, id));
CREATE TABLE target_relation (
    reference_key TEXT NOT NULL, id TEXT NOT NULL,
    right_reference_key TEXT NOT NULL,
    left_target_key TEXT NOT NULL, right_target_key TEXT NOT NULL, record TEXT NOT NULL,
    PRIMARY KEY(reference_key, id),
    FOREIGN KEY(reference_key, left_target_key) REFERENCES target(reference_key, key),
    FOREIGN KEY(right_reference_key, right_target_key) REFERENCES target(reference_key, key));
CREATE TABLE reference_mapping (
    reference_key TEXT NOT NULL, id TEXT NOT NULL, source_target_key TEXT NOT NULL,
    target_reference_key TEXT NOT NULL REFERENCES reference(key), record TEXT NOT NULL,
    PRIMARY KEY(reference_key, id),
    UNIQUE(reference_key, source_target_key, target_reference_key),
    FOREIGN KEY(reference_key, source_target_key) REFERENCES target(reference_key, key));
CREATE TABLE reference_assignment (
    source_reference_key TEXT NOT NULL, mapping_id TEXT NOT NULL,
    target_reference_key TEXT NOT NULL, target_key TEXT NOT NULL,
    PRIMARY KEY(source_reference_key, mapping_id, target_reference_key, target_key),
    FOREIGN KEY(source_reference_key, mapping_id) REFERENCES reference_mapping(reference_key, id),
    FOREIGN KEY(target_reference_key, target_key) REFERENCES target(reference_key, key));
CREATE INDEX reference_assignment_lookup ON reference_assignment(target_reference_key, target_key);
"""


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _target_key(target):
    return hashlib.sha256(_json(target.key).encode()).hexdigest()


_STORED_FIELDS = tuple(f for f in fields(Measurement) if f.name not in {
    "measurement_id", "target", "sample_id", "assay_id", "quantity", "status",
    "dataset_id", "dataset_version"})
_V4_DEFAULTS = {
    "origin": "observed", "scale": "linear", "log_base": None,
    "uncertainty": None, "uncertainty_type": None, "detection_limit": None,
    "input_measurements": (), "input_references": (), "derivation": {}, "provenance": {},
}


def _measurement_record(observation):
    # Schema defaults are fixed here, independent of future model defaults.
    # Do not serialize the target
    # and then discard it: its full reference is already stored once.
    result = {}
    for field in _STORED_FIELDS:
        value = getattr(observation, field.name)
        default = _V4_DEFAULTS.get(field.name, MISSING)
        if default is MISSING or value != default:
            result[field.name] = _plain(value)
    return result


class QuantitativeDataset:
    """A read-only release. Constructing or querying it never downloads data.

    Local IDs are unique only within one dataset release AND reference context.
    Iterators preserve that context even when IDs collide between species.
    Registry-managed files are checked against pinned hashes before opening.
    """

    def __init__(self, path):
        self.path = Path(path).resolve()
        with closing(self._connect()) as connection:
            self._storage_version = connection.execute("PRAGMA user_version").fetchone()[0]
            if self._storage_version not in (1, 2, 3, 4):
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
    def create(cls, path, manifest, *, samples, assays, measurements, mappings=(),
               pairings=(), relations=(), reference_mappings=()):
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
            connection.execute("PRAGMA user_version=4")
            connection.execute("INSERT INTO metadata VALUES ('manifest', ?)",
                               (_json(manifest.to_dict()),))
            reference_ids = {}
            for reference in manifest.references:
                cursor = connection.execute("INSERT INTO reference (key, record) VALUES (?, ?)",
                                            (reference.key, _json(reference.to_dict())))
                reference_ids[reference.key] = cursor.lastrowid

            @lru_cache(maxsize=131072)
            def insert_target(target):
                reference_key = target.reference.key
                if reference_key not in reference_ids:
                    raise sqlite3.IntegrityError("Target reference absent from dataset")
                reference_i = reference_ids[reference_key]
                key = _target_key(target)
                connection.execute("INSERT OR IGNORE INTO target "
                                   "(reference_i, reference_key, key, kind, namespace, identifier) "
                                   "VALUES (?, ?, ?, ?, ?, ?)",
                                   (reference_i, reference_key, key, target.kind, target.namespace, target.identifier))
                target_i = connection.execute("SELECT i FROM target WHERE reference_key=? AND key=?",
                                              (reference_key, key)).fetchone()[0]
                return reference_i, target_i, key
            for sample in samples:
                if not isinstance(sample, Sample):
                    raise TypeError("samples must contain Samples")
                connection.execute("INSERT INTO sample VALUES (?, ?, ?)",
                                   (sample.reference.key, sample.sample_id, _json(sample.to_dict())))
            assay_ids = {}
            for assay in assays:
                if not isinstance(assay, Assay):
                    raise TypeError("assays must contain Assays")
                if assay.kind not in manifest.assay_kinds:
                    raise ValueError("Assay kind absent from manifest coverage")
                reference_i = reference_ids.get(assay.reference.key)
                cursor = connection.execute("INSERT INTO assay "
                                            "(reference_i, reference_key, id, sample_id, record) VALUES (?, ?, ?, ?, ?)",
                                            (reference_i, assay.reference.key, assay.assay_id, assay.sample_id,
                                             _json(assay.to_dict())))
                assay_ids[(reference_i, assay.assay_id, assay.sample_id)] = cursor.lastrowid
            for pair in pairings:
                if not isinstance(pair, AssayPair):
                    raise TypeError("pairings must contain AssayPairs")
                if (pair.dataset_id, pair.dataset_version) != (manifest.dataset_id, manifest.version):
                    raise ValueError("Pairing belongs to a different dataset release")
                if pair.basis in ("same_material", "adjacent_donor_tissue"):
                    units = []
                    for reference, identifier in ((pair.reference, pair.left_assay_id),
                                                  (pair.right_reference, pair.right_assay_id)):
                        row = connection.execute(
                            "SELECT s.record FROM assay a JOIN sample s ON a.reference_key=s.reference_key "
                            "AND a.sample_id=s.id WHERE a.reference_key=? AND a.id=?",
                            (reference.key, identifier)).fetchone()
                        if row:
                            unit = json.loads(row[0]).get("biological_unit_id")
                            if unit is not None:
                                units.append(unit)
                    if len(set(units)) > 1:
                        raise ValueError("Declared pairing contradicts biological-unit identities")
                connection.execute("INSERT INTO assay_pair VALUES (?, ?, ?, ?, ?, ?)",
                                   (pair.reference.key, pair.pair_id, pair.right_reference.key, pair.left_assay_id,
                                    pair.right_assay_id, _json(pair.to_dict())))
            relation_links = ({}, {})
            for relation in relations:
                if not isinstance(relation, TargetRelation):
                    raise TypeError("relations must contain TargetRelations")
                if (relation.dataset_id, relation.dataset_version) != (manifest.dataset_id, manifest.version):
                    raise ValueError("Target relation belongs to a different dataset release")
                for target in (relation.left, relation.right):
                    insert_target(target)
                for index, (first, second) in enumerate(((relation.left, relation.right), (relation.right, relation.left))):
                    _check_one_to_one(relation_links[index],
                                      (relation.relation, first.key, second.reference.key, second.kind, second.namespace), second.key,
                                      relation.cardinality == "one_to_one")
                connection.execute("INSERT INTO target_relation VALUES (?, ?, ?, ?, ?, ?)",
                                   (relation.left.reference.key, relation.relation_id,
                                    relation.right.reference.key,
                                    _target_key(relation.left), _target_key(relation.right),
                                    _json(relation.to_dict())))
            for mapping in mappings:
                if not isinstance(mapping, TargetMapping):
                    raise TypeError("mappings must contain TargetMappings")
                connection.execute("INSERT INTO mapping VALUES (?, ?, ?)",
                                   (mapping.reference.key, mapping.mapping_id,
                                    _json(mapping.to_dict())))
            for mapping in reference_mappings:
                if not isinstance(mapping, ReferenceMapping):
                    raise TypeError("reference_mappings must contain ReferenceMappings")
                if (mapping.dataset_id, mapping.dataset_version) != (manifest.dataset_id, manifest.version):
                    raise ValueError("Reference mapping belongs to a different dataset release")
                for target in (mapping.source, *mapping.targets):
                    insert_target(target)
                connection.execute("INSERT INTO reference_mapping VALUES (?, ?, ?, ?, ?)",
                                   (mapping.source.reference.key, mapping.mapping_id, _target_key(mapping.source),
                                    mapping.target_reference.key, _json(mapping.to_dict())))
                connection.executemany("INSERT INTO reference_assignment VALUES (?, ?, ?, ?)",
                                       ((mapping.source.reference.key, mapping.mapping_id,
                                         target.reference.key, _target_key(target)) for target in mapping.targets))
            source_ids = {s.asset_id for s in manifest.sources}
            for observation in measurements:
                if not isinstance(observation, Measurement):
                    raise TypeError("measurements must contain Measurements")
                if (observation.dataset_id, observation.dataset_version) != (manifest.dataset_id, manifest.version):
                    raise ValueError("Measurement belongs to a different dataset release")
                source = observation.provenance.get("source_asset")
                if source is not None and source not in source_ids:
                    raise ValueError("Measurement provenance names an unknown source asset")
                if source is None and not (observation.input_measurements or observation.input_references):
                    raise ValueError("Measurement requires a source_asset or input measurements")
                reference_i, target_i, _ = insert_target(observation.target)
                inputs = [(reference_i, identifier) for identifier in observation.input_measurements]
                inputs.extend((reference_ids.get(value.reference.key), value.measurement_id)
                              for value in observation.input_references)
                input_ids = []
                for input_reference, identifier in inputs:
                    found = connection.execute(
                        "SELECT i FROM observation WHERE reference_i=? AND id=?",
                        (input_reference, identifier),
                    ).fetchone()
                    if found is None:
                        raise ValueError("Derived inputs must precede outputs in their declared reference contexts")
                    input_ids.append(found[0])
                assay_i = assay_ids.get((reference_i, observation.assay_id, observation.sample_id))
                if assay_i is None:
                    raise sqlite3.IntegrityError("Observation assay/sample absent from its reference context")
                cursor = connection.execute("INSERT INTO observation "
                                            "(reference_i, id, target_i, assay_i, quantity, status, record) "
                                            "VALUES (?, ?, ?, ?, ?, ?, ?)",
                                            (reference_i, observation.measurement_id, target_i, assay_i,
                                             observation.quantity, observation.status,
                                             _json(_measurement_record(observation))))
                connection.executemany("INSERT INTO observation_input VALUES (?, ?)",
                                       ((cursor.lastrowid, input_i) for input_i in input_ids))
            insert_target.cache_clear()
            counts = {table: connection.execute(
                f"SELECT count(*) FROM {'observation' if table == 'measurement' else table}").fetchone()[0]
                      for table in ("sample", "assay", "mapping", "target", "measurement",
                                    "assay_pair", "target_relation", "reference_mapping")}
            connection.execute("INSERT INTO metadata VALUES ('counts', ?)", (_json(counts),))
            # Real cardinalities prevent the planner from scanning a whole
            # reference in ID order for a narrow target query.
            connection.execute("ANALYZE")
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

    def reference_mappings(self, *, source_reference=None, target=None):
        """Iterate original source-to-reference mappings, including unresolved IDs."""
        scope = self._scope(source_reference)
        if target is not None:
            if not isinstance(target, ScopedTarget):
                raise TypeError("target must be a ScopedTarget")
            self._scope(target.reference)
        if self._storage_version < 3:
            return
        query = "SELECT m.record FROM reference_mapping m"
        conditions, params = [], []
        if target is not None:
            query += (" JOIN reference_assignment a ON m.reference_key=a.source_reference_key "
                      "AND m.id=a.mapping_id")
            conditions.extend(("a.target_reference_key=?", "a.target_key=?"))
            params.extend((target.reference.key, _target_key(target)))
        if scope is not None:
            conditions.append("m.reference_key=?")
            params.append(scope)
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        query += " ORDER BY m.reference_key, m.id"
        with closing(self._connect()) as connection:
            for row in connection.execute(query, params):
                yield ReferenceMapping.from_dict(json.loads(row[0]))

    def mapped_measurements(self, target, *, mapping_statuses=("unique",),
                            require_sequence=False, **observation_filters):
        """Return original values plus mappings, never relabeled measurements.

        Ambiguous or unresolved candidate membership requires explicit status
        opt-in. Gene RNA cannot map to ORF abundance. require_sequence excludes
        identifier-only mappings; ORF assignments always require their own
        sequence/frame/position validation when the mapping is constructed.
        """
        allowed = set(mapping_statuses)
        if not allowed or allowed - {"unique", "ambiguous", "unresolved"}:
            raise ValueError("Choose explicit usable mapping statuses")
        if type(require_sequence) is not bool:
            raise TypeError("require_sequence must be a boolean")
        if set(observation_filters) - {"sample_id", "assay_id", "quantity", "status"}:
            raise ValueError("Mapped observation filters cannot override source/target reference identity")
        for mapping in self.reference_mappings(target=target):
            if mapping.status not in allowed or (require_sequence and not mapping.sequence_validated):
                continue
            for observation in self.measurements(target=mapping.source, **observation_filters):
                yield MappedMeasurement(observation, mapping, target)

    def pairings(self, *, reference=None):
        """Iterate source-declared pairings; legacy releases have none."""
        if self._storage_version == 1:
            self._scope(reference)
            return iter(())
        return self._records("assay_pair", AssayPair, reference=reference)

    def relations(self, *, reference=None):
        """Iterate same-species relationships, optionally scoped by their left reference."""
        if self._storage_version == 1:
            self._scope(reference)
            return iter(())
        return self._records("target_relation", TargetRelation, reference=reference)

    def matched_measurements(self, *, allowed_bases, allowed_temporal=("concurrent",),
                             reference=None, pair_ids=None, relation_ids=None):
        """Iterate explicit assay/target matches, preserving all values and missing mates.

        Defaults require concurrent assays; opting into longitudinal or unknown
        timing is explicit. No matching by tissue names or inferred orthology.
        Each side contains its own original observations, never a Cartesian
        multiplication or an aggregate over technical replicates/quantities.
        """
        bases, temporal = set(allowed_bases), set(allowed_temporal)
        if isinstance(pair_ids, str) or isinstance(relation_ids, str):
            raise TypeError("Relationship ID selections must be sequences, not a single string")
        if not bases or bases - set(PAIRING_BASES):
            raise ValueError("Specify valid allowed pairing bases explicitly")
        if not temporal or temporal - set(TEMPORAL_RELATIONS):
            raise ValueError("Specify valid allowed temporal relationships")
        if (pair_ids or relation_ids) and reference is None and len(self._references) > 1:
            raise ValueError("Relationship ID queries across species require a reference")
        relations = tuple(relation for relation in self.relations(reference=reference)
                          if relation_ids is None or relation.relation_id in relation_ids)
        for pair in self.pairings(reference=reference):
            if pair.basis not in bases or pair.temporal_relation not in temporal:
                continue
            if pair_ids is not None and pair.pair_id not in pair_ids:
                continue
            by_target = []
            for assay_reference, assay in ((pair.reference, pair.left_assay_id),
                                           (pair.right_reference, pair.right_assay_id)):
                values = {}
                if assay is not None:
                    for value in self.measurements(reference=assay_reference, assay_id=assay):
                        values.setdefault(value.target.key, []).append(value)
                by_target.append(values)
            for relation in relations:
                if (relation.left.reference.key == pair.reference.key
                        and relation.right.reference.key == pair.right_reference.key):
                    yield MatchedMeasurements(pair, relation,
                                              tuple(by_target[0].get(relation.left.key, ())),
                                              tuple(by_target[1].get(relation.right.key, ())))

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
        query = ("SELECT m.record, m.reference_key, t.kind, t.namespace, t.identifier, "
                 "m.id, m.sample_id, m.assay_id, m.quantity, m.status "
                 "FROM measurement m JOIN target t ON m.reference_key=t.reference_key "
                 "AND m.target_key=t.key")
        if selected:
            query += " WHERE " + " AND ".join(f"m.{key}=?" for key, _ in selected)
        query += " ORDER BY m.reference_key, m.id"
        with closing(self._connect()) as connection:
            for row in connection.execute(query, [value for _, value in selected]):
                values = json.loads(row[0])
                if self._storage_version == 4:
                    values = dict(_V4_DEFAULTS, **values)
                    values.update(measurement_id=row[5], sample_id=row[6], assay_id=row[7],
                                  quantity=row[8], status=row[9], dataset_id=self.manifest.dataset_id,
                                  dataset_version=self.manifest.version)
                values["target"] = ScopedTarget(self._references[row[1]], row[2], row[4], row[3])
                if "input_references" in values:
                    values["input_references"] = tuple(MeasurementReference.from_dict(value)
                                                       for value in values["input_references"])
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
