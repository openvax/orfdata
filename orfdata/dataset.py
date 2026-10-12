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
    MeasurementReference,
    Sample,
    TargetMapping,
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
    reference_key TEXT NOT NULL, output_id TEXT NOT NULL,
    input_reference_key TEXT NOT NULL, input_id TEXT NOT NULL,
    PRIMARY KEY(reference_key, output_id, input_reference_key, input_id),
    FOREIGN KEY(reference_key, output_id) REFERENCES measurement(reference_key, id),
    FOREIGN KEY(input_reference_key, input_id) REFERENCES measurement(reference_key, id));
CREATE INDEX target_lookup ON target(reference_key, namespace, kind, identifier);
CREATE INDEX observation_lookup ON measurement(reference_key, target_key, sample_id, assay_id);
CREATE INDEX assay_lookup ON measurement(reference_key, assay_id, quantity, status);
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
            if self._storage_version not in (1, 2, 3):
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
            connection.execute("PRAGMA user_version=3")
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
                    connection.execute("INSERT OR IGNORE INTO target VALUES (?, ?, ?, ?, ?)",
                                       (target.reference.key, _target_key(target), target.kind,
                                        target.namespace, target.identifier))
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
                    connection.execute("INSERT OR IGNORE INTO target VALUES (?, ?, ?, ?, ?)",
                                       (target.reference.key, _target_key(target), target.kind,
                                        target.namespace, target.identifier))
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
                reference_key = observation.target.reference.key
                inputs = [(reference_key, identifier) for identifier in observation.input_measurements]
                inputs.extend((value.reference.key, value.measurement_id) for value in observation.input_references)
                for input_reference, identifier in inputs:
                    if connection.execute(
                        "SELECT 1 FROM measurement WHERE reference_key=? AND id=?",
                        (input_reference, identifier),
                    ).fetchone() is None:
                        raise ValueError("Derived inputs must precede outputs in their declared reference contexts")
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
                connection.executemany("INSERT INTO input_measurement VALUES (?, ?, ?, ?)",
                                       ((reference_key, observation.measurement_id, input_reference, identifier)
                                        for input_reference, identifier in inputs))
            counts = {table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
                      for table in ("sample", "assay", "mapping", "target", "measurement",
                                    "assay_pair", "target_relation", "reference_mapping")}
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
