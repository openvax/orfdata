"""Indexed common measurements projected from lossless canine evidence vectors.

Storage version 6 keeps numeric arrays and shared sample axes once. It has no
flat observation table: the public dataset API expands requested cells, including
missing cells, and reports the complete logical observation count. All source
facts remain accessible through EvidenceStore; no vendor file or network access
is needed to open or query this representation.
"""

import hashlib
import json
from contextlib import closing

from ..dataset import _json, _target_key
from ..reference import ScopedTarget
from .pfizer_evidence import ENCODING as EVIDENCE_ENCODING
from .pfizer_evidence import SCHEMA as EVIDENCE_SCHEMA
from .pfizer_evidence import EvidenceStore, record_keys
from .pfizer_vendor import DistributionDecoder, NativeQuantification

DRIVER = "orfdata.pfizer-quantification.v1"
METADATA_TABLES = ("DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
                   "EnumDataTypes", "EnumDataTypeValues")
PROJECTION_SCHEMA = """
CREATE TABLE peptide_layout(i INTEGER PRIMARY KEY, identity TEXT NOT NULL UNIQUE);
CREATE TABLE peptide_axis(layout_i INTEGER NOT NULL REFERENCES peptide_layout(i),
    position INTEGER NOT NULL, box_id INTEGER NOT NULL,
    assay_i INTEGER NOT NULL REFERENCES assay(i), PRIMARY KEY(layout_i,position)) WITHOUT ROWID;
CREATE INDEX peptide_axis_assay ON peptide_axis(assay_i,layout_i);
CREATE TABLE peptide_projection(record_i INTEGER NOT NULL, column_name TEXT NOT NULL,
    layout_i INTEGER NOT NULL REFERENCES peptide_layout(i),
    target_i INTEGER NOT NULL REFERENCES target(i), prefix TEXT NOT NULL UNIQUE,
    unit TEXT NOT NULL, PRIMARY KEY(record_i,column_name),
    FOREIGN KEY(record_i,column_name) REFERENCES evidence_vector(record_i,column_name)) WITHOUT ROWID;
CREATE INDEX peptide_projection_target ON peptide_projection(target_i,prefix);
CREATE INDEX peptide_projection_layout ON peptide_projection(layout_i,prefix);
"""


def decoder_from_database(connection):
    metadata = []
    for name in METADATA_TABLES:
        rows = connection.execute(
            "SELECT r.i,r.source_rowid,r.attributes FROM evidence_record r JOIN evidence_table t ON t.i=r.table_i "
            "WHERE t.name=? ORDER BY r.source_rowid", (name,))
        metadata.append([EvidenceStore.restore(connection, row) for row in rows])
    return DistributionDecoder(*metadata)


def copy_evidence(connection, path, *, checkpoint=None):
    """Copy a completed evidence component into a caller-owned bounded build.

    The caller owns disk quotas, transactions and removal of an unsuccessful new
    output. Schema definitions come from this module, never executable source SQL.
    """
    evidence = EvidenceStore(path)
    connection.executescript(EVIDENCE_SCHEMA)
    connection.execute("CREATE TABLE evidence_missing(table_name TEXT NOT NULL, source_key TEXT NOT NULL, "
                       "reason TEXT NOT NULL, PRIMARY KEY(table_name,source_key)) WITHOUT ROWID")
    with closing(evidence.connect()) as original:
        original.execute("BEGIN")
        for table in ("evidence_metadata", "evidence_table", "evidence_record", "evidence_binary",
                      "evidence_vector", "evidence_missing"):
            cursor = original.execute("SELECT * FROM " + table)
            fields = [d[0] for d in cursor.description]
            statement = "INSERT INTO " + table + " VALUES (" + ",".join("?" for _ in fields) + ")"
            count = 0
            while rows := cursor.fetchmany(500):
                connection.executemany(statement, rows)
                connection.commit()
                count += len(rows)
                if checkpoint is not None:
                    checkpoint(dict(table=table, rows=count))


def install_projection(connection, manifest, lineage, *, source_asset, reference=None, checkpoint=None):
    """Install an exact measurement projection in a new dataset with no flat rows.

    A complete normalized evidence component and the checked lineage must already
    be present. Sample/assay/mapping records are the ordinary dataset records.
    No source filtering, protein inference or normalization is performed here.
    The owning builder must publish only after this function succeeds.
    """
    if connection.execute("PRAGMA user_version").fetchone()[0] != 5:
        raise ValueError("Install peptide projection only in a new version-5 metadata dataset")
    if reference is None:
        if len(manifest.references) != 1:
            raise ValueError("Select the native reference explicitly when the manifest includes mapping destinations")
        reference, = manifest.references
    if reference not in manifest.references or reference.taxon_id != 9615:
        raise ValueError("Native canine evidence requires a declared canine reference")
    if source_asset not in {asset.asset_id for asset in manifest.sources}:
        raise ValueError("Native peptide projection source asset is absent from its manifest")
    if connection.execute("SELECT count(*) FROM observation").fetchone()[0]:
        raise ValueError("Native peptide projection cannot replace existing flat observations")
    stored = json.loads(connection.execute("SELECT value FROM metadata WHERE name='manifest'").fetchone()[0])
    if stored != manifest.to_dict():
        raise ValueError("Peptide projection manifest differs from the stored dataset")
    state = connection.execute("SELECT value FROM evidence_metadata WHERE name='state'").fetchone()
    if state is None or json.loads(state[0]) != "complete":
        raise ValueError("Complete the normalized evidence component before projecting measurements")
    decoder = decoder_from_database(connection)
    converter = NativeQuantification(decoder, lineage, reference, manifest.dataset_id, manifest.version,
                                     source_asset=source_asset)
    connection.executescript(PROJECTION_SCHEMA)
    reference_i, = connection.execute("SELECT i FROM reference WHERE key=?", (reference.key,)).fetchone()
    assays = {(identifier, sample): index for index, identifier, sample in connection.execute(
        "SELECT i,id,sample_id FROM assay WHERE reference_key=?", (reference.key,))}
    layouts, logical_count, vector_count = {}, 0, 0
    for table, columns in converter.COLUMNS.items():
        table_id, = connection.execute("SELECT i FROM evidence_table WHERE name=?", (table,)).fetchone()
        for record in connection.execute("SELECT i,source_rowid,attributes FROM evidence_record WHERE table_i=? "
                                         "ORDER BY source_rowid", (table_id,)):
            row = json.loads(record[2])
            target = ScopedTarget(reference, "peptide", row["Sequence"], "amino-acid-sequence")
            key = _target_key(target)
            connection.execute("INSERT OR IGNORE INTO target(reference_i,reference_key,key,kind,namespace,identifier) "
                               "VALUES (?,?,?,?,?,?)", (reference_i, reference.key, key, target.kind,
                                                       target.namespace, target.identifier))
            target_i, = connection.execute("SELECT i FROM target WHERE reference_key=? AND key=?",
                                           (reference.key, key)).fetchone()
            for column in columns:
                axes = converter.peptide_axes(table, row, column)
                axis_rows = tuple((position, axis["box_id"], assays[axis["assay_id"], axis["sample_id"]])
                                  for position, axis in enumerate(axes))
                identity = hashlib.sha256(_json(axis_rows).encode()).hexdigest()
                if identity not in layouts:
                    cursor = connection.execute("INSERT INTO peptide_layout(identity) VALUES (?)", (identity,))
                    layouts[identity] = cursor.lastrowid
                    connection.executemany("INSERT INTO peptide_axis VALUES (?,?,?,?)",
                                           ((cursor.lastrowid, *axis) for axis in axis_rows))
                vector_row = connection.execute(
                    "SELECT value_type,element_count,values_data,presence_data,source_null FROM evidence_vector "
                    "WHERE record_i=? AND column_name=?", (record[0], column)).fetchone()
                if vector_row is None or EvidenceStore._vector(vector_row).count != len(axes):
                    raise ValueError("Native numeric vector does not match its common sample axes")
                source_id = ":".join(str(row[key]) for key in record_keys(table))
                prefix = f"{table}:{source_id}:{column}:"
                connection.execute("INSERT INTO peptide_projection VALUES (?,?,?,?,?,?)",
                                   (record[0], column, layouts[identity], target_i, prefix, f"PD2.4 {table}.{column}"))
                logical_count += len(axes)
                vector_count += 1
            if vector_count % 1500 == 0:
                connection.commit()
                if checkpoint is not None:
                    checkpoint(dict(vectors=vector_count, measurements=logical_count))
    if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
        raise ValueError("Native peptide projection has unresolved identities")
    # Remove unused flat storage so direct SQL fails explicitly instead of
    # suggesting that this dataset has zero measurements. Public queries use
    # the declared version-6 driver; SQLite table layouts are internal.
    for name in ("measurement", "input_measurement"):
        connection.execute("DROP VIEW " + name)
    for name in ("observation_input", "observation"):
        connection.execute("DROP TABLE " + name)
    counts = json.loads(connection.execute("SELECT value FROM metadata WHERE name='counts'").fetchone()[0])
    counts.update(measurement=logical_count, target=connection.execute("SELECT count(*) FROM target").fetchone()[0])
    connection.execute("UPDATE metadata SET value=? WHERE name='counts'", (_json(counts),))
    connection.execute("INSERT INTO metadata VALUES ('measurement_storage',?)", (_json(dict(
        driver=DRIVER, state="complete", source_asset=source_asset, reference_key=reference.key, lineage=lineage,
        vectors=vector_count, measurements=logical_count)),))
    connection.execute("PRAGMA user_version=6")
    connection.execute("ANALYZE")
    connection.commit()
    return dict(vectors=vector_count, measurements=logical_count, layouts=len(layouts))


def validate_storage(connection):
    manifest = json.loads(connection.execute("SELECT value FROM metadata WHERE name='manifest'").fetchone()[0])
    if manifest["processing"].get("importer") == "pfizer-canine-evidence-v1":
        row = connection.execute("SELECT value FROM metadata WHERE name='canine_build'").fetchone()
        if row is None or json.loads(row[0]).get("state") != "complete":
            raise ValueError("Canine dataset build is incomplete")
    row = connection.execute("SELECT value FROM metadata WHERE name='measurement_storage'").fetchone()
    config = json.loads(row[0]) if row else {}
    if config.get("driver") != DRIVER or config.get("state") != "complete":
        raise ValueError("Unsupported or incomplete vector measurement storage")
    evidence = {key: json.loads(value) for key, value in connection.execute(
        "SELECT name,value FROM evidence_metadata WHERE name IN ('state','encoding')")}
    if evidence != dict(state="complete", encoding=EVIDENCE_ENCODING):
        raise ValueError("Unsupported or incomplete normalized evidence component")
    if connection.execute("SELECT 1 FROM reference WHERE key=?", (config["reference_key"],)).fetchone() is None:
        raise ValueError("Native vector reference is absent from its dataset")
    counts = json.loads(connection.execute("SELECT value FROM metadata WHERE name='counts'").fetchone()[0])
    if counts["measurement"] != config["measurements"]:
        raise ValueError("Vector measurement count differs from dataset counts")
    return config


def measurements(dataset, *, scope, target, target_keys, sample_id, assay_id, quantity, status,
                 unit, scale, log_base):
    """Stream the common API projection, with indexed target and assay selection."""
    if quantity not in (None, "peptide_abundance") or scale not in (None, "linear") or log_base is not None:
        return
    if status not in (None, "measured", "missing"):
        return
    with closing(dataset._connect()) as connection:
        converter = getattr(dataset, "_peptide_quantification", None)
        if converter is None:
            config = validate_storage(connection)
            reference = dataset._references[config["reference_key"]]
            converter = NativeQuantification(decoder_from_database(connection), config["lineage"], reference,
                                             dataset.manifest.dataset_id, dataset.manifest.version,
                                             source_asset=config["source_asset"])
            # Dataset releases are immutable. Cache only small source metadata
            # and checked lineage, never rows, vectors or query results.
            dataset._peptide_quantification = converter
        reference = converter.reference
        conditions, params = [], []
        for value, field in ((scope, "t.reference_key"), (_target_key(target) if target else None, "t.key"),
                             (unit, "p.unit")):
            if value is not None:
                conditions.append(field + "=?")
                params.append(value)
        if target_keys is not None:
            conditions.append("t.key IN (SELECT value FROM json_each(?))")
            params.append(_json(target_keys))
        if status == "measured":
            # Every cell of a source-null vector is missing. Skip its expansion
            # in SQL; non-null vectors still need per-cell presence/finite checks.
            # NOT EXISTS only skips an explicitly null vector, preserving the
            # reader's error if a corrupted artifact has lost a required vector.
            conditions.append("NOT EXISTS (SELECT 1 FROM evidence_vector v WHERE v.record_i=p.record_i "
                              "AND v.column_name=p.column_name AND v.source_null=1)")
        if sample_id is not None or assay_id is not None:
            query = "SELECT x.layout_i FROM peptide_axis x JOIN assay a ON a.i=x.assay_i WHERE a.reference_key=?"
            params.append(reference.key)
            for value, field in ((sample_id, "a.sample_id"), (assay_id, "a.id")):
                if value is not None:
                    query += " AND " + field + "=?"
                    params.append(value)
            conditions.append("p.layout_i IN (" + query + ")")
        query = ("SELECT r.i,r.source_rowid,r.attributes,e.name,p.column_name,p.prefix "
                 "FROM peptide_projection p JOIN evidence_record r ON r.i=p.record_i "
                 "JOIN evidence_table e ON e.i=r.table_i JOIN target t ON t.i=p.target_i")
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        query += " ORDER BY p.prefix"
        for record in connection.execute(query, params):
            row = EvidenceStore.restore(connection, record[:3])
            values = converter.peptide_abundances(record[3], row, record[4])
            for value in sorted(values, key=lambda m: m.measurement_id):
                if not value.measurement_id.startswith(record[5]):
                    raise ValueError("Native measurement identity differs from its projection index")
                if ((sample_id is None or value.sample_id == sample_id)
                        and (assay_id is None or value.assay_id == assay_id)
                        and (status is None or value.status == status)):
                    yield value
