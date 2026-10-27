"""Indexed common measurements projected from lossless canine evidence vectors.

Storage version 6 keeps numeric arrays and shared sample axes once. It has no
flat observation table: the public dataset API expands requested cells, including
missing cells, and reports the complete logical observation count. All source
facts remain accessible through EvidenceStore; no vendor file or network access
is needed to open or query this representation.
"""

import hashlib
import json
import sqlite3
from contextlib import closing

from ..dataset import _json, _target_key
from ..reference import ScopedTarget
from .pfizer_evidence import SCHEMA as EVIDENCE_SCHEMA
from .pfizer_evidence import (
    SUPPORTED_ENCODINGS,
    EvidenceStore,
    column_layouts,
    decode_scalars,
    record_keys,
)
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
    layouts = column_layouts(connection)
    for name in METADATA_TABLES:
        rows = connection.execute(
            "SELECT r.i,r.source_rowid,r.attributes FROM evidence_record r JOIN evidence_table t ON t.i=r.table_i "
            "WHERE t.name=? ORDER BY r.source_rowid", (name,))
        metadata.append([EvidenceStore.restore(connection, row, columns=layouts[name]) for row in rows])
    return DistributionDecoder(*metadata)


def copy_evidence(connection, path, *, checkpoint=None, defer_checks=False):
    """Seed a new dataset with verified evidence pages and existing metadata.

    Only the known normalized schema can be copied. The destination must hold
    at most 16 MiB of fresh common metadata, saved as SQL from our own database;
    no SQL from the source is executed. The caller owns the disposable output,
    its limits and failure cleanup. Copying pages retains existing indexes.
    defer_checks requires the owning builder to check all foreign keys before
    publishing the completed database.
    """
    evidence = EvidenceStore(path)
    page_size = connection.execute("PRAGMA page_size").fetchone()[0]
    cap = page_size * connection.execute("PRAGMA max_page_count").fetchone()[0]
    if page_size * connection.execute("PRAGMA page_count").fetchone()[0] > 16 * 1024**2:
        raise ValueError("Canine evidence seeding requires bounded fresh metadata")
    # ANALYZE's implementation table may already exist in the evidence seed.
    metadata = "\n".join(statement for statement in connection.iterdump()
                         if not statement.startswith(('ANALYZE ', 'INSERT INTO "sqlite_stat1"')))
    if len(metadata.encode()) > 32 * 1024**2:
        raise ValueError("Canine metadata SQL exceeds its memory bound")
    journal = connection.execute("PRAGMA journal_mode").fetchone()[0]
    user_version = connection.execute("PRAGMA user_version").fetchone()[0]
    foreign_keys = connection.execute("PRAGMA foreign_keys").fetchone()[0]
    if connection.in_transaction:
        raise ValueError("Finish metadata writes before seeding native evidence")

    def schema(db):
        return {(kind, name): "".join(sql.split()) if sql else None
                for kind, name, sql in db.execute("SELECT type,name,sql FROM sqlite_master")}

    with closing(evidence.connect()) as original, closing(sqlite3.connect(":memory:")) as expected:
        expected.executescript(EVIDENCE_SCHEMA)
        expected.execute("CREATE TABLE evidence_missing(table_name TEXT NOT NULL, source_key TEXT NOT NULL, "
                         "reason TEXT NOT NULL, PRIMARY KEY(table_name,source_key)) WITHOUT ROWID")
        expected.execute("ANALYZE")
        original.execute("BEGIN")
        if schema(original) != schema(expected):
            raise ValueError("Normalized evidence schema differs from the supported schema")
        source_pages = original.execute("PRAGMA page_count").fetchone()[0]
        source_page_size = original.execute("PRAGMA page_size").fetchone()[0]
        if source_page_size != page_size or source_pages * source_page_size > cap:
            raise ValueError("Normalized evidence pages exceed the output limit or differ in page size")

        def progress(status, remaining, total):
            if total != source_pages:
                raise ValueError("Normalized evidence snapshot changed during backup")
            if checkpoint is not None:
                checkpoint(dict(stage="evidence_backup", pages=total-remaining, total_pages=total))

        connection.execute("PRAGMA foreign_keys=OFF")
        original.backup(connection, pages=4096, progress=progress)
    # Backup restores the source header; reapply the caller's size/storage policy
    # before adding the bounded common metadata and derived indexes.
    connection.execute(f"PRAGMA max_page_count={cap // page_size}")
    connection.execute("PRAGMA journal_mode=" + journal)
    connection.executescript(metadata)
    connection.execute(f"PRAGMA user_version={user_version}")
    connection.execute(f"PRAGMA foreign_keys={foreign_keys}")
    if not defer_checks and connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
        raise ValueError("Seeded canine metadata/evidence has unresolved identities")


def install_projection(connection, manifest, lineage, *, source_asset, reference=None, checkpoint=None, defer_checks=False):
    """Install an exact measurement projection in a new dataset with no flat rows.

    A complete normalized evidence component and the checked lineage must already
    be present. Sample/assay/mapping records are the ordinary dataset records.
    No source filtering, protein inference or normalization is performed here.
    The owning builder must publish only after this function succeeds. With
    defer_checks it must run the final foreign-key check and ANALYZE itself.
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
    scalar_layouts = column_layouts(connection)
    for table, columns in converter.COLUMNS.items():
        table_id, = connection.execute("SELECT i FROM evidence_table WHERE name=?", (table,)).fetchone()
        for record in connection.execute("SELECT i,source_rowid,attributes FROM evidence_record WHERE table_i=? "
                                         "ORDER BY source_rowid", (table_id,)):
            row = decode_scalars(record[2], scalar_layouts[table])
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
    if not defer_checks and connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
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
    if not defer_checks:
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
    if evidence.get("state") != "complete" or evidence.get("encoding") not in SUPPORTED_ENCODINGS:
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
        scalar_layouts = column_layouts(connection)
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
            row = EvidenceStore.restore(connection, record[:3], columns=scalar_layouts[record[3]])
            values = converter.peptide_abundances(record[3], row, record[4])
            for value in sorted(values, key=lambda m: m.measurement_id):
                if not value.measurement_id.startswith(record[5]):
                    raise ValueError("Native measurement identity differs from its projection index")
                if ((sample_id is None or value.sample_id == sample_id)
                        and (assay_id is None or value.assay_id == assay_id)
                        and (status is None or value.status == status)):
                    yield value
