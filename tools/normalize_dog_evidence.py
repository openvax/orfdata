"""Normalize a completed bounded native export, retaining all selected facts."""

import argparse
import hashlib
import json
import resource
import sqlite3
import sys
from pathlib import Path

from orfdata.sources.pfizer_evidence import EvidenceStore, EvidenceWriter
from orfdata.sources.pfizer_vendor import DistributionDecoder
from tools.export_sqlite_rows import canonical_row, quoted
from tools.extraction_budget import ScratchBudget

PLAN_PATH = Path(__file__).with_name("dog_result_plan.json")
CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
SOURCE_SHA256 = {name: hashlib.sha256((Path(__file__).resolve().parents[1] / name).read_bytes()).hexdigest()
                 for name in ("tools/normalize_dog_evidence.py", "tools/export_sqlite_rows.py", "tools/extraction_budget.py",
                              "orfdata/numeric_vector.py", "orfdata/sources/pfizer_vendor.py",
                              "orfdata/sources/pfizer_evidence.py")}


def digest_file(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def normalize(source, output, scratch, output_limit, *, batch_rows=500):
    """Build a complete, bit-checked evidence component; never acquire source bytes.

    The input must have finished every table in the pinned extraction plan.
    Original rows and independently restored normalized rows get matching typed
    digests in source-rowid order. Acquisition tables additionally retain their
    extractor digest, whose composite-key order differs from rowid order.
    """
    source, output = Path(source).resolve(), Path(output).resolve()
    if type(batch_rows) is not int or batch_rows < 1:
        raise ValueError("Normalization batch size must be positive")
    plan = json.loads(PLAN_PATH.read_text())
    original = sqlite3.connect(source.as_uri() + "?mode=ro", uri=True)
    destination = None
    created = False
    try:
        original.execute("PRAGMA query_only=ON")
        original.execute("PRAGMA cache_size=-8192")
        original.execute("PRAGMA hard_heap_limit=536870912")
        original.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
        original.execute("BEGIN")
        configuration = json.loads(original.execute("SELECT value FROM _export_metadata WHERE key='configuration'").fetchone()[0])
        history_row = original.execute("SELECT value FROM _export_metadata WHERE key='selection_history'").fetchone()
        selection_history = json.loads(history_row[0]) if history_row else []
        if configuration["plan"] != plan:
            raise ValueError("Native source selection differs from the pinned dog extraction plan")
        progress = {name: dict(rows=count, digest=digest, complete=complete, missing_keys=missing)
                    for name, count, digest, complete, missing in original.execute(
                        "SELECT name,rows,digest,complete,missing_keys FROM _export_progress")}
        if set(progress) != {item["table"] for item in plan["tables"]} or not all(r["complete"] == 1 for r in progress.values()):
            raise ValueError("Complete all selected native tables before normalization")
        source_hash = digest_file(source)
        original.row_factory = sqlite3.Row
        metadata = {name: [dict(row) for row in original.execute(f"SELECT * FROM {quoted(name)} ORDER BY _source_rowid")]
                    for name in ("DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
                                 "EnumDataTypes", "EnumDataTypeValues")}
        decoder = DistributionDecoder(*(metadata[name] for name in metadata))
        scratch.reserve_sqlite(output, output_limit)
        # Never clobber a user's existing evidence component or interrupted build.
        output.touch(exist_ok=False)
        created = True
        destination = sqlite3.connect(output)
        destination.execute("PRAGMA page_size=4096")
        destination.execute(f"PRAGMA max_page_count={output_limit // 4096}")
        destination.execute("PRAGMA foreign_keys=ON")
        destination.execute("PRAGMA synchronous=FULL")
        destination.execute("PRAGMA journal_mode=DELETE")
        destination.execute("PRAGMA temp_store=MEMORY")
        destination.execute("PRAGMA cache_size=-8192")
        destination.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
        writer = EvidenceWriter(destination, plan, decoder)
        destination.execute("INSERT INTO evidence_metadata VALUES ('source',?)", (json.dumps(dict(
            source_file_sha256=source_hash, selection=configuration, selection_history=selection_history, extraction_progress=progress,
            integrity_scope="Lossless selected native facts; the vendor database has no authoritative full-file checksum"),
            sort_keys=True, separators=(",", ":")),))
        destination.execute("CREATE TABLE evidence_missing(table_name TEXT NOT NULL, source_key TEXT NOT NULL, "
                            "reason TEXT NOT NULL, PRIMARY KEY(table_name,source_key)) WITHOUT ROWID")
        missing_counts = {name: 0 for name in progress}
        for name, key, reason in original.execute("SELECT name,source_key,reason FROM _export_missing ORDER BY name,source_key"):
            if name not in missing_counts:
                raise ValueError("Missing-link report names an unselected source table")
            missing_counts[name] += 1
            destination.execute("INSERT INTO evidence_missing VALUES (?,?,?)", (name, key, reason))
        if any(count != progress[name]["missing_keys"] for name, count in missing_counts.items()):
            raise ValueError("Native missing-link reports differ from extraction checkpoints")
        destination.commit()
        counts, hashes = {}, {}
        for item in plan["tables"]:
            name, columns = item["table"], item["columns"]
            count, digest = 0, bytes(32)
            ordering_keys = item.get("keys", ()) if item.get("select_from") else ()
            if ordering_keys:
                keys = ["k" + str(i) for i in range(len(ordering_keys))]
                fields = ",".join(k + " INTEGER NOT NULL" for k in keys)
                destination.execute("CREATE TABLE normalization_order(" + fields + ",source_rowid INTEGER NOT NULL, "
                                    "PRIMARY KEY(" + ",".join(keys) + ",source_rowid)) WITHOUT ROWID")
            for record in original.execute(f"SELECT * FROM {quoted(name)} ORDER BY _source_rowid"):
                row = dict(record)
                writer.insert(name, row)
                if ordering_keys:
                    values = tuple(row[key] for key in ordering_keys) + (row["_source_rowid"],)
                    if any(type(value) is not int for value in values):
                        raise ValueError("Native acquisition sort keys must remain int64 values")
                    destination.execute("INSERT INTO normalization_order VALUES (" + ",".join("?" for _ in values) + ")", values)
                digest = hashlib.sha256(digest + canonical_row((row["_source_rowid"], *(row[c] for c in columns)))).digest()
                count += 1
                if count % batch_rows == 0:
                    destination.commit()
                    scratch.usage()
            destination.commit()
            if count != progress[name]["rows"]:
                raise ValueError("Native checkpoint count differs from actual selected rows: " + name)
            if not item.get("select_from") and digest.hex() != progress[name]["digest"]:
                raise ValueError("Native checkpoint digest differs from actual selected rows: " + name)
            restored_digest = bytes(32)
            for row in destination.execute("SELECT i,source_rowid,attributes FROM evidence_record WHERE table_i=? ORDER BY source_rowid",
                                           (writer.tables[name][0],)):
                restored = EvidenceStore.restore(destination, row, columns=columns)
                restored_digest = hashlib.sha256(restored_digest + canonical_row(
                    (restored["_source_rowid"], *(restored[c] for c in columns)))).digest()
            if restored_digest != digest:
                raise ValueError("Normalized evidence does not reconstruct the original typed rows: " + name)
            if ordering_keys:
                # A disk-bounded b-tree avoids a large unbounded temp sort. These
                # are the exact numeric key ordering rules used by NativeExport.
                extraction_digest = bytes(32)
                for (rowid,) in destination.execute("SELECT source_rowid FROM normalization_order ORDER BY " +
                                                     ",".join(keys) + ",source_rowid"):
                    stored = destination.execute("SELECT i,source_rowid,attributes FROM evidence_record "
                                                 "WHERE table_i=? AND source_rowid=?", (writer.tables[name][0], rowid)).fetchone()
                    restored = EvidenceStore.restore(destination, stored, columns=columns)
                    extraction_digest = hashlib.sha256(extraction_digest + canonical_row(
                        (rowid, *(restored[c] for c in columns)))).digest()
                if extraction_digest.hex() != progress[name]["digest"]:
                    raise ValueError("Native acquisition checkpoint digest differs from restored source facts: " + name)
                destination.execute("DROP TABLE normalization_order")
                destination.commit()
            counts[name], hashes[name] = count, digest.hex()
            print(json.dumps(dict(table=name, rows=count, restored_native_sha256=digest.hex())), flush=True)
        if destination.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("Normalized evidence has invalid references")
        if destination.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Normalized evidence failed SQLite integrity checks")
        if digest_file(source) != source_hash:
            raise ValueError("Selected native input changed during normalization")
        result = dict(counts=counts, native_rowid_order_sha256=hashes, source_file_sha256=source_hash,
                      logical_sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest())
        destination.execute("INSERT INTO evidence_metadata VALUES ('reconciliation',?)", (json.dumps(result, sort_keys=True),))
        destination.execute("INSERT INTO evidence_metadata VALUES ('tool_sha256',?)", (json.dumps(CODE_SHA256),))
        destination.execute("INSERT INTO evidence_metadata VALUES ('software_sha256',?)", (json.dumps(SOURCE_SHA256, sort_keys=True),))
        destination.execute("UPDATE evidence_metadata SET value='\"complete\"' WHERE name='state'")
        destination.execute("ANALYZE")
        destination.commit()
        scratch.usage()
        return result
    except BaseException:
        if destination is not None:
            destination.close()
            destination = None
        if created:
            output.unlink()
        raise
    finally:
        if destination is not None:
            destination.close()
        original.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--output-mib", type=int, default=2048)
    parser.add_argument("--scratch-mib", type=int, default=8192)
    args = parser.parse_args()
    root = args.output_directory.resolve()
    scratch = ScratchBudget(root, args.scratch_mib * 1024 * 1024)
    result = normalize(args.source, root / "evidence.sqlite", scratch, args.output_mib * 1024 * 1024)
    result["output_file_sha256"] = digest_file(root / "evidence.sqlite")
    result["scratch"] = scratch.receipt()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["peak_process_rss_bytes"] = peak if sys.platform == "darwin" else peak * 1024
    result["runtime"] = dict(python=sys.version, sqlite=sqlite3.sqlite_version)
    result["tool_sha256"] = CODE_SHA256
    result["software_sha256"] = SOURCE_SHA256
    scratch.atomic_write(root / "normalization-receipt.json", (json.dumps(result, sort_keys=True, indent=2) + "\n").encode())
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
