"""All selected source facts survive compact evidence storage and completion gates."""

import base64
import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from orfdata.sources.pfizer_evidence import EvidenceStore, source_key
from tools import normalize_dog_evidence
from tools.export_sqlite_rows import canonical_row, quoted, source_schema
from tools.extraction_budget import ScratchBudget

FIXTURE = Path(__file__).parent / "data/pfizer/vendor/native.json"


@pytest.fixture
def source(tmp_path, monkeypatch):
    original = json.loads(FIXTURE.read_text())["tables"]
    names = ("DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
             "EnumDataTypes", "EnumDataTypeValues", "TargetProteins", "TargetPeptideGroups", "TargetPsms")
    tables = {name: [{key: base64.b64decode(value["base64"]) if isinstance(value, dict)
                     and value.get("type") == "blob" else value for key, value in row.items()}
                    for row in original[name]] for name in names}
    # Explicit synthetic acquisitions exercise composite order different from
    # source-rowid order, plus nonfinite scalars, opaque bytes and signed IDs.
    tables["AcquisitionFixture"] = [
        dict(workflow=-2, spectrum=2, value=float("inf"), opaque=b"\xff\0", native_id=-9222281927842617047),
        dict(workflow=-3, spectrum=1, value=0.0, opaque=b"", native_id=9),
        dict(workflow=-2, spectrum=1, value=None, opaque=None, native_id=10),
    ]
    path = tmp_path / "source.sqlite"
    db = sqlite3.connect(path)
    plan = dict(scope="Explicit native regression slice with synthetic acquisition ordering", tables=[])
    progress = []
    for name, rows in tables.items():
        columns = list(rows[0])
        item = dict(table=name, columns=columns, role="target" if name.startswith("Target") else "metadata")
        if name == "AcquisitionFixture":
            item.update(keys=["workflow", "spectrum"], select_from=[dict(table="fixture selection")], role="referenced_acquisition")
        plan["tables"].append(item)
        db.execute(f"CREATE TABLE {quoted(name)}(_source_rowid INTEGER PRIMARY KEY," + ",".join(quoted(c) for c in columns) + ")")
        for index, row in enumerate(rows, 1):
            row["_source_rowid"] = index
            db.execute(f"INSERT INTO {quoted(name)} VALUES (" + ",".join("?" for _ in range(len(columns) + 1)) + ")",
                       (index, *(row[c] for c in columns)))
        ordered = sorted(rows, key=lambda r: (r["workflow"], r["spectrum"])) if item.get("select_from") else rows
        digest = bytes(32)
        for row in ordered:
            digest = hashlib.sha256(digest + canonical_row((row["_source_rowid"], *(row[c] for c in columns)))).digest()
        progress.append((name, len(rows), digest.hex(), 1, 0))
    plan["source_schema_sha256"] = source_schema(db)[1]
    db.execute("CREATE TABLE _export_metadata(key,value)")
    db.execute("INSERT INTO _export_metadata VALUES ('configuration',?)", (json.dumps(dict(plan=plan, identity="explicit fixture")),))
    db.execute("CREATE TABLE _export_progress(name,rows,digest,complete,missing_keys)")
    db.executemany("INSERT INTO _export_progress VALUES (?,?,?,?,?)", progress)
    db.execute("CREATE TABLE _export_missing(name,source_key,reason)")
    db.commit()
    db.close()
    selection = tmp_path / "fixture-selection.json"
    selection.write_text(json.dumps(plan))
    monkeypatch.setattr(normalize_dog_evidence, "PLAN_PATH", selection)
    return path, tables


def build(source, tmp_path, *, cap=2 * 1024 * 1024):
    root = tmp_path / "normalized"
    budget = ScratchBudget(root, 8 * 1024 * 1024)
    output = root / "evidence.sqlite"
    result = normalize_dog_evidence.normalize(source, output, budget, cap, batch_rows=2)
    return output, result


def test_normalized_evidence_reconstructs_every_selected_field_and_native_checkpoint(source, tmp_path):
    path, tables = source
    output, result = build(path, tmp_path)
    evidence = EvidenceStore(output)
    assert evidence.metadata["state"] == "complete"
    assert result["counts"] == {name: len(rows) for name, rows in tables.items()}
    for name, rows in tables.items():
        assert list(evidence.records(name)) == rows
    assert list(evidence.records("TargetPsms", key=(-160, 1)))[0]["Sequence"] == "GHEEEHHHR"
    assert len(list(evidence.records("TargetPsms", sequence="GHEEEHHHR"))) == 1
    vector = evidence.vector("TargetPsms", (-160, 1), "Abundances")
    assert vector.value_type == "float64" and vector.count == 11
    assert list(vector.cells())[0]["value"] == 7496.059512428747
    assert list(vector.cells())[1]["status"] == "missing"
    assert source_key((-9222281927842617047,)) == '["-9222281927842617047"]'
    with sqlite3.connect(output) as db:
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
        assert db.execute("SELECT count(*) FROM evidence_vector WHERE source_null=1").fetchone()[0] > 0
        assert db.execute("SELECT count(*) FROM evidence_binary WHERE value_kind='float64'").fetchone()[0] == 1
        assert db.execute("SELECT count(*) FROM sqlite_master WHERE name='normalization_order'").fetchone()[0] == 0


@pytest.mark.parametrize("table", ["TargetPsms", "AcquisitionFixture"])
def test_changed_checkpoint_digest_is_rejected_for_rowid_and_composite_orders(source, tmp_path, table):
    path, _ = source
    with sqlite3.connect(path) as db:
        db.execute("UPDATE _export_progress SET digest=? WHERE name=?", ("f" * 64, table))
    with pytest.raises(ValueError, match="checkpoint digest"):
        build(path, tmp_path)
    assert not (tmp_path / "normalized/evidence.sqlite").exists()


def test_incomplete_selection_missing_link_disagreement_and_small_disk_cap_are_rejected(source, tmp_path):
    path, _ = source
    with sqlite3.connect(path) as db:
        db.execute("UPDATE _export_progress SET complete=0 WHERE name='TargetPsms'")
    with pytest.raises(ValueError, match="Complete all selected"):
        build(path, tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE _export_progress SET complete=1")
        db.execute("UPDATE _export_progress SET missing_keys=1 WHERE name='TargetPsms'")
    with pytest.raises(ValueError, match="missing-link reports"):
        build(path, tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE _export_progress SET missing_keys=0")
    with pytest.raises(sqlite3.OperationalError, match="full"):
        build(path, tmp_path, cap=32768)
    assert not (tmp_path / "normalized/evidence.sqlite").exists()


def test_reader_refuses_building_component_and_normalizer_never_overwrites_existing_output(source, tmp_path):
    path, _ = source
    output, _ = build(path, tmp_path)
    original_hash = hashlib.sha256(output.read_bytes()).hexdigest()
    with pytest.raises(FileExistsError):
        build(path, tmp_path)
    assert hashlib.sha256(output.read_bytes()).hexdigest() == original_hash
    with sqlite3.connect(output) as db:
        db.execute("UPDATE evidence_metadata SET value='\"building\"' WHERE name='state'")
    with pytest.raises(ValueError, match="incomplete"):
        EvidenceStore(output)
