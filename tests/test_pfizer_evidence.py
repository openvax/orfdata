"""All selected source facts survive compact evidence storage and completion gates."""

import base64
import hashlib
import json
import sqlite3
import struct
from pathlib import Path

import pytest

from orfdata.sources.pfizer_evidence import (
    EvidenceStore,
    decode_scalars,
    encode_scalars,
    source_key,
)
from tools import normalize_dog_evidence
from tools.export_sqlite_rows import canonical_row, quoted, source_schema
from tools.extraction_budget import ScratchBudget

FIXTURE = Path(__file__).parent / "data/pfizer/vendor/native.json"


@pytest.fixture
def source(tmp_path, monkeypatch, request):
    original = json.loads(FIXTURE.read_text())["tables"]
    memberships = json.loads(FIXTURE.with_name("memberships.json").read_text())["tables"]
    for name in ("TargetProteinGroups", "TargetProteinGroupsTargetProteins"):
        original[name] = memberships[name]
    names = ("DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
             "EnumDataTypes", "EnumDataTypeValues", "WorkflowInputFiles", "TargetProteins", "TargetPeptideGroups", "TargetPsms",
             "TargetProteinGroups", "TargetProteinGroupsTargetProteins")
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
    if getattr(request, "param", None) == "synthetic-memberships":
        # Deliberately synthetic shared/dangling assignments for common-API
        # regression. These links are not asserted as biological source facts.
        peptide, psm = tables["TargetPeptideGroups"][0], tables["TargetPsms"][0]
        proteins = [row["UniqueSequenceID"] for row in tables["TargetProteins"][:2]]
        tables["TargetPeptideGroupsTargetProteins"] = [
            dict(TargetPeptideGroupsPeptideGroupID=peptide["PeptideGroupID"], TargetProteinsUniqueSequenceID=identifier)
            for identifier in [*proteins, 999999]]
        tables["TargetPsmsTargetProteins"] = [
            dict(TargetPsmsWorkflowID=psm["WorkflowID"], TargetPsmsPeptideID=psm["PeptideID"],
                 TargetProteinsUniqueSequenceID=identifier) for identifier in proteins]
    path = tmp_path / "source.sqlite"
    db = sqlite3.connect(path)
    plan = dict(scope="Explicit regression slice with synthetic acquisition ordering and optional synthetic memberships", tables=[])
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


def test_native_absent_reporters_and_present_null_vectors_remain_distinct(source, tmp_path):
    fixture = json.loads(FIXTURE.with_name("acquisition-missingness.json").read_text())
    native_plan = json.loads((Path(__file__).parents[1] / "tools/dog_result_plan.json").read_text())
    selected = {item["table"]: item for item in native_plan["tables"]}
    path, _ = source
    originals = {}
    with sqlite3.connect(path) as db:
        configuration = json.loads(db.execute("SELECT value FROM _export_metadata WHERE key='configuration'").fetchone()[0])
        for name, records in fixture["metadata_additions"].items():
            item = next(t for t in configuration["plan"]["tables"] if t["table"] == name)
            # The old probe fixture also contained UI fields outside the actual
            # extraction selection. Use the shared selected metadata columns.
            for column in set(item["columns"]) - set(records[0]):
                db.execute(f"ALTER TABLE {quoted(name)} DROP COLUMN {quoted(column)}")
            columns = item["columns"] = [c for c in item["columns"] if c in records[0]]
            last_rowid = db.execute(f"SELECT max(_source_rowid) FROM {quoted(name)}").fetchone()[0]
            for index, row in enumerate(records, last_rowid + 1):
                values = index, *(row[c] for c in columns)
                db.execute(f"INSERT INTO {quoted(name)} VALUES (" + ",".join("?" for _ in values) + ")", values)
            digest, count = bytes(32), 0
            for row in db.execute(f"SELECT * FROM {quoted(name)} ORDER BY _source_rowid"):
                digest = hashlib.sha256(digest + canonical_row(row)).digest()
                count += 1
            db.execute("UPDATE _export_progress SET rows=?,digest=? WHERE name=?", (count, digest.hex(), name))
        for name, records in fixture["tables"].items():
            item = dict(selected[name], select_from=[dict(table="explicit regression selection")])
            columns = item["columns"]
            rows = [{k: base64.b64decode(v["base64"]) if isinstance(v, dict) and v.get("type") == "blob"
                     else v for k, v in row.items()} for row in records]
            originals[name] = rows
            configuration["plan"]["tables"].append(item)
            db.execute(f"CREATE TABLE {quoted(name)}(_source_rowid INTEGER PRIMARY KEY," +
                       ",".join(quoted(c) for c in columns) + ")")
            digest = bytes(32)
            for row in sorted(rows, key=lambda r: (*[r[k] for k in item["keys"]], r["_source_rowid"])):
                values = row["_source_rowid"], *(row[c] for c in columns)
                digest = hashlib.sha256(digest + canonical_row(values)).digest()
                db.execute(f"INSERT INTO {quoted(name)} VALUES (" + ",".join("?" for _ in values) + ")", values)
            assert digest.hex() == fixture["checkpoints"][name]["typed_sha256"]
            missing = [r for r in fixture["missing_acquisitions"] if r["name"] == name]
            db.execute("INSERT INTO _export_progress VALUES (?,?,?,?,?)", (name, len(rows), digest.hex(), 1, len(missing)))
            db.executemany("INSERT INTO _export_missing VALUES (?,?,?)",
                           ((r["name"], r["source_key"], r["reason"]) for r in missing))
        db.execute("UPDATE _export_metadata SET value=? WHERE key='configuration'", (json.dumps(configuration),))
    normalize_dog_evidence.PLAN_PATH.write_text(json.dumps(configuration["plan"]))
    output, _ = build(path, tmp_path)
    evidence = EvidenceStore(output)
    for name, rows in originals.items():
        assert list(evidence.records(name)) == rows
    table = "QuanSpectrumInfoReporterQuanResults"
    expected = [dict(table=table, key=(-160, 8546281), reason="No native row for the referenced key")]
    assert list(evidence.missing_acquisitions()) == expected
    assert list(evidence.missing_acquisitions(table, key=(-160, 8546281))) == expected
    assert list(evidence.missing_acquisitions("ReporterQuanResults")) == []
    assert list(evidence.missing_acquisitions(table, key=(-160, 8546311))) == []
    # Both null cases retain the native No Quan Values flag. Only one lacks links.
    for identifier in (8546281, 8546311):
        record, = evidence.records("QuanSpectrumInfo", key=(-160, identifier))
        assert record["QuanInfo"] == 1
        assert record["Abundances"] is record["AbundancesNormalized"] is None
        assert evidence.vector("QuanSpectrumInfo", (-160, identifier), "Abundances").source_null
    links = list(evidence.records(table))
    assert {r["QuanSpectrumInfoSpectrumID"] for r in links} == {8546263, 8546311}
    assert len([r for r in links if r["QuanSpectrumInfoSpectrumID"] == 8546311]) == 11
    assert not evidence.vector("QuanSpectrumInfo", (-160, 8546263), "Abundances").source_null
    with pytest.raises(ValueError, match="Unknown evidence table"):
        list(evidence.missing_acquisitions("unknown"))
    with pytest.raises(TypeError, match="native integers"):
        list(evidence.missing_acquisitions(table, key=(-160, 8546281.0)))


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


def test_compact_scalar_presence_preserves_null_signed_zero_and_int64():
    columns = ("native_id", "zero", "missing", "external_bytes")
    attributes = dict(native_id=-9222281927842617047, zero=-0.0, missing=None)
    payload = encode_scalars(attributes, columns)
    # Check the specified wire layout independently of the production decoder.
    assert json.loads(payload) == ["7", [-9222281927842617047, -0.0, None]]
    restored = decode_scalars(payload, columns)
    assert restored == attributes and "external_bytes" not in restored
    assert struct.pack("<d", restored["zero"]) == struct.pack("<d", -0.0)
    assert decode_scalars('["0",[]]', columns) == {}
    assert decode_scalars(json.dumps(attributes), columns) == restored


@pytest.mark.parametrize("payload", [
    '["10",[1]]', '["3",[1]]', '["1",[1,2]]', '["00",[]]', '["-1",[]]',
    '[1,[2]]', '["1",{}]', '["1",[[1]]]', '["1",[NaN]]',
    '["1",[Infinity]]', 'null', '{"unknown":1}',
])
def test_malformed_scalar_layout_cannot_change_or_hide_fields(payload):
    with pytest.raises(ValueError, match="[Ss]calar"):
        decode_scalars(payload, ("one", "two", "three", "four"))


def test_legacy_object_component_replays_every_native_field(source, tmp_path):
    output, _ = build(source[0], tmp_path)
    expected = {name: list(EvidenceStore(output).records(name)) for name in source[1]}
    with sqlite3.connect(output) as db:
        layouts = {index: json.loads(columns) for index, columns in db.execute("SELECT i,columns_json FROM evidence_table")}
        for identifier, table, payload in db.execute("SELECT i,table_i,attributes FROM evidence_record").fetchall():
            attributes = decode_scalars(payload, layouts[table])
            db.execute("UPDATE evidence_record SET attributes=? WHERE i=?",
                       (json.dumps(attributes, sort_keys=True, separators=(",", ":")), identifier))
        db.execute("UPDATE evidence_metadata SET value='\"orfdata.pfizer-evidence.v1\"' WHERE name='encoding'")
    legacy = EvidenceStore(output)
    for name, rows in expected.items():
        assert list(legacy.records(name)) == rows


@pytest.mark.parametrize("corruption", ["missing_vector", "duplicate_vector"])
def test_scalar_presence_cannot_conceal_missing_or_duplicate_external_fields(source, tmp_path, corruption):
    output, _ = build(source[0], tmp_path)
    with sqlite3.connect(output) as db:
        identifier, payload, layout = db.execute(
            "SELECT r.i,r.attributes,t.columns_json FROM evidence_record r JOIN evidence_table t ON t.i=r.table_i "
            "WHERE t.name='TargetPsms' ORDER BY r.source_rowid LIMIT 1").fetchone()
        if corruption == "missing_vector":
            db.execute("DELETE FROM evidence_vector WHERE record_i=? AND column_name='Abundances'", (identifier,))
        else:
            columns = json.loads(layout)
            attributes = decode_scalars(payload, columns)
            attributes["Abundances"] = None
            db.execute("UPDATE evidence_record SET attributes=? WHERE i=?", (encode_scalars(attributes, columns), identifier))
    with pytest.raises(ValueError, match="Stored evidence field"):
        list(EvidenceStore(output).records("TargetPsms"))
