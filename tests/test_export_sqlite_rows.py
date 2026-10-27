"""Native export checkpoints preserve exact values, scope and source identities."""

import copy
import json
import sqlite3

import pytest

from tools.export_sqlite_rows import NativeExport, source_schema
from tools.extraction_budget import LimitExceeded, ScratchBudget


@pytest.fixture
def native(tmp_path):
    database = sqlite3.connect(tmp_path / "native.sqlite")
    database.executescript("""
        CREATE TABLE Peptides(id,sequence,quantity,vector,unused_spectrum);
        CREATE TABLE Membership(workflow,spectrum);
        CREATE TABLE Acquisition(workflow,spectrum,channel,height,unused_spectrum);
        CREATE INDEX acquisition_key ON Acquisition(workflow,spectrum);
    """)
    database.executemany("INSERT INTO Peptides(rowid,id,sequence,quantity,vector,unused_spectrum) VALUES (?,?,?,?,?,?)", [
        (-7, -9222281927842617047, "SHARED", None, b"\0\xff", b"excluded" * 2000),
        (2, 13, "PEPTIDE", 0.0, None, b"excluded" * 2000),
        (50, 14, "OTHER", float("inf"), b"", b"excluded" * 2000),
    ])
    database.executemany("INSERT INTO Membership VALUES (?,?)", [(1, 1), (1, 1), (1, 2), (2, 99)])
    database.executemany("INSERT INTO Acquisition VALUES (?,?,?,?,?)", [
        (1, 1, 1, 7.0, b"not exported"), (1, 1, 2, None, b"not exported"),
        (1, 2, 1, 0.0, b"not exported"), (1, 3, 1, 999.0, b"unreferenced"),
    ])
    database.commit()
    plan = dict(source_schema_sha256=source_schema(database)[1], tables=[
        dict(table="Peptides", columns=["id", "sequence", "quantity", "vector"]),
        dict(table="Membership", columns=["workflow", "spectrum"]),
        dict(table="Acquisition", columns=["workflow", "spectrum", "channel", "height"],
             keys=["workflow", "spectrum"], select_from=[dict(table="Membership", columns=["workflow", "spectrum"])])])
    yield database, plan
    database.close()


def exporter(native, directory, **kwargs):
    source, plan = native
    directory.mkdir(exist_ok=True)
    return NativeExport(source, directory / "output.sqlite", plan, {"asset": "pinned"},
                        ScratchBudget(directory, 16 * 1024 * 1024), output_limit=1024 * 1024, **kwargs)


def test_resume_retains_negative_rowids_int64_missing_and_binary_without_spectra(native, tmp_path):
    output = exporter(native, tmp_path / "export", batch_rows=1)
    state = output.copy("Peptides", max_batches=1)
    assert state["rows"] == 1 and state["cursor"] == -7 and not state["complete"]
    output.close()
    resumed = exporter(native, tmp_path / "export", batch_rows=2)
    assert resumed.copy("Peptides")["complete"]
    expected = list(native[0].execute("SELECT rowid,id,sequence,quantity,vector FROM Peptides ORDER BY rowid"))
    assert list(resumed.database.execute("SELECT * FROM Peptides ORDER BY _source_rowid")) == expected
    assert "unused_spectrum" not in {row[1] for row in resumed.database.execute("PRAGMA table_info(Peptides)")}
    assert resumed.database.execute("SELECT typeof(id) FROM Peptides WHERE _source_rowid=-7").fetchone()[0] == "integer"
    assert resumed.database.execute("SELECT typeof(quantity) FROM Peptides WHERE _source_rowid=2").fetchone()[0] == "real"
    resumed.close()


def test_acquisition_is_selected_by_all_distinct_keys_and_missing_links_are_explicit(native, tmp_path):
    output = exporter(native, tmp_path / "export", lookup_keys=1)
    with pytest.raises(ValueError, match="Complete native link"):
        output.copy("Acquisition")
    output.copy("Membership")
    output.copy("Acquisition", max_batches=1)
    output.close()
    resumed = exporter(native, tmp_path / "export", lookup_keys=2)
    result = resumed.copy("Acquisition")
    assert result["complete"] and result["rows"] == 3 and result["selected_keys"] == 3
    assert result["missing_keys"] == 1
    assert list(resumed.database.execute("SELECT workflow,spectrum,channel,height FROM Acquisition")) == [
        (1, 1, 1, 7.0), (1, 1, 2, None), (1, 2, 1, 0.0)]
    assert list(resumed.database.execute("SELECT source_key FROM _export_missing")) == [("[2,99]",)]
    resumed.close()


def test_interruption_between_rows_and_checkpoint_rolls_back_both(native, tmp_path):
    output = exporter(native, tmp_path / "export", batch_rows=1)
    output.database.execute("CREATE TRIGGER interruption BEFORE UPDATE ON _export_progress "
                            "BEGIN SELECT RAISE(ABORT,'simulated interruption'); END")
    with pytest.raises(sqlite3.IntegrityError, match="simulated interruption"):
        output.copy("Peptides")
    assert output.progress("Peptides")["rows"] == 0
    assert output.database.execute("SELECT count(*) FROM Peptides").fetchone()[0] == 0
    output.database.execute("DROP TRIGGER interruption")
    output.close()
    resumed = exporter(native, tmp_path / "export")
    assert resumed.copy("Peptides")["rows"] == 3
    resumed.close()


def test_logical_digest_is_independent_of_batch_sizes_and_restart(native, tmp_path):
    first = exporter(native, tmp_path / "one", batch_rows=1, lookup_keys=1)
    second = exporter(native, tmp_path / "two", batch_rows=500, lookup_keys=32)
    for table in first.tables:
        first.copy(table)
        second.copy(table)
    assert first.summary() == second.summary()
    assert first.summary()["selection_complete"]
    first.close()
    second.close()


def test_changed_source_selection_and_missing_columns_are_rejected(native, tmp_path):
    output = exporter(native, tmp_path / "export")
    output.copy("Peptides", max_batches=1)
    output.close()
    source, plan = native
    changed = copy.deepcopy(plan)
    changed["tables"][0]["columns"].remove("vector")
    with pytest.raises(ValueError, match="another source or selection"):
        exporter((source, changed), tmp_path / "export")
    changed["tables"][0]["columns"].append("absent")
    with pytest.raises(ValueError, match="column is absent"):
        exporter((source, changed), tmp_path / "other")
    source.execute("CREATE TABLE changed_schema(new_field)")
    with pytest.raises(ValueError, match="schema differs"):
        exporter(native, tmp_path / "changed")


def test_batch_and_sqlite_output_budgets_stop_without_false_completion(native, tmp_path):
    output = exporter(native, tmp_path / "small", batch_bytes=1)
    with pytest.raises(LimitExceeded, match="row exceeds"):
        output.copy("Peptides")
    assert not output.progress("Peptides")["complete"]
    output.close()
    source, plan = native
    source.executemany("INSERT INTO Peptides(id,vector) VALUES (?,?)", [(i, bytes(150000)) for i in range(20)])
    source.commit()
    output = exporter((source, plan), tmp_path / "limited", batch_rows=1)
    with pytest.raises(sqlite3.OperationalError, match="full"):
        output.copy("Peptides")
    state = output.progress("Peptides")
    assert not state["complete"] and 3 < state["rows"] < 23
    assert output.path.stat().st_size <= 1024 * 1024
    output.close()


def extension_plans(native):
    database, original = native
    database.executescript("CREATE TABLE ExtraMembership(workflow,spectrum); CREATE TABLE ExtraMetadata(name);")
    database.executemany("INSERT INTO ExtraMembership VALUES (?,?)", [(1, 3), (2, 99)])
    database.execute("INSERT INTO ExtraMetadata VALUES ('additional source definition')")
    database.commit()
    old = copy.deepcopy(original)
    old["source_schema_sha256"] = source_schema(database)[1]
    new = copy.deepcopy(old)
    new["tables"].insert(0, dict(table="ExtraMetadata", columns=["name"]))
    new["tables"].insert(-1, dict(table="ExtraMembership", columns=["workflow", "spectrum"]))
    new["tables"][-1]["select_from"].append(dict(table="ExtraMembership", columns=["workflow", "spectrum"]))
    return database, old, new


def test_explicit_additive_selection_keeps_checkpoints_and_records_previous_scope(native, tmp_path):
    database, old, new = extension_plans(native)
    output = exporter((database, old), tmp_path / "extend", batch_rows=1)
    previous = output.copy("Peptides", max_batches=1)
    output.copy("Membership")
    output.close()
    with pytest.raises(ValueError, match="another source or selection"):
        exporter((database, new), tmp_path / "extend")
    resumed = exporter((database, new), tmp_path / "extend", extend_selection=True)
    assert resumed.progress("Peptides") == previous
    history = json.loads(resumed.database.execute("SELECT value FROM _export_metadata WHERE key='selection_history'").fetchone()[0])
    assert len(history) == 1 and history[0]["previous_configuration"]["plan"] == old
    assert history[0]["previous_progress"]["Peptides"] == previous
    assert history[0]["added_tables"] == ["ExtraMetadata", "ExtraMembership"]
    assert history[0]["expanded_unstarted_selections"] == ["Acquisition"]
    for name in resumed.tables:
        resumed.copy(name)
    state = resumed.progress("Acquisition")
    assert state["rows"] == 4 and state["selected_keys"] == 4 and state["missing_keys"] == 1
    assert resumed.summary()["selection_complete"]
    resumed.close()
    again = exporter((database, new), tmp_path / "extend")
    assert again.summary()["selection_complete"]
    again.close()


@pytest.mark.parametrize("progress", ["queued", "started", "complete"])
def test_expanding_started_acquisitions_never_relabels_existing_scope(native, tmp_path, progress):
    database, old, new = extension_plans(native)
    output = exporter((database, old), tmp_path / "started", lookup_keys=1)
    output.copy("Membership")
    if progress == "queued":
        output._queue(output.tables["Acquisition"])
    else:
        output.copy("Acquisition", max_batches=1 if progress == "started" else None)
    previous = output.progress("Acquisition")
    output.close()
    with pytest.raises(ValueError, match="after it has started"):
        exporter((database, new), tmp_path / "started", extend_selection=True)
    reopened = exporter((database, old), tmp_path / "started")
    assert reopened.progress("Acquisition") == previous
    assert reopened.database.execute("SELECT value FROM _export_metadata WHERE key='selection_history'").fetchone() is None
    reopened.close()


def test_selection_extension_failure_is_atomic_and_cannot_remove_prior_facts(native, tmp_path):
    database, old, new = extension_plans(native)
    output = exporter((database, old), tmp_path / "atomic")
    previous = output.copy("Peptides")
    output.database.execute("CREATE TRIGGER fail_extension BEFORE INSERT ON _export_progress "
                            "WHEN NEW.name='ExtraMetadata' BEGIN SELECT RAISE(ABORT,'simulated extension failure'); END")
    output.close()
    with pytest.raises(sqlite3.IntegrityError, match="extension failure"):
        exporter((database, new), tmp_path / "atomic", extend_selection=True)
    reopened = exporter((database, old), tmp_path / "atomic")
    assert reopened.progress("Peptides") == previous
    assert reopened.database.execute("SELECT 1 FROM sqlite_master WHERE name='ExtraMetadata'").fetchone() is None
    assert reopened.database.execute("SELECT value FROM _export_metadata WHERE key='selection_history'").fetchone() is None
    reopened.database.execute("DROP TRIGGER fail_extension")
    reopened.close()
    removed_column = copy.deepcopy(new)
    removed_column["tables"][1]["columns"].remove("vector")
    with pytest.raises(ValueError, match="existing columns"):
        exporter((database, removed_column), tmp_path / "atomic", extend_selection=True)
    removed_table = copy.deepcopy(new)
    removed_table["tables"] = [t for t in removed_table["tables"] if t["table"] != "Peptides"]
    with pytest.raises(ValueError, match="remove or reorder"):
        exporter((database, removed_table), tmp_path / "atomic", extend_selection=True)
