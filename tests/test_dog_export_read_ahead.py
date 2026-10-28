"""Indexed read-ahead changes bounded transport, not selected native facts."""

import json
import sqlite3
import sys

import pytest
from test_ftps_sqlite import LocalSource

from tools import export_dog_results
from tools.export_sqlite_rows import source_schema


def test_indexed_window_keeps_selection_and_digest_with_fewer_requests(tmp_path, monkeypatch):
    pytest.importorskip("apsw")
    source_path = tmp_path / "source.sqlite"
    with sqlite3.connect(source_path) as native:
        native.executescript("""
            PRAGMA page_size=4096;
            CREATE TABLE Membership(workflow,spectrum);
            CREATE TABLE Acquisition(workflow,spectrum,height,unused_spectrum);
            CREATE INDEX acquisition_key ON Acquisition(workflow,spectrum);
        """)
        rows = [(i // 200, 2**60 + i, None if i % 7 == 0 else float(i),
                 bytes([i % 251]) * (20000 if i % 37 == 0 else 100)) for i in range(800)]
        native.executemany("INSERT INTO Acquisition VALUES (?,?,?,?)", rows)
        selected = [row[:2] for row in rows[::3]]
        native.executemany("INSERT INTO Membership VALUES (?,?)", [*selected, selected[0], (77, 2**60 + 9999)])
        plan = dict(source_schema_sha256=source_schema(native)[1], tables=[
            dict(table="Membership", columns=["workflow", "spectrum"]),
            dict(table="Acquisition", columns=["workflow", "spectrum", "height"],
                 keys=["workflow", "spectrum"],
                 select_from=[dict(table="Membership", columns=["workflow", "spectrum"])])])
    data = source_path.read_bytes()
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr(export_dog_results, "PLAN_PATH", plan_path)
    monkeypatch.setattr(export_dog_results, "FTPSRanges", lambda identity, **kwargs: LocalSource(data, **kwargs))
    receipts = []
    for name, options in (("pages", []), ("window", ["--lookup-read-ahead-kib", "256"])):
        output = tmp_path / name
        monkeypatch.setattr(sys, "argv", ["export_dog_results", str(output), "--transfer-mib", "16",
                                         "--cache-mib", "4", "--scratch-mib", "16", "--output-mib", "1",
                                         "--read-ahead-kib", "4", "--lookup-keys", "32", *options])
        export_dog_results.main()
        receipt = json.loads((output / "export-receipt-0001.json").read_text())
        receipts.append(receipt)
        with sqlite3.connect(output / "native-evidence.partial.sqlite") as result:
            assert result.execute("SELECT workflow,spectrum,height FROM Acquisition ORDER BY workflow,spectrum").fetchall() == [
                row[:3] for row in rows[::3]]
            assert result.execute("SELECT source_key FROM _export_missing").fetchall() == [
                (json.dumps((77, 2**60 + 9999), separators=(",", ":")),)]
            assert "unused_spectrum" not in {r[1] for r in result.execute("PRAGMA table_info(Acquisition)")}
        assert receipt["failure"] is None and receipt["result"]["selection_complete"]
        assert receipt["acquisition"]["cached_payload_bytes"] <= len(data)
        assert receipt["acquisition"]["payload_bytes"] <= 16 * 1024**2
        assert receipt["export_options"]["read_ahead_kib"] == 4
    assert receipts[0]["result"] == receipts[1]["result"]
    assert receipts[0]["export_options"]["lookup_read_ahead_kib"] == 4
    assert receipts[1]["export_options"]["lookup_read_ahead_kib"] == 256
    assert receipts[1]["acquisition"]["range_requests"] < receipts[0]["acquisition"]["range_requests"]


@pytest.mark.parametrize("option", ["--read-ahead-kib", "--lookup-read-ahead-kib"])
@pytest.mark.parametrize("invalid", ["0", "-4", "6", "16388"])
def test_invalid_window_is_rejected_before_creating_cache(tmp_path, monkeypatch, option, invalid):
    output = tmp_path / "untouched"
    monkeypatch.setattr(sys, "argv", ["export_dog_results", str(output), option, invalid])
    with pytest.raises(SystemExit) as error:
        export_dog_results.main()
    assert error.value.code == 2
    assert not output.exists()
