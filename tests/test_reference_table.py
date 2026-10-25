"""Every packaged reference table is read the same way, and only when it checks out."""

import hashlib

import pandas as pd
import pytest

from orfdata.reference_table import PackagedTable, _packaged, read_table, write_table

COLUMNS = ("index", "abbreviation", "note")
SCHEMA = "orfdata.test_table.v1"
ROWS = pd.DataFrame([{"index": 1, "abbreviation": "TV1", "note": "first"},
                     {"index": 2, "abbreviation": "TV2", "note": ""}])


def _content(schema=SCHEMA, columns=COLUMNS):
    return write_table(ROWS, schema=schema, columns=columns,
                       metadata={"panel": "TEST", "checked_utc": "2026-01-01T00:00:00Z"})


def test_a_table_round_trips_through_its_own_format():
    frame = read_table(_content(), schema=SCHEMA, columns=COLUMNS, numeric=("index",))
    assert frame.attrs == {"schema": SCHEMA, "panel": "TEST",
                           "checked_utc": "2026-01-01T00:00:00Z"}
    assert frame["abbreviation"].to_list() == ["TV1", "TV2"]
    # A blank stays a blank, and a declared numeric column reads as a number.
    assert frame.at[1, "note"] == "" and frame["index"].to_list() == [1, 2]
    assert str(frame["index"].dtype) == "Int64"


def test_a_table_that_is_not_what_it_claims_is_refused():
    with pytest.raises(ValueError, match="schema"):
        read_table(_content(schema="orfdata.something_else.v1"), schema=SCHEMA, columns=COLUMNS)
    with pytest.raises(ValueError, match="columns"):
        read_table(_content(), schema=SCHEMA, columns=COLUMNS + ("extra",))
    with pytest.raises(ValueError, match="Unreadable .* 'abbreviation'"):
        read_table(_content(), schema=SCHEMA, columns=COLUMNS, numeric=("abbreviation",))


def _packaged_as(monkeypatch, tmp_path, content):
    resource = tmp_path / "table.tsv"
    resource.write_bytes(content)

    class _Anchor:
        def joinpath(self, _name):
            return resource

    monkeypatch.setattr("orfdata.reference_table.files", lambda _package: _Anchor())
    _packaged.cache_clear()


def test_a_packaged_table_is_read_only_when_its_bytes_are_the_pinned_ones(tmp_path, monkeypatch):
    content = _content()
    _packaged_as(monkeypatch, tmp_path, content)
    table = PackagedTable(resource="data/test/table.tsv",
                          sha256=hashlib.sha256(content).hexdigest(),
                          schema=SCHEMA, columns=COLUMNS)
    assert table.name == "table.tsv"
    assert len(table.read()) == 2 and table.read().attrs["panel"] == "TEST"
    assert table.to_dict()["columns"] == list(COLUMNS)
    # Each caller owns what it reads, so one caller's edits are nobody else's.
    mine = table.read()
    mine.loc[0, "note"] = "scribbled on"
    assert table.read().at[0, "note"] == "first"

    _packaged.cache_clear()
    wrong = PackagedTable(resource="data/test/table.tsv", sha256="0" * 64,
                          schema=SCHEMA, columns=COLUMNS)
    with pytest.raises(ValueError, match="checksum mismatch"):
        wrong.read()
    _packaged.cache_clear()


def test_the_published_collections_declare_every_table_they_read():
    from orfdata import viral_collections

    for collection in viral_collections():
        assert len(collection.inventory.read())
        for source in collection.sources():
            table = collection.ledgers[source]
            assert table.schema and table.sha256 and len(table.read())
