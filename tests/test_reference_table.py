"""Every packaged reference table is read the same way, and only when it checks out."""

import hashlib

import pytest

from orfdata.reference_table import PackagedTable, _packaged, parse_table, write_table

COLUMNS = ("index", "abbreviation", "note")
SCHEMA = "orfdata.test_table.v1"
ROWS = ({"index": 1, "abbreviation": "TV1", "note": "first"},
        {"index": 2, "abbreviation": "TV2", "note": ""})


def _content(schema=SCHEMA, columns=COLUMNS):
    return write_table(ROWS, schema=schema, columns=columns,
                       metadata={"panel": "TEST", "checked_utc": "2026-01-01T00:00:00Z"})


def test_a_table_round_trips_through_its_own_format():
    metadata, rows = parse_table(_content(), schema=SCHEMA, columns=COLUMNS)
    assert metadata == {"schema": SCHEMA, "panel": "TEST", "checked_utc": "2026-01-01T00:00:00Z"}
    assert [row["abbreviation"] for row in rows] == ["TV1", "TV2"]
    assert rows[1]["note"] == ""


def test_a_table_that_is_not_what_it_claims_is_refused():
    with pytest.raises(ValueError, match="schema"):
        parse_table(_content(schema="orfdata.something_else.v1"), schema=SCHEMA, columns=COLUMNS)
    with pytest.raises(ValueError, match="columns"):
        parse_table(_content(), schema=SCHEMA, columns=COLUMNS + ("extra",))


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
    assert len(table.read()) == 2 and table.metadata()["panel"] == "TEST"
    assert table.to_dict()["columns"] == list(COLUMNS)

    _packaged.cache_clear()
    wrong = PackagedTable(resource="data/test/table.tsv", sha256="0" * 64,
                          schema=SCHEMA, columns=COLUMNS)
    with pytest.raises(ValueError, match="checksum mismatch"):
        wrong.read()
    _packaged.cache_clear()


def test_the_published_collections_declare_every_table_they_read():
    from orfdata import viral_collections

    for collection in viral_collections():
        assert collection.inventory.read()
        for source in collection.sources():
            table = collection.ledgers[source]
            assert table.schema and table.sha256 and table.read()
