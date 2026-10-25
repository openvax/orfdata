"""An inventory can name an entry at several sources; a ledger says what each one is."""

import hashlib

import pytest

from orfdata import viral_collection, viral_collections, viral_panel_coverage
from orfdata.record_sources import SourceRecord
from orfdata.reference_ledger import (
    COLUMNS,
    FRAME_COLUMNS,
    INCONSISTENT,
    INDEX,
    MEASURED,
    NUMERIC,
    PUBLISHED,
    UNAVAILABLE,
    Anchor,
    empty_ledger,
    ledger_from_rows,
    ledger_row,
    published,
    read_ledger_bytes,
    same_molecule,
    verify_record,
    write_ledger,
)

# Pinned in tools/build_hpv16_cds.py from GDC's published reference files, and
# reproduced exactly by the alternate reference the inventory names for HPV16.
GDC_HPV16_MD5 = "f2463721818c29998ef5ab917e6df74a"


def _row(sequence="ACGTACGT", **overrides):
    return dict(index=1, abbreviation="TV1", source="pave", identifier="TV1REF",
                state=PUBLISHED, length_nt=len(sequence),
                sequence_md5=hashlib.md5(sequence.encode()).hexdigest(),
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                coding_features=2, taxon_id=0, title="") | overrides


def _one(sequence="ACGTACGT", **overrides):
    """One ledger row, read back the way a packaged ledger hands it over."""
    row = _row(sequence, **overrides)
    return ledger_from_rows([row]).loc[row["abbreviation"]]


def _record(sequence="ACGTACGT", features=2):
    return SourceRecord(source="pave", identifier="TV1REF", sequence=sequence,
                        coding_features=tuple(range(features)))


def test_a_row_recognises_its_own_molecule():
    row = _one()
    assert published(row) and same_molecule(row, "ACGTACGT") and same_molecule(row, "acgtacgt")
    assert not same_molecule(row, "ACGTACGA")
    assert not same_molecule(_one(state=UNAVAILABLE), "ACGTACGT")
    assert not published(_one(state=INCONSISTENT)) and not published(None)


def test_a_fetched_record_is_checked_against_what_the_ledger_pins():
    row = _one()
    assert verify_record(row, _record()) is not None
    with pytest.raises(ValueError, match="no longer serves the molecule"):
        verify_record(row, _record("ACGTACGA"))
    with pytest.raises(ValueError, match="no longer serves the molecule"):
        verify_record(row, _record("ACGT"))
    with pytest.raises(ValueError, match="serves 3 coding features, not the 2"):
        verify_record(row, _record(features=3))
    with pytest.raises(ValueError, match="does not publish"):
        verify_record(_one(state=UNAVAILABLE), _record())


def test_a_ledger_round_trips_through_its_own_format():
    content = write_ledger(
        ledger_from_rows([_row(taxon_id=10001, title="Test virus 1"),
                          _row(index=2, abbreviation="TV2", state=UNAVAILABLE,
                               title="Test virus 2")]),
        panel="TEST", source="pave", source_url="https://example.invalid",
        checked="2026-01-01T00:00:00Z")
    lines = content.decode().splitlines()
    assert lines[0] == "# schema: orfdata.reference_ledger.v2"
    assert lines[4] == "\t".join(COLUMNS)

    ledger = read_ledger_bytes(content)
    assert set(ledger.index) == {"TV1", "TV2"}
    assert same_molecule(ledger.loc["TV1"], "ACGTACGT")
    assert ledger.at["TV1", "taxon_id"] == 10001
    # An entry the source does not publish pins nothing, whatever was measured.
    assert not published(ledger.loc["TV2"])
    assert ledger["length_nt"].isna()["TV2"] and ledger.at["TV2", "sequence_md5"] == ""
    assert ledger.at["TV2", "title"] == "Test virus 2"


def test_an_unreadable_ledger_is_refused_rather_than_half_read():
    content = write_ledger(ledger_from_rows([_row()]), panel="TEST", source="pave",
                           source_url="", checked="2026-01-01T00:00:00Z")
    with pytest.raises(ValueError, match="Unknown ledger state"):
        read_ledger_bytes(content.replace(b"\tpublished\t", b"\tmaybe\t"))
    doubled = content + content.decode().splitlines()[-1].encode() + b"\n"
    with pytest.raises(ValueError, match="Repeated ledger entry"):
        read_ledger_bytes(doubled)


def test_an_anchor_is_the_claim_a_ledger_rests_on():
    anchor = Anchor(source="pave", abbreviation="TV1",
                    sequence_md5=_row()["sequence_md5"], note="the molecule in use")
    assert anchor.check(ledger_from_rows([_row()])) is not None
    with pytest.raises(ValueError, match="no longer reproduces"):
        anchor.check(ledger_from_rows([_row("ACGTACGA")]))
    with pytest.raises(ValueError, match="no longer reproduces"):
        anchor.check(ledger_from_rows([]))


def test_the_published_collection_carries_a_ledger_for_each_of_its_sources():
    collection, = viral_collections()
    assert collection.sources() == ("genbank", "pave")

    primary = collection.ledger(collection.primary_source)
    inventory = viral_panel_coverage()
    assert set(primary.index) == set(inventory["abbreviation"])
    # Every row's accession resolves to the version the ledger pins.
    named = inventory.set_index("abbreviation")["genbank"].str.split(".").str[0]
    assert (primary["identifier"].str.split(".").str[0] == named).all()
    assert (primary["taxon_id"] > 0).all() and (primary["title"] != "").all()

    alternate = collection.ledger()
    with_pave = set(inventory.loc[inventory["pave_id"] != "", "abbreviation"])
    assert set(alternate.index) == with_pave and len(with_pave) == 189
    assert set(alternate["source"]) == {collection.alternate_source}
    # The anchor for the whole ledger: GDC's HPV16 contig is the alternate's molecule.
    assert alternate.at["HPV16", "sequence_md5"] == GDC_HPV16_MD5
    assert alternate.at["HPV16", "length_nt"] == 7906
    assert alternate.at["HPV16", "coding_features"] == 10
    for anchor in collection.anchors:
        anchor.check(collection.ledger(anchor.source))
    assert viral_collection("gdc-viruses") is collection


def test_every_ledger_carries_exactly_the_columns_it_documents():
    """Packaged, freshly built and empty ledgers are all the same frame."""
    collection, = viral_collections()
    built = ledger_from_rows([_row()])
    for ledger in (collection.ledger("genbank"), collection.ledger("pave"),
                   built, empty_ledger(), read_ledger_bytes(write_ledger(
                       built, panel="TEST", source="pave", source_url="",
                       checked="2026-01-01T00:00:00Z"))):
        assert ledger.index.name == INDEX
        assert tuple(ledger.columns) == FRAME_COLUMNS
    assert set(COLUMNS) == set(FRAME_COLUMNS) | {INDEX}
    # A row the source could not serve pins nothing in any measured column.
    pave = collection.ledger("pave")
    unserved = pave[pave["state"] != PUBLISHED]
    assert len(unserved) == 2
    for column in MEASURED:
        blank = unserved[column].isna() if column in NUMERIC else unserved[column] == ""
        assert blank.all(), column


def test_entries_that_did_not_check_out_are_named_as_such():
    ledger = viral_collections()[0].ledger()
    assert set(ledger["state"]) <= {PUBLISHED, UNAVAILABLE, INCONSISTENT}
    # An alternate whose own annotation does not reproduce its proteins is not used.
    ruled_out = ledger[ledger["state"].isin([UNAVAILABLE, INCONSISTENT])]
    assert not any(published(row) for _, row in ruled_out.iterrows())
    # A source with nothing packaged reads as an empty ledger, not a missing one.
    assert ledger_row(ledger, "NOT-IN-THE-PANEL") is None
