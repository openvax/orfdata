"""An inventory can name an entry at several sources; a ledger says what each one is."""

import hashlib

import pytest

from orfdata import viral_collection, viral_collections, viral_panel_coverage
from orfdata.record_sources import SourceRecord
from orfdata.reference_ledger import (
    COLUMNS,
    INCONSISTENT,
    PUBLISHED,
    UNAVAILABLE,
    LedgerAnchor,
    LedgerEntry,
    read_ledger_bytes,
    write_ledger,
)

# Pinned in tools/build_hpv16_cds.py from GDC's published reference files, and
# reproduced exactly by the alternate reference the inventory names for HPV16.
GDC_HPV16_MD5 = "f2463721818c29998ef5ab917e6df74a"


def _entry(sequence="ACGTACGT", **overrides):
    digest = hashlib.md5(sequence.encode()).hexdigest()
    fields = dict(index=1, abbreviation="TV1", source="pave", identifier="TV1REF",
                  state=PUBLISHED, length=len(sequence), sequence_md5=digest,
                  sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                  coding_features=2)
    return LedgerEntry(**dict(fields, **overrides))


def _record(sequence="ACGTACGT", features=2):
    return SourceRecord(source="pave", identifier="TV1REF", sequence=sequence,
                        coding_features=tuple(range(features)))


def test_an_entry_recognises_its_own_molecule():
    entry = _entry()
    assert entry.published and entry.matches("ACGTACGT") and entry.matches("acgtacgt")
    assert not entry.matches("ACGTACGA")
    assert not _entry(state=UNAVAILABLE).matches("ACGTACGT")
    assert not _entry(state=INCONSISTENT).published


def test_a_fetched_record_is_checked_against_what_the_ledger_pins():
    entry = _entry()
    assert entry.verify(_record()) is not None
    with pytest.raises(ValueError, match="no longer serves the molecule"):
        entry.verify(_record("ACGTACGA"))
    with pytest.raises(ValueError, match="no longer serves the molecule"):
        entry.verify(_record("ACGT"))
    with pytest.raises(ValueError, match="serves 3 coding features, not the 2"):
        entry.verify(_record(features=3))


def test_a_ledger_round_trips_through_its_own_format():
    content = write_ledger(
        [_entry(taxon_id=10001, title="Test virus 1"),
         _entry(index=2, abbreviation="TV2", state=UNAVAILABLE, title="Test virus 2")],
        panel="TEST", source="pave", source_url="https://example.invalid",
        checked="2026-01-01T00:00:00Z")
    lines = content.decode().splitlines()
    assert lines[0] == "# schema: orfdata.reference_ledger.v2"
    assert lines[4] == "\t".join(COLUMNS)

    entries = read_ledger_bytes(content)
    assert set(entries) == {"TV1", "TV2"}
    assert entries["TV1"].matches("ACGTACGT") and entries["TV1"].taxon_id == 10001
    # An entry the source does not publish pins nothing, whatever was measured.
    assert not entries["TV2"].published
    assert entries["TV2"].length == 0 and entries["TV2"].sequence_md5 == ""
    assert entries["TV2"].title == "Test virus 2"


def test_an_unreadable_ledger_is_refused_rather_than_half_read():
    content = write_ledger([_entry()], panel="TEST", source="pave", source_url="",
                           checked="2026-01-01T00:00:00Z")
    with pytest.raises(ValueError, match="Unknown ledger state"):
        read_ledger_bytes(content.replace(b"\tpublished\t", b"\tmaybe\t"))
    doubled = content + content.decode().splitlines()[-1].encode() + b"\n"
    with pytest.raises(ValueError, match="Repeated ledger entry"):
        read_ledger_bytes(doubled)


def test_an_anchor_is_the_claim_a_ledger_rests_on():
    anchor = LedgerAnchor(source="pave", abbreviation="TV1",
                          sequence_md5=_entry().sequence_md5, note="the molecule in use")
    assert anchor.check({"TV1": _entry()}) is not None
    with pytest.raises(ValueError, match="no longer reproduces"):
        anchor.check({"TV1": _entry("ACGTACGA")})
    with pytest.raises(ValueError, match="no longer reproduces"):
        anchor.check({})


def test_the_published_collection_carries_a_ledger_for_each_of_its_sources():
    collection, = viral_collections()
    assert collection.sources() == ("genbank", "pave")

    primary = collection.ledger(collection.primary_source)
    inventory = viral_panel_coverage()
    assert set(primary) == {entry.abbreviation for entry in inventory}
    # Every row's accession resolves to the version the ledger pins.
    assert all(primary[entry.abbreviation].identifier.split(".")[0] == entry.genbank.split(".")[0]
               for entry in inventory)
    assert all(entry.taxon_id and entry.title for entry in primary.values())

    alternate = collection.ledger()
    named = {entry.abbreviation for entry in inventory if entry.pave_id}
    assert set(alternate) == named and len(named) == 189
    assert {entry.source for entry in alternate.values()} == {collection.alternate_source}
    # The anchor for the whole ledger: GDC's HPV16 contig is the alternate's molecule.
    assert alternate["HPV16"].sequence_md5 == GDC_HPV16_MD5
    assert alternate["HPV16"].length == 7906 and alternate["HPV16"].coding_features == 10
    for anchor in collection.anchors:
        anchor.check(collection.ledger(anchor.source))
    assert viral_collection("gdc-viruses") is collection


def test_entries_that_did_not_check_out_are_named_as_such():
    ledger = viral_collections()[0].ledger()
    states = {entry.state for entry in ledger.values()}
    assert states <= {PUBLISHED, UNAVAILABLE, INCONSISTENT}
    # An alternate whose own annotation does not reproduce its proteins is not used.
    assert all(not entry.published for entry in ledger.values()
               if entry.state in (UNAVAILABLE, INCONSISTENT))
