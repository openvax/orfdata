"""How a coding sequence begins, and what the source says the ribosome did to read it."""

import csv

import pytest
from conftest import Unquoted, serve_records, viral_documents

from orfdata import download_viral_references, viral_reference_summary
from orfdata.genbank import parse
from orfdata.viral_annotation import (
    ALTERNATIVE_START,
    ATG_START,
    DECLARED_START,
    NON_CANONICAL_START,
    PARTIAL_START,
    UNKNOWN_START,
    coding_features,
    declares_initiator,
    initiation,
    translate,
)


@pytest.mark.parametrize("codon, residue", [
    ("ATY", "I"),   # ATC and ATT are both isoleucine.
    ("TGY", "C"),   # TGC and TGT are both cysteine.
    ("TTR", "L"),   # TTA and TTG are both leucine.
    ("CTN", "L"),   # Every CTN codon is leucine.
    ("GGN", "G"),
])
def test_an_ambiguous_codon_that_cannot_change_the_residue_still_translates(codon, residue):
    """A published sequence with one ambiguity code still has a determined protein."""
    assert translate(codon, table=1) == residue


@pytest.mark.parametrize("codon", ["RAT", "MGN", "NNN", "AT?"])
def test_a_codon_whose_residue_really_is_ambiguous_reads_as_X(codon):
    assert translate(codon, table=1) == "X"


@pytest.mark.parametrize("codon, table, declared, expected", [
    ("ATG", 1, False, ATG_START),
    ("CTG", 1, False, ALTERNATIVE_START),   # An initiation codon of the standard code.
    ("TTG", 1, False, ALTERNATIVE_START),
    ("AAT", 1, False, NON_CANONICAL_START),  # Not an initiation codon in any sense.
    ("ACG", 1, True, DECLARED_START),        # The source declares the initiator itself.
    ("CTG", 11, False, UNKNOWN_START),       # A genetic code ORFdata does not carry.
    ("ATG", 11, False, ATG_START),           # ATG needs no table to be read.
])
def test_initiation_names_what_the_source_recorded(codon, table, declared, expected):
    assert initiation(codon, table=table, declared=declared) == expected


def test_a_truncated_five_prime_end_has_no_initiation_codon_to_classify():
    assert initiation("CTG", table=1, partial_start=True) == PARTIAL_START
    assert initiation("AT", table=1) == PARTIAL_START


def test_transl_except_counts_only_over_the_first_codon():
    """The same qualifier marks a read-through residue, which is not an initiator."""
    exons = ((100, 150, "+"),)
    assert declares_initiator("(pos:100..102,aa:Met)", exons, 1)
    assert not declares_initiator("(pos:130..132,aa:Sec)", exons, 1)
    assert declares_initiator("(pos:101..103,aa:Met)", exons, 2)  # codon_start shifts it.
    assert not declares_initiator("", exons, 1)
    assert not declares_initiator("(pos:1..3,aa:Met)", exons, 1)
    # On the minus strand the first codon is at the high end of the block.
    assert declares_initiator("(pos:complement(148..150),aa:Met)", ((100, 150, "-"),), 1)


def test_an_unreadable_exception_is_not_mistaken_for_a_declaration():
    assert not declares_initiator("(pos:one..three,aa:Met)", ((1, 30, "+"),), 1)
    assert not declares_initiator("Met", ((1, 30, "+"),), 1)


ALTERNATIVE_RECORD = {
    # An ACG start the source declares, with ribosome-profiling evidence, as CMV's UL30A does.
    "declared": ("1..30", {"gene": "UL30A", "product": "protein UL30A", "codon_start": 1,
                           "transl_except": Unquoted("(pos:1..3,aa:Met)"),
                           "translation": "MAAAAAAAA"}),
    # A CTG start, which the standard code lists as an initiation codon.
    "alternative": ("1..30", {"gene": "E4", "product": "protein E4", "codon_start": 1,
                              "translation": "MAAAAAAAA"}),
}


def _record(location, qualifiers, sequence):
    from conftest import _genbank

    return parse(_genbank("VTEST009.1", sequence, [(location, qualifiers)], name="VTEST009"))


def test_a_declared_start_and_an_alternative_one_are_told_apart():
    sequence = "ACG" + "GCT" * 8 + "TAA" + "ACGT" * 5
    location, qualifiers = ALTERNATIVE_RECORD["declared"]
    feature, = coding_features(_record(location, qualifiers, sequence))
    assert (feature["start_codon"], feature["initiation"]) == ("ACG", DECLARED_START)

    location, qualifiers = ALTERNATIVE_RECORD["alternative"]
    feature, = coding_features(_record(location, qualifiers, "CTG" + sequence[3:]))
    assert (feature["start_codon"], feature["initiation"]) == ("CTG", ALTERNATIVE_START)


def test_a_valueless_qualifier_records_the_recoding_it_claims():
    """HTLV-1 writes its frameshift products as a join that reads one base twice."""
    sequence = "ATG" + "GCT" * 8 + "TAA" + "ACGT" * 10
    feature, = coding_features(_record(
        "join(1..12,12..27)", {"gene": "gag-pro", "product": "Pr gag-pro",
                               "ribosomal_slippage": None, "codon_start": 1,
                               "translation": "MAAAXAAAA"}, sequence))
    assert feature["recoding"] == "ribosomal_slippage"
    assert feature["start_codon"] == "ATG"


def test_a_cds_with_nothing_unusual_says_so():
    sequence = "ATG" + "GCT" * 8 + "TAA" + "ACGT" * 5
    feature, = coding_features(_record(
        "1..30", {"gene": "E1", "product": "protein E1", "codon_start": 1,
                  "translation": "MAAAAAAAA"}, sequence))
    assert (feature["start_codon"], feature["initiation"], feature["recoding"]) == (
        "ATG", ATG_START, "")


def test_the_transcript_table_records_how_each_transcript_begins(viral_panel, tmp_path,
                                                                 monkeypatch):
    serve_records(monkeypatch, viral_documents())
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    rows = list(csv.DictReader(open(bundle.transcript_table), delimiter="\t"))
    assert {"start_codon", "initiation", "recoding"} <= set(rows[0])
    by_id = {row["transcript_id"]: row for row in rows}
    whole = by_id["VTEST001.1_cds1"]
    assert (whole["start_codon"], whole["initiation"]) == ("ATG", ATG_START)
    # The second record's spliced CDS begins one base in and publishes an
    # M-leading protein over a codon that does not initiate anything — the
    # shape HPV63 and HPV65 have, where the submitter says as much in a note.
    odd = by_id["VTEST002.1_cds1"]
    assert odd["codon_start"] == "2" and odd["initiation"] == NON_CANONICAL_START

    summary = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert summary["initiations"][ATG_START] >= 1
    assert summary["non_canonical_starts"] == ["VTEST002.1_cds1", "VTEST002.1_cds2"]
    assert summary["recoded_transcripts"] == {}
