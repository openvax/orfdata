"""Source membership and complete sequence compatibility remain distinct facts."""

import base64
import copy
import json
import sqlite3
from pathlib import Path

import pytest
from test_pfizer_vendor import decoder, fixture_search

from orfdata.sources.pfizer import protein_manifest
from orfdata.sources.pfizer_attribution import CanineEvidence, build_attribution
from orfdata.sources.pfizer_evidence import EvidenceWriter

DIRECTORY = Path(__file__).parent / "data/pfizer/vendor"
MEMBERSHIPS = json.loads((DIRECTORY / "memberships.json").read_text())


def component(path, tables):
    """Explicit local regression selection; rowids below are fixture-only IDs."""
    tables = {name: [{key: base64.b64decode(value["base64"]) if isinstance(value, dict) and value.get("type") == "blob"
                     else value for key, value in row.items()} for row in rows] for name, rows in tables.items()}
    plan = dict(scope="Explicit native/synthetic regression slice; no source completeness claim", tables=[
        dict(table=name, columns=list(rows[0]), role="decoy" if name.startswith("Decoy") else "target")
        for name, rows in tables.items()])
    db = sqlite3.connect(path)
    db.execute("PRAGMA foreign_keys=ON")
    writer = EvidenceWriter(db, plan, decoder())
    for name, rows in tables.items():
        for index, row in enumerate(rows, 1):
            writer.insert(name, dict(row, _source_rowid=index))
    db.execute("UPDATE evidence_metadata SET value='\"complete\"' WHERE name='state'")
    db.commit()
    return db


@pytest.fixture
def native(tmp_path):
    path = tmp_path / "native.sqlite"
    searched = fixture_search(DIRECTORY / "membership-search-slice.fasta", MEMBERSHIPS["fixture_fasta_records"])
    with component(path, MEMBERSHIPS["tables"]) as db:
        counts = build_attribution(db, searched, max_residues=100)
    return CanineEvidence(path), counts


def test_native_memberships_preserve_ids_positions_and_group_resolution(native):
    evidence, counts = native
    assert counts["native_links"] == 6 and counts["unresolved_links"] == 0
    assert counts["target_proteins"] == counts["exact_native_proteins"] == 6
    assert counts["unmatched_native_proteins"] == 0
    group = list(evidence.linked("TargetProteinGroups", (1,), via="TargetProteinGroupsTargetProteins"))
    assert len(group) == 1 and group[0]["status"] == "resolved"
    original = MEMBERSHIPS["tables"]["TargetProteinGroupsTargetProteins"][0]
    assert group[0]["source_link"]["record"] == dict(original, _source_rowid=1)
    assert group[0]["left"]["record"]["ProteinGroupID"] == 1
    assert group[0]["right"]["record"]["UniqueSequenceID"] == -9222281927842617047
    assert list(evidence.linked("TargetProteins", (-9222281927842617047,),
                                via="TargetProteinGroupsTargetProteins", direction="incoming")) == group
    peptide = list(evidence.linked("TargetPeptideGroups", (1,), via="TargetPeptideGroupsTargetProteins"))
    assert len(peptide) == 1
    assert peptide[0]["source_link"]["record"]["PositionsinProteins"] == "[7-59]"
    assert peptide[0]["source_link"]["record"]["QuanUsage"] == 1
    native_protein = peptide[0]["right"]["record"]
    match = evidence.protein_search_match("TargetProteins", (native_protein["UniqueSequenceID"],))
    assert match["status"] == "exact_sequence" and match["accession"] == native_protein["Accession"]
    assert evidence.search_protein(match["accession"]).sequence == native_protein["Sequence"]
    assert evidence.search_protein("unknown") is None


def test_reported_parent_list_and_full_fasta_candidates_do_not_collapse(native):
    evidence, counts = native
    shared, = evidence.records("TargetPeptideGroups", key=(8,))
    reported = {value.strip() for value in shared["ParentProteinAccessions"].split(";")}
    assert reported == {"A0A5F4BP72", "E2RMY1"}
    candidates = list(evidence.candidates(shared["Sequence"]))
    assert {hit["accession"] for hit in candidates} == reported | {"E2RPC3", "J9NWC2"}
    expected = MEMBERSHIPS["candidate_placements"]
    observed = [hit for sequence in sorted({hit["peptide"] for hit in expected}) for hit in evidence.candidates(sequence)]
    assert observed == expected
    assert counts["peptide_sequences"] == 4 and counts["exact_placements"] == 7
    assert counts["no_search_match"] == 0
    assert list(evidence.candidates("UNSEARCHEDEXAMPLE"))[0]["status"] == "not_searched"


def synthetic_tables():
    return dict(
        TargetProteins=[dict(UniqueSequenceID=-9223372036854775808, Accession="ONE", Sequence="MIIILLLAAA"),
                        dict(UniqueSequenceID=9, Accession="TWO", Sequence="MISMATCH"),
                        dict(UniqueSequenceID=10, Accession="ABSENT", Sequence="AAAA")],
        DecoyProteins=[dict(UniqueSequenceID=1, Accession="ONE", Sequence="AAALLLIIIM")],
        TargetPeptideGroups=[dict(PeptideGroupID=i, Sequence=s) for i, s in enumerate(("II", "AA", "ZZZ"), 1)],
        TargetPsms=[dict(WorkflowID=-160, PeptideID=1, Sequence="II")],
        TargetPeptideGroupsTargetProteins=[
            dict(TargetPeptideGroupsPeptideGroupID=1, TargetProteinsUniqueSequenceID=-9223372036854775808),
            dict(TargetPeptideGroupsPeptideGroupID=1, TargetProteinsUniqueSequenceID=999999)],
        DecoyPeptideGroups=[dict(PeptideGroupID=1, Sequence="MIIILLLAAA")],
        DecoyPeptideGroupsDecoyProteins=[dict(DecoyPeptideGroupsPeptideGroupID=1, DecoyProteinsUniqueSequenceID=1)],
        QuanSpectrumInfo=[dict(WorkflowID=-160, SpectrumID=8546263)],
        QuanSpectrumInfoTargetPeptideGroups=[dict(QuanSpectrumInfoWorkflowID=-160, QuanSpectrumInfoSpectrumID=8546263,
                                                 TargetPeptideGroupsPeptideGroupID=1, QuanUsage=1)])


def test_repeated_il_candidates_missing_matches_dangling_links_and_decoys(tmp_path):
    fasta = tmp_path / "search.fasta"
    fasta.write_text(">sp|ONE|ONE_DOG synthetic OX=9615\nMIIILLLAAA\n>tr|TWO|TWO_DOG synthetic OX=9615\nMAAAA\n")
    searched = fixture_search(fasta, 2)
    digests = []
    for size in (3, 250000):
        path = tmp_path / f"batch-{size}.sqlite"
        with component(path, synthetic_tables()) as db:
            counts = build_attribution(db, searched, max_residues=size)
        evidence = CanineEvidence(path)
        assert counts["peptide_sequences"] == 3  # Duplicate PSM sequence searched once.
        assert counts["no_search_match"] == 1 and counts["unresolved_links"] == 1
        assert counts["exact_native_proteins"] == 1 and counts["unmatched_native_proteins"] == 2
        assert counts["decoy_proteins"] == 1
        ii = list(evidence.candidates("II"))
        assert [(r["start"], r["end"]) for r in ii] == [(1, 3), (2, 4), (3, 5), (4, 6), (5, 7)]
        assert [r["basis"] for r in ii] == ["exact_residues"] * 2 + ["isoleucine_leucine_equivalent"] * 3
        assert len(list(evidence.candidates("AA"))) == 5
        assert list(evidence.candidates("ZZZ"))[0]["status"] == "no_search_match"
        assert list(evidence.candidates("MIIILLLAAA"))[0]["status"] == "not_searched"  # Decoy-only sequence.
        links = list(evidence.linked("TargetPeptideGroups", (1,)))
        assert [link["status"] for link in links] == ["resolved", "unresolved"]
        assert links[1]["right"] is None
        assert links[1]["source_link"]["record"]["TargetProteinsUniqueSequenceID"] == 999999
        assert evidence.protein_search_match("TargetProteins", (9,))["status"] == "sequence_mismatch"
        assert evidence.protein_search_match("TargetProteins", (10,))["status"] == "not_in_search"
        assert evidence.protein_search_match("DecoyProteins", (1,))["status"] == "decoy_partition"
        decoy = next(evidence.linked("DecoyPeptideGroups", (1,)))
        assert decoy["right"]["role"] == "decoy"
        quant = next(evidence.linked("TargetPeptideGroups", (1,), via="QuanSpectrumInfoTargetPeptideGroups", direction="incoming"))
        assert quant["source_link"]["record"]["QuanUsage"] == 1
        assert quant["left"]["record"]["SpectrumID"] == 8546263
        digests.append(evidence.metadata["candidate_digest"])
    assert digests[0] == digests[1]  # Resource batching cannot change the logical hash.


def test_invalid_or_incomplete_attribution_is_not_exposed_as_empty_evidence(tmp_path):
    path = tmp_path / "incomplete.sqlite"
    searched = fixture_search(DIRECTORY / "membership-search-slice.fasta", MEMBERSHIPS["fixture_fasta_records"])
    tables = copy.deepcopy(MEMBERSHIPS["tables"])
    tables["TargetProteins"].append(tables["TargetProteins"][0])
    with component(path, tables) as db:
        with pytest.raises(ValueError, match="not unique"):
            build_attribution(db, searched)
    with pytest.raises(ValueError, match="incomplete"):
        CanineEvidence(path)


def test_link_filters_are_explicit(native):
    evidence, _ = native
    assert list(evidence.linked("TargetProteinGroups", (999999,))) == []
    with pytest.raises(ValueError, match="direction"):
        list(evidence.linked("TargetProteinGroups", (1,), direction="both"))
    with pytest.raises(ValueError, match="relationship"):
        list(evidence.linked("TargetProteinGroups", (1,), via="inferred protein detection"))
    with pytest.raises(ValueError, match="protein record"):
        evidence.protein_search_match("TargetPeptideGroups", (1,))


def test_attribution_cannot_join_a_different_dataset_reference(tmp_path):
    searched = fixture_search(DIRECTORY / "membership-search-slice.fasta", MEMBERSHIPS["fixture_fasta_records"])
    with component(tmp_path / "wrong-reference.sqlite", MEMBERSHIPS["tables"]) as db:
        db.execute("CREATE TABLE metadata(name PRIMARY KEY,value)")
        db.execute("INSERT INTO metadata VALUES ('manifest',?)", (json.dumps(protein_manifest("dog").to_dict()),))
        with pytest.raises(ValueError, match="differs from the quantitative"):
            build_attribution(db, searched)
