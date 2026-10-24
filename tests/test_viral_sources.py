"""A collection carries the molecule its reference was built from, and says where each came from."""

import hashlib
import json

import pytest
from conftest import serve_alternate_ledger, serve_records, viral_documents

from orfdata import download_viral_references, viral_reference_summary
from orfdata.reference_ledger import INCONSISTENT, PUBLISHED, UNAVAILABLE, LedgerEntry

# The fourth synthetic record annotates nothing; the first is a different molecule
# from the alternate reference its inventory row names.
UNANNOTATED = "ATG" + "GCT" * 5 + "TAA" + "ACGT" * 10
ALTERNATE = "ATG" + "GCT" * 8 + "TAA" + "ACGT" * 12
PROTEIN = "MAAAAA"


def _document(sequence, protein=PROTEIN, locations=((1, 21),), name="E6"):
    return {
        "itemSequence": {"sequence": sequence},
        "features": [
            {"featureType": "URR", "featureName": "URR",
             "featureLocations": [{"start": 1, "end": 12}]},
            {"featureType": "CDS", "featureName": name, "proteinID": f"TVREF-{name}",
             "featureDescription": f"Test virus, {name} gene", "proteinSequence": protein,
             "featureLocations": [{"start": start, "end": end} for start, end in locations]},
        ],
    }


def _ledger_entry(abbreviation, sequence, *, index=1, state=PUBLISHED, features=1):
    return LedgerEntry(
        index=index, abbreviation=abbreviation, source="pave", identifier=f"{abbreviation}REF",
        state=state, length=len(sequence),
        sequence_md5=hashlib.md5(sequence.encode()).hexdigest(),
        sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(), coding_features=features)


@pytest.fixture
def collection(viral_panel, monkeypatch):
    """The synthetic panel, with an alternate reference named for two of its rows."""
    published = {
        # TV4's record annotates nothing, and the alternate is the same molecule.
        "TV4": _ledger_entry("TV4", UNANNOTATED, index=4),
        # TV1's record is a different molecule from the alternate its row names.
        "TV1": _ledger_entry("TV1", ALTERNATE, index=1),
    }
    serve_alternate_ledger(monkeypatch, published)
    return published


def _rows(table):
    lines = table.read_text().splitlines()
    columns = lines[0].split("\t")
    return {row["abbreviation"]: row for row in
            (dict(zip(columns, line.split("\t"))) for line in lines[1:])}


def test_the_molecule_the_reference_was_built_from_is_the_one_carried(collection, viral_panel,
                                                                      tmp_path, monkeypatch):
    serve_records(monkeypatch, _genbank(viral_panel), {
        "TV1REF": _document(ALTERNATE, protein="MAAAAAAAA", locations=((1, 30),)),
        "TV4REF": _document(UNANNOTATED),
    })
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    rows = _rows(bundle.genome_table)

    # TV1's row names an alternate that is a different molecule, so that one is carried.
    assert rows["TV1"]["genome_source"] == "pave"
    assert rows["TV1"]["genome_id"] == "TV1REF"
    assert rows["TV1"]["alternate_sequence"] == "same"
    assert rows["TV1"]["accession_version"] == "VTEST001.1"  # The inventory's own accession.
    # TV4's alternate is the same molecule, so its own record is kept and annotated.
    assert rows["TV4"]["genome_source"] == "genbank"
    assert rows["TV4"]["annotation_source"] == "pave"
    # A row with no alternate named is untouched.
    assert rows["TV2"]["genome_source"] == "genbank" and rows["TV2"]["alternate_source"] == ""

    summary = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert summary["genome_sources"] == {"genbank": 3, "pave": 1}
    assert summary["genomes_unlike_alternate"] == []


def test_keeping_the_inventory_record_is_a_choice_that_is_recorded(collection, viral_panel,
                                                                   tmp_path, monkeypatch):
    serve_records(monkeypatch, _genbank(viral_panel), {"TV4REF": _document(UNANNOTATED)})
    bundle = download_viral_references(species="human", cache_root=str(tmp_path),
                                       prefer_reference_molecule=False)
    rows = _rows(bundle.genome_table)
    assert rows["TV1"]["genome_source"] == "genbank" and rows["TV1"]["genome_id"] == "VTEST001.1"
    assert rows["TV1"]["alternate_sequence"] == "differs"
    summary = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert summary["genomes_unlike_alternate"] == ["TV1"]


def test_an_alternate_that_cannot_be_had_leaves_the_download_standing(collection, viral_panel,
                                                                      tmp_path, monkeypatch):
    # Nothing is served for either alternate: the 200-record download is not lost.
    serve_records(monkeypatch, _genbank(viral_panel))
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    rows = _rows(bundle.genome_table)
    assert rows["TV1"]["genome_source"] == "genbank"
    assert rows["TV1"]["alternate_sequence"] == "differs"
    assert rows["TV4"]["transcript_count"] == "0"
    notes = json.loads((bundle.path / "metadata.json").read_text())["source_notes"]
    assert "TV1" in notes and "TV4" in notes
    assert "404" in notes["TV1"] or "refused" in notes["TV1"]

    summary = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert summary["genomes_unlike_alternate"] == ["TV1"]
    assert summary["genomes_without_transcripts"] == ["TV4"]


def test_an_alternate_that_does_not_match_the_pinned_ledger_is_refused(collection, viral_panel,
                                                                       tmp_path, monkeypatch):
    from orfdata.viral_annotation import translate

    # A self-consistent document, but not the molecule the ledger pins.
    served = "ACGT" * 20
    serve_records(monkeypatch, _genbank(viral_panel), {
        "TV1REF": _document(served, protein=translate(served[:9], table=1),
                            locations=((1, 9),)),
        "TV4REF": _document(UNANNOTATED),
    })
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    rows = _rows(bundle.genome_table)
    assert rows["TV1"]["genome_source"] == "genbank"
    notes = json.loads((bundle.path / "metadata.json").read_text())["source_notes"]
    assert "no longer serves the molecule" in notes["TV1"]


@pytest.mark.parametrize("state", [UNAVAILABLE, INCONSISTENT])
def test_an_alternate_the_ledger_rules_out_is_never_fetched(state, viral_panel, tmp_path,
                                                            monkeypatch):
    entries = {"TV1": _ledger_entry("TV1", ALTERNATE, state=state)}
    serve_alternate_ledger(monkeypatch, entries)
    asked = serve_records(monkeypatch, _genbank(viral_panel))
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    assert not any("pave" in url for url in asked)
    rows = _rows(bundle.genome_table)
    assert rows["TV1"]["genome_source"] == "genbank"
    assert rows["TV1"]["alternate_sequence"] == "unavailable"


def _genbank(entries):
    """The synthetic GenBank documents the collection's primary source serves."""
    return viral_documents()


def test_sources_are_a_published_registry():
    from orfdata import reference_source, reference_sources

    sources = reference_sources()
    assert set(sources) == {"genbank", "pave"}
    assert all(source.name and source.description for source in sources.values())
    assert reference_source("pave") is sources["pave"]
    with pytest.raises(ValueError, match="Unknown reference source"):
        reference_source("ensembl")


def test_collections_are_a_published_registry():
    from orfdata import viral_collection, viral_collections

    collection, = viral_collections()
    assert collection.name == "gdc-viruses"
    assert collection.to_dict()["target_taxon_id"] == 9606
    assert collection.source().name == collection.primary_source
    assert collection.source(collection.alternate_source).name == "pave"
    assert viral_collection("GDC_VIRUSES") is collection  # Names normalize.
    with pytest.raises(ValueError, match="Unknown viral collection"):
        viral_collection("gdc-bacteria")
    with pytest.raises(ValueError, match="is a different target"):
        viral_collection("gdc-viruses", species="mouse")


def test_a_source_record_names_its_own_document():
    from orfdata.record_sources import SourceRecord

    record = SourceRecord(source="pave", identifier="TV1REF", sequence="ACGT",
                          coding_features=(), document=b"{}", suffix=".json")
    assert record.document_name("TV 1") == "TV_1__pave__TV1REF.json.gz"
    assert record.length == 4
    md5, sha256 = record.checksums()
    assert len(md5) == 32 and len(sha256) == 64
