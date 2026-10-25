"""How completely an installed collection is annotated, and what a second source would add."""


import pytest
from conftest import (
    VIRAL_CONTEXT,
    ledger_entry,
    serve_alternate_ledger,
    serve_records,
    viral_documents,
)

from orfdata import download_viral_references, viral_annotation_survey

UNANNOTATED = "ATG" + "GCT" * 5 + "TAA" + "ACGT" * 10
FIRST = "ATG" + "GCT" * 8 + "TAA" + "ACGT" * 10


def _document(sequence, features):
    return {
        "itemSequence": {"sequence": sequence},
        "features": [
            {"featureType": "CDS", "featureName": name, "proteinID": f"TVREF-{name}",
             "featureDescription": f"Test virus, {name}", "proteinSequence": protein,
             "featureLocations": [{"start": start, "end": end} for start, end in locations]}
            for name, protein, locations in features
        ],
    }


@pytest.fixture
def collection(viral_panel, monkeypatch):
    """An alternate source that publishes the same molecules, annotated further."""
    published = [
        # The same molecule as TV1's own record, but with a second, spliced model.
        ledger_entry("TV1", FIRST, index=1, features=2),
        # A molecule whose own record annotates nothing at all.
        ledger_entry("TV4", UNANNOTATED, index=4, features=1),
    ]
    return serve_alternate_ledger(monkeypatch, published)


def _served():
    return {
        "TV1REF": _document(FIRST, [
            ("E1", "MAAAAAAAA", ((1, 30),)),
            ("E1^E4", "MAA", ((1, 6), (25, 30))),
        ]),
        "TV4REF": _document(UNANNOTATED, [("E6", "MAAAAA", ((1, 21),))]),
    }


def test_the_survey_measures_each_genome_against_the_alternate_source(collection, viral_panel,
                                                                      tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents(), _served())
    download_viral_references(species="human", cache_root=str(tmp_path))
    survey = viral_annotation_survey(species="human", cache_root=str(tmp_path))
    rows = {row["abbreviation"]: row for row in survey["genomes"]}

    # TV1's own record annotates one of the two coding features the alternate has.
    assert rows["TV1"]["transcripts"] == 1
    assert rows["TV1"]["alternate_features"] == 2
    assert rows["TV1"]["shortfall"] == 1
    # TV4 annotated nothing and was filled from the alternate, so it is level.
    assert rows["TV4"]["annotation_source"] == "pave"
    assert rows["TV4"]["shortfall"] == 0
    # A row with no alternate has nothing to be measured against, and says so.
    assert rows["TV2"]["alternate_features"] is None and rows["TV2"]["shortfall"] is None

    totals = survey["totals"]
    assert totals["genomes"] == 4
    assert totals["genomes_behind_alternate"] == 1
    assert totals["coding_features_behind_alternate"] == 1
    # The survey carries the whole summary, so one call answers both questions.
    assert survey["genome_count"] == 4 and "layers" in survey


def test_the_survey_counts_partial_unnamed_and_spliced_annotation(collection, viral_panel,
                                                                   tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents(), _served())
    download_viral_references(species="human", cache_root=str(tmp_path))
    survey = viral_annotation_survey(species="human", cache_root=str(tmp_path))
    rows = {row["abbreviation"]: row for row in survey["genomes"]}
    # TV2's spliced CDS begins mid-codon and does not close its reading frame.
    assert rows["TV2"]["incomplete"] == 1 and rows["TV2"]["complete"] == 1
    assert rows["TV2"]["spliced"] == 1
    assert all(row["unnamed_genes"] == 0 for row in survey["genomes"])
    assert survey["totals"]["genomes_with_incomplete_transcripts"] == 1
    assert survey["totals"]["genomes_with_spliced_transcripts"] == 1


def test_a_selection_narrows_the_survey(collection, viral_panel, tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents(), _served())
    download_viral_references(species="human", cache_root=str(tmp_path))
    survey = viral_annotation_survey(species="human", cache_root=str(tmp_path),
                                     contexts=[VIRAL_CONTEXT])
    assert [row["abbreviation"] for row in survey["genomes"]] == ["TV1", "TV2"]
    assert survey["totals"]["genomes"] == 2
    assert survey["totals"]["genomes_behind_alternate"] == 1


def test_preferring_the_alternate_annotation_keeps_the_molecule(collection, viral_panel,
                                                                tmp_path, monkeypatch):
    """Taking a second source's gene models changes how the bases are read, not the bases."""
    serve_records(monkeypatch, viral_documents(), _served())
    kept = download_viral_references(species="human", cache_root=str(tmp_path))
    before = kept.genome_fasta.read_text()
    plain = viral_annotation_survey(species="human", cache_root=str(tmp_path))

    swapped = download_viral_references(species="human", cache_root=str(tmp_path / "wider"),
                                        prefer_alternate_annotation=True)
    assert swapped.genome_fasta.read_text() == before
    wider = viral_annotation_survey(species="human", cache_root=str(tmp_path / "wider"))
    rows = {row["abbreviation"]: row for row in wider["genomes"]}
    assert rows["TV1"]["annotation_source"] == "pave"
    assert rows["TV1"]["transcripts"] == 2 and rows["TV1"]["spliced"] == 1
    assert wider["totals"]["genomes_behind_alternate"] == 0
    # The rows the alternate does not improve on are left exactly as they were.
    assert rows["TV2"]["annotation_source"] == plain["genomes"][1]["annotation_source"]


def test_the_default_leaves_each_submission_its_own_annotation(collection, viral_panel,
                                                               tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents(), _served())
    download_viral_references(species="human", cache_root=str(tmp_path))
    survey = viral_annotation_survey(species="human", cache_root=str(tmp_path))
    rows = {row["abbreviation"]: row for row in survey["genomes"]}
    assert rows["TV1"]["annotation_source"] == "genbank"
