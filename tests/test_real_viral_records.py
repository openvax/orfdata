"""Real GenBank and PaVE records, through the same download and binding path as a collection.

Synthetic records only hold the cases their author thought of. Each record here
is a real source record (three trimmed to the features that matter) that once
hit a defect the synthetic suite missed: what a defect looked like on real data
is pinned by the expectations below. The records live in ``tests/data/viral_records``.
"""

import gzip
from pathlib import Path

import pandas as pd
import pytest
from conftest import VIRAL_CONTEXT, serve_records

from orfdata import download_viral_references, viral_reference, viral_references
from orfdata import viral_panel as viral_panel_module
from orfdata.record_sources import GENBANK_SOURCE, reference_source
from orfdata.viral_annotation import name_transcripts, unique_source_proteins

DATA = Path(__file__).parent / "data" / "viral_records"

# abbreviation, accession: every one a record of the GDC collection.
PANEL = (("HBV", "X04615.1"), ("HTLV-1", "AF033817.1"), ("HPV63", "X70828.1"),
         ("EBV", "AJ507799.2"), ("HPV171", "KF006398.1"), ("HPV175", "KC108721.1"))


def _record(name):
    return gzip.open(DATA / f"{name}.gb.gz", "rt").read()


@pytest.fixture(scope="module")
def collection(tmp_path_factory):
    """Download the real records as a collection and return the cache root."""
    patch = pytest.MonkeyPatch()
    inventory = pd.DataFrame([
        {"index": index, "name": abbreviation, "abbreviation": abbreviation,
         "refseq": "", "genbank": accession, "pave_id": "", "pave_notes": "",
         "bundled_reference_id": "", "coding_models": pd.NA, "status": "not_bundled",
         "available": False, "contexts": (VIRAL_CONTEXT,)}
        for index, (abbreviation, accession) in enumerate(PANEL, 1)])
    inventory["coding_models"] = inventory["coding_models"].astype("Int64")
    patch.setattr(viral_panel_module, "viral_panel_coverage", lambda **options: inventory)
    serve_records(patch, {accession: _record(accession) for _, accession in PANEL})
    ledger = reference_source(GENBANK_SOURCE).build_ledger(inventory, timeout=1)
    patch.setattr(viral_panel_module, "_primary_ledger",
                  lambda collection, rows: ledger.loc[rows["abbreviation"]])
    root = tmp_path_factory.mktemp("real-records")
    download_viral_references(species="human", cache_root=str(root))
    patch.undo()
    return str(root)


def _offered(collection, accession):
    reference = viral_reference(accession, cache_root=collection)
    return reference, list(reference.orfs)


def test_every_real_record_becomes_an_installed_reference(collection):
    found = {item.id for item in viral_references(cache_root=collection)}
    assert {accession for _, accession in PANEL} <= found


def test_the_group_is_the_genus_not_the_last_lineage_name(collection):
    groups = {accession: viral_reference(accession, cache_root=collection).group
              for _, accession in PANEL}
    assert groups == {"X04615.1": "Orthohepadnavirus", "AF033817.1": "Deltaretrovirus",
                      "X70828.1": "Mupapillomavirus", "AJ507799.2": "Lymphocryptovirus",
                      "KF006398.1": "Gammapapillomavirus", "KC108721.1": "Gammapapillomavirus"}


def test_the_isolation_host_is_taken_from_the_record_when_it_names_one(collection):
    hosts = {accession: viral_reference(accession, cache_root=collection).host_taxon_ids
             for _, accession in PANEL}
    assert hosts["KF006398.1"] == hosts["KC108721.1"] == (9606,)
    assert hosts["X04615.1"] == hosts["AF033817.1"] == hosts["X70828.1"] == ()
    named = viral_reference("KF006398.1", cache_root=collection)
    assert "X04615.1" not in {item.id for item in viral_references(
        host="human", cache_root=collection)}
    assert "KF006398.1" in named.host_source


def test_hbv_polymerase_wraps_the_circular_origin_and_is_offered(collection):
    reference, offered = _offered(collection, "X04615.1")
    assert reference.coding_models == 4 and "not offered" not in reference.scope
    wrapped = [orf for orf in offered
               if any(a <= 10 for a, _ in orf.coding_blocks)
               and any(b >= reference.length - 10 for _, b in orf.coding_blocks)]
    assert len(wrapped) == 1
    (orf,) = wrapped
    assert len(orf.coding_sequence) % 3 == 0
    assert sum(b - a + 1 for a, b in orf.coding_blocks) == len(orf.coding_sequence)


def test_htlv1_slippage_cds_that_read_a_base_twice_are_left_out(collection):
    reference, offered = _offered(collection, "AF033817.1")
    assert reference.coding_models == 4
    assert "2 of 6 source CDS not offered as ORFs" in reference.scope
    for orf in offered:  # What is offered accounts for every base once.
        assert sum(b - a + 1 for a, b in orf.coding_blocks) == len(orf.coding_sequence)


def test_hpv63_cds_whose_start_does_not_match_its_source_protein_is_left_out(collection):
    reference, _ = _offered(collection, "X70828.1")
    assert reference.coding_models == 6
    assert "1 of 7 source CDS not offered as ORFs" in reference.scope


@pytest.mark.parametrize("accession", ["KF006398.1", "KC108721.1"])
def test_a_cds_with_an_ambiguous_base_is_left_out_not_guessed(collection, accession):
    reference, offered = _offered(collection, accession)
    assert reference.coding_models == 6
    assert "1 of 7 source CDS not offered as ORFs" in reference.scope
    assert all(set(orf.coding_sequence) <= set("ACGT") for orf in offered)


def test_ebv_lmp2_spliced_back_across_the_terminal_repeat_is_left_out(collection):
    """The trimmed record keeps LMP2's two CDSs: one inside the repeat, one spliced across it."""
    reference, offered = _offered(collection, "AJ507799.2")
    assert reference.coding_models == len(offered) == 1
    assert "1 of 2 source CDS not offered as ORFs" in reference.scope
    # The CDS that begins at 166103 and continues from base 58 is the one left out.
    assert all(a < 100000 for orf in offered for a, _ in orf.coding_blocks)


def test_hpv16_e6_and_e6_star_stay_distinct_through_their_sanitized_names():
    """The PaVE HPV16 record names E6 and E6*, which sanitize to the same token."""
    from orfdata.record_sources import PAVE_SOURCE

    pave = (DATA / "HPV16REF.json.gz").read_bytes()
    patch = pytest.MonkeyPatch()
    try:
        serve_records(patch, alternate={"HPV16REF": gzip.decompress(pave)})
        genome = reference_source(PAVE_SOURCE).record("HPV16REF", timeout=1)
    finally:
        patch.undo()
    features = genome.coding_features
    proteins = [feature["source_protein_id"] for feature in features]
    assert any(name.endswith("E6*") for name in proteins) and "HPV16-E6" in proteins
    named = name_transcripts("HPV16", features, unique_protein_ids=unique_source_proteins(features))
    ids = [row["protein_id"] for row in named]
    assert len(set(ids)) == len(ids)
