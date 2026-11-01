"""Installed collections answer viral_references, viral_orfs and viral_panel_coverage."""

import pytest

from orfdata import (
    viral_orfs,
    viral_panel_coverage,
    viral_reference,
    viral_reference_context,
    viral_reference_sequence,
    viral_references,
)
from orfdata.viral_panel import download_viral_references, viral_genome


def _install(tmp_path):
    return download_viral_references(species="human", cache_root=str(tmp_path))


def test_nothing_installed_leaves_the_bundled_catalogue_alone(tmp_path):
    assert [item.id for item in viral_references(cache_root=str(tmp_path))] == \
        [item.id for item in viral_references(installed=False)]
    assert viral_panel_coverage(cache_root=str(tmp_path)).equals(
        viral_panel_coverage(installed=False))


def test_an_installed_collection_lists_its_references(viral_panel, tmp_path):
    _install(tmp_path)
    found = {item.id: item for item in viral_references(cache_root=str(tmp_path))}
    assert {"VTEST001.1", "VTEST002.1", "VTEST003.1", "VTEST004.1"} <= set(found)
    first = found["VTEST001.1"]
    assert first.bundle is not None and first.virus_type == "TV1"
    assert first.group == "" and first.host_taxon_ids == ()
    assert first.coding_models == 1
    assert [item.protein_id for item in first.products] == ["VTEST001.1_cds1.p"]
    assert viral_reference("VTEST001.1", cache_root=str(tmp_path)) == first
    assert viral_references(virus_type="tv1", cache_root=str(tmp_path)) == (first,)
    # Host and bundled-only queries are unchanged by what is installed.
    assert viral_references(host="human", cache_root=str(tmp_path)) == \
        viral_references(host="human", installed=False)
    assert all(item.bundle is None for item in viral_references(installed=False))
    with pytest.raises(ValueError, match="No unique"):
        viral_reference("VTEST001.1", installed=False)


def test_installed_dna_and_orfs_load_from_the_bundle(viral_panel, tmp_path):
    _install(tmp_path)
    root = str(tmp_path)
    first = viral_reference("VTEST001.1", cache_root=root)
    assert viral_reference_sequence("VTEST001.1", cache_root=root) == first.dna
    assert len(first.dna) == first.length

    (orf,) = first.orfs
    assert orf.protein_sequence == "MAAAAAAAA"
    assert orf.coding_sequence == "ATG" + "GCT" * 8 + "TAA"
    # The same ORFs load by exact sequence, with the caller's own contig name.
    (named,) = viral_orfs(first.dna, contig="chrV", cache_root=root)
    assert named.protein_sequence == orf.protein_sequence
    assert named.transcript.contig == "chrV"
    # A CDS that is partial, or whose source protein is not what its bases
    # translate to, stays in the transcriptome layer and is counted as not offered.
    second = viral_reference("VTEST002.1", cache_root=root)
    assert second.coding_models == 0 and second.orfs == ()
    assert "2 of 2 source CDS not offered" in second.scope
    assert viral_orfs(reference_id="VTEST001.1", cache_root=root)[0].orf_id == orf.orf_id
    with pytest.raises(ValueError, match="No unique released viral"):
        viral_orfs(first.dna, installed=False)


def test_minus_strand_cds_is_offered_and_an_origin_wrap_is_counted_as_not(
        viral_panel, tmp_path, monkeypatch):
    from conftest import _genbank, serve_records, viral_documents

    from orfdata import viral_panel as module
    from orfdata.record_sources import GENBANK_SOURCE, reference_source

    coding = "ATG" + "GCT" * 8 + "TAA"
    complement = coding[::-1].translate(str.maketrans("ACGT", "TGCA"))
    records = viral_documents()
    records["VTEST001.1"] = _genbank("VTEST001.1", complement + "ACGT" * 10, [
        ("complement(1..30)", {"gene": "E1", "product": "protein E1", "codon_start": 1,
                               "translation": "MAAAAAAAA"})], name="VTEST001", taxon=10001)
    # The reading frame runs off the end of the record and wraps to its start,
    # which SourceRNA cannot yet hold (#66); the transcriptome layer still has it.
    records["VTEST002.1"] = _genbank("VTEST002.1", coding[15:] + "ACGT" * 7 + "AC" + coding[:15], [
        ("join(46..60,1..15)", {"gene": "P", "product": "protein P", "codon_start": 1,
                                "translation": "MAAAAAAAA"})], name="VTEST002", taxon=10002)
    serve_records(monkeypatch, records)
    ledger = reference_source(GENBANK_SOURCE).build_ledger(viral_panel, timeout=1)
    monkeypatch.setattr(module, "_primary_ledger",
                        lambda collection, rows: ledger.loc[rows["abbreviation"]])
    _install(tmp_path)
    root = str(tmp_path)
    (reverse,) = viral_reference("VTEST001.1", cache_root=root).orfs
    assert reverse.transcript.strand == "-" and reverse.protein_sequence == "MAAAAAAAA"
    assert reverse.coding_sequence == coding
    wrapped = viral_reference("VTEST002.1", cache_root=root)
    assert wrapped.coding_models == 0 and "1 of 1 source CDS not offered" in wrapped.scope


def test_the_coverage_ledger_reports_installed_references(viral_panel, tmp_path, monkeypatch):
    from orfdata import viral_installed
    from orfdata.viral_registry import GDC_VIRUSES

    _install(tmp_path)
    root = str(tmp_path)
    held = viral_installed.installed_coverage(GDC_VIRUSES, root)
    assert held == {"TV1": ("VTEST001.1", 1), "TV2": ("VTEST002.1", 0),
                    # TV3 declares NCBI table 11, which ORFdata cannot yet translate.
                    "TV3": ("VTEST003.1", 0), "TV4": ("VTEST004.1", 0)}

    # The fixture's abbreviations are not GDC rows, so answer for two real ones.
    monkeypatch.setattr(viral_installed, "installed_coverage", lambda *a, **k: {
        "HPV18": ("VTEST001.1", 1), "HPV31": ("VTEST002.1", 0), "HPV16": ("VTEST003.1", 1)})
    frame = viral_panel_coverage(cache_root=root).set_index("abbreviation")
    assert frame.loc["HPV18", ["status", "bundled_reference_id", "coding_models", "available"]].tolist() == \
        ["source_cds_installed", "VTEST001.1", 1, True]
    assert frame.loc["HPV31", ["status", "bundled_reference_id", "coding_models", "available"]].tolist() == \
        ["installed_without_cds", "", 0, False]
    # A row the bundled catalogue already covers is left as it was.
    assert frame.loc["HPV16", "status"] == "source_cds_bundled"
    assert frame.loc["HPV16", "bundled_reference_id"] != "VTEST003.1"
    assert frame["status"].value_counts().to_dict() == {
        "not_bundled": 197, "source_cds_bundled": 1, "source_cds_installed": 1,
        "installed_without_cds": 1}
    assert viral_panel_coverage(cache_root=root, installed=False)["status"].value_counts().to_dict() == {
        "not_bundled": 199, "source_cds_bundled": 1}


def test_the_viral_genome_is_declared_for_downstream_consumers(viral_panel, tmp_path):
    pytest.importorskip("pyensembl")
    from orfdata import ReferenceContext

    bundle = _install(tmp_path)
    genome = viral_genome(bundle, cache_root=str(tmp_path))
    context = viral_reference_context(bundle)
    assert (context.taxon_id, context.species) == (9606, "homo_sapiens")
    assert context.genome_build == genome.reference_name
    assert context.annotation_name == genome.annotation_name
    assert context.annotation_version == genome.annotation_version
    assert ReferenceContext.from_genome(genome, context=context) == context
