"""Qualify real, versioned viral CDS paths independently of sample assembly."""

import gzip
import hashlib
import io
import json
from copy import deepcopy
from importlib.resources import files
from pathlib import Path

import pytest
from Bio import SeqIO

import orfdata.viral as viral_module
from orfdata import SourceRNA, viral_orfs, viral_reference_sequence, viral_references
from tools.build_hpv16_cds import build

FIXTURE = Path(__file__).parent / "data/viral/NC_001526.4.gb.gz"
EXPECTED = {
    "YP_009268708.1": 92, "NP_041327.2": 649, "NP_041328.1": 365,
    "NP_041330.2": 83, "NP_041331.2": 473, "NP_041332.2": 505,
    "NP_041325.1": 158, "YP_009268707.1": 50, "NP_041326.1": 98,
}


@pytest.fixture
def source():
    raw = gzip.decompress(FIXTURE.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == "282be442f226cc809fefbda616ac136b8d20dfde23fad09558f3f653a3de96a4"
    return raw, SeqIO.read(io.StringIO(raw.decode()), "genbank")


def test_catalogue_rebuilds_from_the_pinned_public_annotation(source):
    raw, _ = source
    packaged = json.loads(gzip.decompress(files("orfdata").joinpath(
        "data/viral/hpv16-refseq-CDS.json.gz").read_bytes()))
    assert build(raw) == packaged
    assert len(packaged["cds"]) == 9


def test_native_source_paths_match_all_nine_source_proteins(source):
    _, record = source
    orfs = viral_orfs(str(record.seq))
    definitions = {f.qualifiers["protein_id"][0]: f for f in record.features if f.type == "CDS"}
    assert {orf.orf_id: len(orf.protein_sequence) for orf in orfs} == EXPECTED
    for orf in orfs:
        feature = definitions[orf.orf_id]
        assert isinstance(orf.transcript, SourceRNA)
        assert orf.transcript.biotype == "source_CDS_fragment"
        assert orf.transcript.genome is None
        assert orf.protein_sequence == feature.qualifiers["translation"][0]
        assert orf.coding_sequence == str(feature.extract(record.seq))
        assert orf.start_complete and orf.stop_complete
        assert orf.context.annotation_version == "NC_001526.4"
        assert orf.context.taxon_id == 333760
        assert not orf.evidence.has_evidence("ms")
        assert not orf.evidence.has_evidence("ribosome")


def test_gdc_rotation_is_exact_and_retains_spliced_coding_positions(source):
    _, record = source
    native = str(record.seq)
    gdc = native[-864:] + native[:-864]
    assert hashlib.md5(gdc.encode()).hexdigest() == "f2463721818c29998ef5ab917e6df74a"
    orfs = {orf.orf_id: orf for orf in viral_orfs(gdc, contig="HPV16")}
    assert orfs["YP_009268708.1"].coding_blocks == ((865, 880), (3358, 3620))
    assert orfs["YP_009268707.1"].coding_blocks == ((83, 226), (409, 417))
    assert orfs["NP_041326.1"].coding_blocks == ((562, 858),)
    for orf in orfs.values():
        assert orf.transcript.contig == "HPV16"
        assert orf.context.genome_build == "GRCh38.d1.vd1:HPV16"
        assert orf.protein_sequence == next(o.protein_sequence for o in viral_orfs(native)
                                            if o.orf_id == orf.orf_id)


def test_mutated_same_name_reference_and_partial_dna_are_not_qualified(source):
    _, record = source
    native = str(record.seq)
    changed = ("A" if native[0] != "A" else "C") + native[1:]
    for sequence in (changed, native[1:], native[::-1]):
        with pytest.raises(ValueError, match="No unique released"):
            viral_orfs(sequence, contig="HPV16")


def test_discovery_is_offline_and_returned_metadata_cannot_change_loading(source):
    _, record = source
    references = viral_references()
    assert {row["id"] for row in references} == {"NC_001526.4", "GRCh38.d1.vd1:HPV16"}
    assert all(row["coding_models"] == 9 for row in references)
    assert all(row["catalogue_sha256"] and row["source"]["sha256"] for row in references)
    assert all("reference_sequence" not in row for row in references)
    references[0]["source"]["accession"] = "changed"
    references[0]["host_taxon_ids"].clear()
    assert viral_references()[0]["source"]["accession"] == "NC_001526.4"
    assert viral_references()[0]["host_taxon_ids"] == [9606]
    assert viral_orfs(str(record.seq))[0].context.annotation_version == "NC_001526.4"


def test_host_group_and_type_select_exact_references_with_distinct_taxa():
    references = viral_references(host="human", group="HPV", virus_type="HPV16")
    assert {row["id"] for row in references} == {"NC_001526.4", "GRCh38.d1.vd1:HPV16"}
    assert viral_references(host=9606, group="hpv", virus_type="hpv16") == references
    assert viral_references(host="mouse") == viral_references(host=10090) == ()
    assert viral_references(host="human", virus_type="HPV18") == ()
    for row in references:
        assert row["taxon_id"] == 333760 and row["host_taxon_ids"] == [9606]
        assert "PaVE HPV16REF" in row["host_source"]
        assert row["virus"]["ictv_species"] == "Alphapapillomavirus 9"
        assert row["virus"]["group"] == "HPV" and row["virus"]["type"] == "HPV16"
        assert row["virus"]["taxonomy_source"].startswith("https://ictv.global/")
    for invalid in ("not_a_host", 0, True):
        with pytest.raises(ValueError, match="host"):
            viral_references(host=invalid)


def test_explicit_multiple_host_associations_can_select_one_reference_from_each_host(monkeypatch):
    # Synthetic metadata checks selection only; it makes no HPV16 mouse-host claim.
    catalog = deepcopy(viral_module._catalog())
    catalog["references"][0]["host_taxon_ids"].append(10090)
    catalog["references"][0]["host_source"] = "synthetic multiple-host fixture"
    monkeypatch.setattr(viral_module, "_catalog", lambda: catalog)
    assert {row["id"] for row in viral_references(host="human")} == {
        "NC_001526.4", "GRCh38.d1.vd1:HPV16"}
    mouse, = viral_references(host="mouse", virus_type="HPV16")
    assert mouse["id"] == "NC_001526.4" and mouse["taxon_id"] == 333760


def test_bundled_reference_sequences_and_named_orfs_match_both_source_placements(source):
    _, record = source
    native = str(record.seq)
    gdc = native[-864:] + native[:-864]
    for reference_id, expected in (("NC_001526.4", native),
                                   ("GRCh38.d1.vd1:HPV16", gdc)):
        sequence = viral_reference_sequence(reference_id)
        assert sequence == expected
        row, = (r for r in viral_references() if r["id"] == reference_id)
        assert hashlib.sha256(sequence.encode()).hexdigest() == row["sequence_sha256"]
        by_name = viral_orfs(reference_id=reference_id, contig="selected")
        by_dna = viral_orfs(sequence, contig="selected")
        assert [(o.orf_id, o.coding_blocks, o.protein_sequence) for o in by_name] == [
            (o.orf_id, o.coding_blocks, o.protein_sequence) for o in by_dna]
        assert {o.transcript.contig for o in by_name} == {"selected"}
    with pytest.raises(ValueError, match="reference ID"):
        viral_reference_sequence("unknown")
    with pytest.raises(ValueError, match="exactly one"):
        viral_orfs()
    with pytest.raises(ValueError, match="exactly one"):
        viral_orfs(native, reference_id="NC_001526.4")


def test_corrupt_packaged_catalogue_fails_before_returning_models(monkeypatch):
    import orfdata.viral as module

    class CorruptResource:
        def joinpath(self, path):
            return self

        def read_bytes(self):
            return b"corrupt catalogue"

    module._catalog.cache_clear()
    try:
        monkeypatch.setattr(module, "files", lambda name: CorruptResource())
        with pytest.raises(ValueError, match="checksum mismatch"):
            viral_references()
    finally:
        module._catalog.cache_clear()
