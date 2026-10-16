"""Qualify real, versioned viral CDS paths independently of sample assembly."""

import gzip
import hashlib
import io
import json
from importlib.resources import files
from pathlib import Path

import pytest
from Bio import SeqIO

from orfdata import SourceRNA, viral_orfs, viral_references
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
    references[0]["source"]["accession"] = "changed"
    assert viral_references()[0]["source"]["accession"] == "NC_001526.4"
    assert viral_orfs(str(record.seq))[0].context.annotation_version == "NC_001526.4"


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
