"""Carrier loci survive export; synthetic mappings are not biological observations."""

import json
from pathlib import Path

import pytest
from test_orf_cases import genome as genome

from orfdata import ORFModel, ORFRecord, ReferenceContext, SourceRNA, TranscriptORF


class Carrier:
    """A second RNA with the same pinned sequence/path and a contig alias."""

    def __init__(self, transcript):
        self.original = transcript
        self.id = "carrier-" + transcript.id
        self.versioned_id = self.id + "." + str(transcript.version)
        self.contig = "chr" + transcript.contig
        self.gene_id = "carrier-gene"
        self.gene_name = "CARRIER"

    def __getattr__(self, name):
        return getattr(self.original, name)


@pytest.fixture(params=("ENST00000644866", "ENST00000395826"), ids=("plus", "minus"))
def placement(request, genome):
    transcript = genome.transcript_by_id(request.param)
    annotated = TranscriptORF.from_annotation(transcript)
    record = ORFRecord.from_model(ORFModel(
        "synthetic-record", transcript.id, transcript.gene_id, transcript.gene_name,
        transcript.contig, transcript.strand, annotated.coding_blocks,
        annotated.protein_sequence, "synthetic", "GRCh38", 9606, "homo_sapiens",
        source="synthetic"))
    return transcript, annotated, record


@pytest.mark.parametrize("route", ("at", "annotated", "mapped_annotation", "record", "mapped_record"))
def test_every_constructor_exports_the_actual_carrier_locus(placement, route):
    transcript, annotated, record = placement
    carrier = Carrier(transcript) if route.startswith("mapped_") else transcript
    if route == "at":
        orf = TranscriptORF.at(carrier, annotated.start, orf_id="synthetic-at",
                               source="synthetic", protein=annotated.protein_sequence)
    elif route in ("annotated", "mapped_annotation"):
        orf = TranscriptORF.from_annotation(transcript, on=carrier)
    else:
        orf = TranscriptORF(record, carrier, any_transcript=route == "mapped_record")
    snapshot = json.loads(json.dumps(orf.to_dict()))
    resolved = snapshot["resolved_transcript"]
    assert resolved["id"] == carrier.id
    assert resolved["contig"] == carrier.contig
    assert resolved["strand"] == carrier.strand
    if route.startswith("mapped_"):
        assert not snapshot["on_source_transcript"]
        assert snapshot["source_transcript_id"] != resolved["id"]
    if route == "mapped_record":
        assert snapshot["model"]["contig"] != resolved["contig"]
    if route in ("at", "annotated", "mapped_annotation"):
        assert snapshot["model"] is None and snapshot["record"] is None
    restored = TranscriptORF.from_dict(snapshot, carrier)
    assert json.loads(json.dumps(restored.to_dict())) == snapshot
    # Earlier v1 exports omitted both fields for every construction path.
    legacy = json.loads(json.dumps(snapshot))
    del legacy["resolved_transcript"]["strand"]
    del legacy["resolved_transcript"]["contig"]
    assert TranscriptORF.from_dict(legacy, carrier).to_dict() == orf.to_dict()
    assert "strand" not in legacy["resolved_transcript"]  # Replay does not rewrite its input.


def source_rna(**locus):
    return SourceRNA("snapshot-fixture", "ATGAAATAA",
                     ReferenceContext(9606, "human", "GRCh38", "synthetic", "1"),
                     "synthetic", **locus)


@pytest.mark.parametrize("strand", (None, "+", "-"))
def test_source_rna_locus_is_explicit_even_when_unknown(strand):
    locus = {} if strand is None else dict(contig="2", strand=strand, exons=((100, 108),))
    rna = source_rna(**locus)
    orf = TranscriptORF.at(rna, 0, orf_id="source-rna", source="synthetic")
    snapshot = json.loads(json.dumps(orf.to_dict()))
    assert snapshot["resolved_transcript"]["strand"] == strand
    assert snapshot["resolved_transcript"]["contig"] == (None if strand is None else "2")
    assert TranscriptORF.from_dict(snapshot, rna).to_dict() == orf.to_dict()


def test_published_042_snapshot_remains_readable():
    # Frozen using the published 0.42.0 wheel, before carrier locus fields existed.
    path = Path(__file__).parent / "data/references/transcript_orf_v1_0.42.0.json"
    snapshot = json.loads(path.read_text())
    restored = TranscriptORF.from_dict(snapshot, source_rna())
    assert restored.protein_sequence == "MK" and restored.record is None
    assert restored.to_dict()["resolved_transcript"]["strand"] is None
    assert restored.to_dict()["resolved_transcript"]["contig"] is None
