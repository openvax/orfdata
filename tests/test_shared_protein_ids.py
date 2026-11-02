"""A source protein accession shared by two CDS features must not break reading the bundle."""

from conftest import serve_records, viral_documents

from orfdata import download_viral_references
from orfdata.fasta import read_fasta


def test_a_shared_source_protein_id_is_written_once_per_record_and_reads_back(
        viral_panel, tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents())
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    raw = [line[1:].split()[0] for line in open(bundle.protein_fasta) if line.startswith(">")]
    assert len(raw) == len(set(raw)), f"writer emitted repeated protein ids: {raw}"
    assert not any(name.startswith("AAA00001") for name in raw)  # Shared id was replaced.
    assert len(read_fasta(bundle.protein_fasta)) == len(raw)


def test_accessions_that_sanitize_to_one_token_are_not_both_kept():
    """HPV16's E6 and E6* have accessions HPV16-E6 and HPV16-E6*; the files strip the asterisk."""
    from orfdata.viral_annotation import name_transcripts, unique_source_proteins

    def feature(ordinal, protein):
        return {"ordinal": ordinal, "source_protein_id": protein, "exons": ((1, 9, "+"),),
                "gene_label": "E6"}

    features = [feature(1, "HPV16-E6*"), feature(2, "HPV16-E6"), feature(3, "OTHER.1")]
    named = name_transcripts("G", features, unique_protein_ids=unique_source_proteins(features))
    ids = [row["protein_id"] for row in named]
    assert ids == ["G_cds1.p", "G_cds2.p", "OTHER.1"]
    assert len(set(ids)) == len(ids)
