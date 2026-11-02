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
