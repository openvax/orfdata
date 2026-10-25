"""The first FASTA word maps source records back to an offline proteome."""

import gzip
import io
import json
from urllib.parse import unquote

from orfdata import ProducerInfo, ProteinSequence, Proteome, ReferenceContext


def _tokens(text):
    return [line[1:].split(maxsplit=1)[0] for line in text.splitlines()
            if line.startswith(">")]


def _fields(token):
    return {name: unquote(value) for name, value in
            (part.split("=", 1) for part in token.split("|"))}


def test_readable_tokens_resolve_after_offline_load_and_preserve_sources(tmp_path):
    reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")
    cds = ProteinSequence("ENSP1.1", "MPEPTIDE", "annotated", 9606, "ensembl:115")
    nc = ProteinSequence("ORF1", "MPEP", "extra_orf", 9606, "catalogue:fixture",
                       orf_type="uORF")
    varcode = ProteinSequence.from_producer("MVAR", taxon_id=9606, producer=ProducerInfo(
        "reference_cds", "dna_prediction", "varcode", "V 1|x", reference=reference,
        sample_id="tumor A", gene_names=("GENE|A",)))
    products = tuple(ProteinSequence.from_producer(
        "MPEP", taxon_id=9606, producer=ProducerInfo(
            "nc_orf", "rna_assembly", "isovar", "ASM1", reference=reference,
            sample_id=sample, gene_names=("GENE 1", "GENE|2"), orf_kind="uORF",
            start_complete=False, stop_complete=None)) for sample in ("A", "B"))
    missing = ProteinSequence.from_producer("M", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "no gene", reference=reference))
    proteome = Proteome(reference, (cds, nc)).with_proteins((varcode, *products, missing))
    path = tmp_path / "joint.proteome.json.gz"
    proteome.save(path)
    loaded = Proteome.load(path)
    output = io.StringIO()
    loaded.write_fasta(output)
    tokens = _tokens(output.getvalue())
    index = loaded.fasta_index()
    assert len(tokens) == len(index) == len(loaded) == 6
    assert set(tokens) == set(index)
    assert all(len(group) == 1 for group in index.values())
    assert all(index[token][0].accession in loaded.by_accession() for token in tokens)
    assert all(_fields(token)["aa"] == str(len(index[token][0].sequence)) for token in tokens)
    assert _fields(cds.fasta_token())["origin"] == "reference_cds"
    assert _fields(nc.fasta_token())["orf"] == "uORF"
    assert _fields(varcode.fasta_token())["tool"] == "varcode"
    assert _fields(varcode.fasta_token())["id"] == "V 1|x"
    assert "%20" in varcode.fasta_token() and "%7C" in varcode.fasta_token()
    assert _fields(products[0].fasta_token())["gene"] == "GENE 1,GENE|2"
    assert _fields(products[0].fasta_token())["gene_count"] == "2"
    assert _fields(products[0].fasta_token())["sample"] == "A"
    assert _fields(products[1].fasta_token())["sample"] == "B"
    assert _fields(missing.fasta_token())["gene"] == "unknown"
    assert products[0].sequence == "MPEP" and products[0].producer.start_complete is False
    assert loaded.fasta_index() == proteome.fasta_index()


def test_exact_dedup_is_opt_in_and_indexes_every_source(tmp_path):
    reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")
    first = ProteinSequence.from_producer("MPEP", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P1", reference=reference, sample_id="A"))
    second = ProteinSequence.from_producer("MPEP", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P1", reference=reference, sample_id="B"))
    partial = ProteinSequence.from_producer("PEP", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P2", reference=reference,
        start_complete=False))
    proteome = Proteome.from_proteins((first, second, partial), reference=reference)
    normal = io.StringIO()
    proteome.write_fasta(normal)
    assert len(_tokens(normal.getvalue())) == 3
    grouped = io.StringIO()
    proteome.write_fasta(grouped, deduplicate=True)
    tokens = _tokens(grouped.getvalue())
    assert len(tokens) == 2
    index = proteome.fasta_index(deduplicate=True)
    assert {tuple(entry.accession for entry in group) for group in index.values()} == {
        (first.accession, second.accession), (partial.accession,)}
    assert any("|records=2|group=" in token for token in tokens)
    assert "\nPEP\n" in grouped.getvalue()
    path = tmp_path / "joint.json.gz"
    proteome.save(path)
    assert Proteome.load(path).fasta_index(deduplicate=True) == index


def test_v3_external_entry_snapshot_keeps_its_original_identity(tmp_path):
    reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")
    entry = ProteinSequence.from_producer("MAA", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P1", reference=reference))
    collection = Proteome.from_proteins([entry], reference=reference)
    path = tmp_path / "v3.json.gz"
    collection.save(path)
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    snapshot["schema"] = "orfdata.proteome.v3"
    assert "orf_kind" not in snapshot["content"]["entries"][0]["producer"]
    path.write_bytes(gzip.compress(json.dumps(snapshot).encode(), mtime=0))
    loaded = Proteome.load(path)
    assert loaded.identity == collection.identity
