"""Catalogue ORFs export as GTF, so PyEnsembl-based tools can annotate variants in them."""

import pytest

from orfdata import load_catalog
from orfdata.cli import run
from orfdata.gtf import orf_gtf_lines

CATALOG = load_catalog()


def _rows(lines):
    return [line.rstrip("\n").split("\t") for line in lines]


def _by_feature(rows, feature):
    return [row for row in rows if row[2] == feature]


def _orf(orf_id):
    record, = (r for r in CATALOG.query() if r.orf_id == orf_id)
    return record


def test_each_orf_becomes_one_coding_transcript_on_its_host_gene():
    record = _orf("c10norep1")
    rows = _rows(orf_gtf_lines([record]))
    assert [row[2] for row in rows] == ["transcript", "exon", "CDS", "start_codon", "stop_codon"]
    assert {row[0] for row in rows} == {record.contig}
    assert {row[6] for row in rows} == {"-"}
    assert all('transcript_id "c10norep1";' in row[8] for row in rows)
    assert 'gene_id "ENSG00000107929"' in rows[0][8]
    assert 'host_transcript_id "ENST00000316157"' in rows[0][8]
    assert 'orf_biotype "intORF"' in rows[0][8]
    assert 'protein_id "c10norep1.p"' in _by_feature(rows, "CDS")[0][8]
    # No gene rows: the file concatenates onto a reference annotation that has them.
    assert not _by_feature(rows, "gene")


@pytest.mark.parametrize("orf_id", ["c10norep119", "c10norep169", "c10norep1"])
def test_exons_cds_and_frames_describe_the_whole_reading_frame(orf_id):
    record = _orf(orf_id)
    rows = _rows(orf_gtf_lines([record]))
    exons = _by_feature(rows, "exon")
    assert len(exons) == len(record.coding_blocks)
    # Exons are written in transcript order, which runs backwards on the minus strand.
    expected = sorted(record.coding_blocks, reverse=record.strand == "-")
    assert [(int(row[3]), int(row[4])) for row in exons] == [tuple(b) for b in expected]

    coding = _by_feature(rows, "CDS")
    assert sum(int(row[4]) - int(row[3]) + 1 for row in coding) == 3 * len(record.protein_sequence)
    offset = 0
    for row in coding:  # Frames continue the reading frame across every block.
        assert row[7] == str(-offset % 3)
        offset += int(row[4]) - int(row[3]) + 1

    start, = [row for row in rows if row[2] == "start_codon"][:1]
    stop, = [row for row in rows if row[2] == "stop_codon"][:1]
    first, last = expected[0], expected[-1]
    if record.strand == "+":
        assert int(start[3]) == first[0] and int(stop[4]) == last[1]
    else:
        assert int(start[4]) == first[1] and int(stop[3]) == last[0]


def test_cli_writes_gtf_for_a_query(capsys):
    run(["--gene", "LARP4B", "--format", "gtf"])
    rows = _rows(capsys.readouterr().out.splitlines())
    assert rows and {row[1] for row in rows} == {"orfdata"}
    assert {row[2] for row in rows} <= {"transcript", "exon", "CDS", "start_codon", "stop_codon"}


def test_pyensembl_reads_the_exported_orfs(tmp_path):
    pytest.importorskip("pyensembl")
    from pyensembl import Genome

    records = [_orf("c10norep119"), _orf("c10norep169")]
    path = tmp_path / "ncorfs.gtf"
    path.write_text("".join(orf_gtf_lines(records)))
    genome = Genome(reference_name="GRCh38", annotation_name="orfdata-ncorfs",
                    annotation_version="test", gtf_path_or_url=str(path),
                    cache_directory_path=str(tmp_path / "index"))
    genome.index()

    transcript = genome.transcript_by_id("c10norep119")
    assert transcript.contig == "10" and transcript.strand == "-"
    assert transcript.contains_start_codon and transcript.contains_stop_codon
    assert [(exon.start, exon.end) for exon in transcript.exons] == [
        tuple(block) for block in sorted(records[0].coding_blocks, reverse=True)]
    assert genome.gene_names_at_locus("10", records[0].coding_blocks[0][0]) == [
        records[0].gene_name]
    assert genome.transcript_by_id("c10norep169").strand == "+"


def test_proteome_writes_a_matching_annotation_for_its_nc_orfs(tmp_path):
    pytest.importorskip("pyensembl")
    from pyensembl import Genome
    from test_proteome import make_catalog, write_genome

    from orfdata import Proteome

    genome = write_genome(tmp_path / "genome")
    proteome = Proteome.build(genome, nc_orf_sources=make_catalog())
    written = proteome.write_annotation(tmp_path / "ncorfs")
    assert written["transcript_count"] == written["protein_count"] == 3
    assert written["omitted_entries"] == sum(e.kind != "extra_orf" for e in proteome.proteins)

    annotated = Genome(reference_name="GRCm39", annotation_name="orfdata-ncorfs",
                       annotation_version="1", gtf_path_or_url=written["annotation_gtf"],
                       transcript_fasta_paths_or_urls=[written["transcript_fasta"]],
                       protein_fasta_paths_or_urls=[written["protein_fasta"]],
                       cache_directory_path=str(tmp_path / "index"))
    annotated.index()
    transcript = annotated.transcript_by_id("ms-orf")
    # The spliced ncORF reads back as a complete coding transcript.
    assert [(exon.start, exon.end) for exon in transcript.exons] == [(103, 104), (205, 217)]
    assert transcript.complete and transcript.biotype == "protein_coding"
    assert transcript.sequence == "ATGAAACCCGGGTAA"
    assert transcript.protein_sequence == "MKPG"
    genome.close()


def test_varcode_gives_a_variant_in_an_nc_orf_a_coding_consequence(tmp_path):
    pytest.importorskip("varcode")
    from pyensembl import Genome
    from test_proteome import make_catalog, write_genome
    from varcode import Variant

    from orfdata import Proteome

    genome = write_genome(tmp_path / "genome")
    written = Proteome.build(genome, nc_orf_sources=make_catalog()).write_annotation(
        tmp_path / "ncorfs")
    annotated = Genome(reference_name="GRCm39", annotation_name="orfdata-ncorfs",
                       annotation_version="1", gtf_path_or_url=written["annotation_gtf"],
                       transcript_fasta_paths_or_urls=[written["transcript_fasta"]],
                       protein_fasta_paths_or_urls=[written["protein_fasta"]],
                       cache_directory_path=str(tmp_path / "index"))
    annotated.index()
    effects = Variant(contig="1", start=206, ref="A", alt="T", genome=annotated).effects()
    assert effects.top_priority_effect().short_description == "p.K2*"
    # The fixture's ORFs share a coding path, so every ncORF carrying it is annotated.
    assert {effect.transcript_id: effect.short_description for effect in effects} == {
        "ms-orf": "p.K2*", "absent-orf": "p.K2*"}
    genome.close()
