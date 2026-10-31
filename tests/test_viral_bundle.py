"""Viral reference files keep genome, transcriptome and proteome layers separate."""

import gzip
import json

import pytest

from orfdata.gtf import transcript_blocks
from orfdata.viral_annotation import split_genbank
from orfdata.viral_panel import (
    VIRAL_LAYERS,
    ViralReferenceFiles,
    download_viral_references,
    open_viral_references,
    select_viral_references,
    viral_reference_summary,
)

CONTEXT = "cancer-transcriptome-screening"


def _install(tmp_path):
    return download_viral_references(species="human", cache_root=str(tmp_path))


def test_bundle_separates_the_three_layers(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    assert isinstance(bundle, ViralReferenceFiles)
    assert (bundle.genome_count, bundle.transcript_count, bundle.protein_count) == (4, 4, 4)
    assert [layer.name for layer in bundle.layers] == ["genome", "transcriptome", "proteome"]
    assert bundle.layer("genome").molecule == "dna"
    assert bundle.layer("transcriptome").annotation == bundle.annotation_gtf
    assert bundle.layer("proteome").records == 4  # One record annotates no CDS at all.
    with pytest.raises(ValueError, match="Unknown viral layer"):
        bundle.layer("exome")

    genomes = bundle.genome_fasta.read_text().splitlines()
    assert genomes[0].startswith(">VTEST001.1 abbreviation=TV1 taxon=10001 molecule=genome")
    assert bundle.to_dict()["layers"][1]["fasta"] == str(bundle.transcript_fasta)


def test_transcripts_and_proteins_share_identifiers(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    transcripts = {line.split()[0][1:]: line for line in
                   bundle.transcript_fasta.read_text().splitlines() if line.startswith(">")}
    proteins = {line.split()[0][1:] for line in
                bundle.protein_fasta.read_text().splitlines() if line.startswith(">")}
    assert set(transcripts) == {"VTEST001.1_cds1", "VTEST002.1_cds1", "VTEST002.1_cds2",
                                "VTEST003.1_cds1"}
    # A GenBank protein accession names a protein only when it names exactly one.
    assert proteins == {"VTEST001.1_cds1.p", "VTEST002.1_cds1.p", "VTEST002.1_cds2.p",
                        "AAA00003.1"}

    rows = {row["transcript_id"]: row for row in
            [dict(zip(*(line.split("\t") for line in
                        [bundle.transcript_table.read_text().splitlines()[0], body])))
             for body in bundle.transcript_table.read_text().splitlines()[1:]]}
    assert rows["VTEST001.1_cds1"]["annotation_complete"] == "yes"
    assert rows["VTEST002.1_cds1"]["codon_start"] == "2"
    assert rows["VTEST002.1_cds1"]["annotation_complete"] == "no"  # 59 coding bases, not whole codons
    assert rows["VTEST002.1_cds1"]["exon_count"] == "2"
    assert rows["VTEST003.1_cds1"]["translation_table"] == "11"
    assert rows["VTEST001.1_cds1"]["source_protein_id"] == "AAA00001.1"


def test_bundle_keeps_the_source_genbank_text_byte_for_byte(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    stored = sorted(bundle.record_directory.glob("*.gz"))
    # Every document is named for the entry, the source that served it and its identifier.
    assert [path.name for path in stored] == [
        "TV1__genbank__VTEST001.1.gb.gz", "TV2__genbank__VTEST002.1.gb.gz",
        "TV3__genbank__VTEST003.1.gb.gz", "TV4__genbank__VTEST004.1.gb.gz"]
    text = gzip.decompress(stored[0].read_bytes()).decode("utf-8")
    assert text.startswith("LOCUS") and text.rstrip().endswith("//")
    # Each stored record is what EFetch returns for that accession alone, blank
    # separator line and all, so splitting one back out is the identity.
    assert split_genbank(text.encode("utf-8")) == (text,)
    batch = b"".join(gzip.decompress(path.read_bytes()) for path in stored)
    assert split_genbank(batch) == tuple(
        gzip.decompress(path.read_bytes()).decode("utf-8") for path in stored)


def test_annotation_places_transcripts_on_the_viral_genomes(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    rows = [line.split("\t") for line in bundle.annotation_gtf.read_text().splitlines()
            if not line.startswith("#")]
    assert {row[0] for row in rows} == {"VTEST001.1", "VTEST002.1", "VTEST003.1"}
    assert {row[1] for row in rows} == {"orfdata"}

    first = [row for row in rows if 'transcript_id "VTEST001.1_cds1"' in row[8]]
    kinds = [row[2] for row in first]
    assert kinds == ["transcript", "exon", "CDS", "start_codon", "stop_codon"]
    transcript, exon, cds, start, stop = first
    assert (transcript[3], transcript[4], transcript[6]) == ("1", "30", "+")
    assert (cds[3], cds[4], cds[7]) == ("1", "27", "0")  # CDS stops before the stop codon.
    assert (start[3], start[4]) == ("1", "3") and (stop[3], stop[4]) == ("28", "30")
    assert 'transl_table "1";' in transcript[8]

    # The spliced, minus-strand, 5'-partial CDS: exons in transcript order, no codon rows.
    second = [row for row in rows if 'transcript_id "VTEST002.1_cds1"' in row[8]]
    exons = [(row[3], row[4]) for row in second if row[2] == "exon"]
    assert exons == [("61", "90"), ("11", "40")]
    assert all(row[6] == "-" for row in second)
    assert [row[2] for row in second].count("start_codon") == 0
    coding = [(row[3], row[4], row[7]) for row in second if row[2] == "CDS"]
    # codon_start=2 drops the first transcript base, which is the highest coordinate;
    # 29 coding bases then precede the second block, so one of its bases finishes a codon.
    assert coding == [("61", "89", "0"), ("11", "40", "1")]
    assert 'transl_table "11";' in [row[8] for row in rows
                                    if 'transcript_id "VTEST003.1_cds1"' in row[8]][0]

    genes = [row for row in rows if row[2] == "gene"]
    assert len(genes) == 4  # E6 on two strands of VTEST002 stays two genes.


def test_empty_trailing_table_values_stay_quoted(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    for table in (bundle.genome_table, bundle.transcript_table):
        assert not any(line.endswith("\t") for line in table.read_text().splitlines())
    assert bundle.genome_table.read_text().splitlines()[1].endswith('\t""')


def test_reuse_is_offline_and_tampering_is_caught(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    reopened = open_viral_references(species="human", cache_root=str(tmp_path))
    assert reopened == bundle

    bundle.protein_fasta.write_text(bundle.protein_fasta.read_text() + ">EXTRA\nMM\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        open_viral_references(species="human", cache_root=str(tmp_path))
    with pytest.raises(FileNotFoundError, match="rebuild them with force=True"):
        download_viral_references(species="human", cache_root=str(tmp_path))


def test_missing_bundle_and_unknown_collection_explain_themselves(viral_panel, tmp_path):
    with pytest.raises(FileNotFoundError, match="are not installed"):
        open_viral_references(species="human", cache_root=str(tmp_path))
    with pytest.raises(ValueError, match="Unknown viral collection"):
        open_viral_references(species="human", collection="gdc", cache_root=str(tmp_path))
    with pytest.raises(ValueError, match="is a different target"):
        open_viral_references(species="mouse", cache_root=str(tmp_path))


def test_index_selects_whole_genomes_with_their_transcripts(viral_panel, tmp_path):
    _install(tmp_path)
    index = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert isinstance(index, ViralReferenceFiles)
    assert (index.genome_count, index.transcript_count, index.protein_count) == (2, 3, 3)
    assert index.subsets == () and index.contexts == (CONTEXT,)
    genomes = {line.split()[0][1:] for line in index.genome_fasta.read_text().splitlines()
               if line.startswith(">")}
    assert genomes == {"VTEST001.1", "VTEST002.1"}
    assert {line.split("\t")[0] for line in index.annotation_gtf.read_text().splitlines()
            if not line.startswith("#")} == genomes
    assert not any(line.endswith("\t") for line in index.genome_table.read_text().splitlines())

    assert select_viral_references(contexts=[CONTEXT], species="human",
                                  cache_root=str(tmp_path)).path == index.path
    with pytest.raises(ValueError, match="No viral references match"):
        select_viral_references(contexts=[CONTEXT], exclude_contexts=[CONTEXT],
                               species="human", cache_root=str(tmp_path))


def test_a_damaged_index_is_rebuilt_from_the_local_bundle(viral_panel, tmp_path):
    _install(tmp_path)
    index = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    index.protein_fasta.write_text("")
    rebuilt = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert rebuilt == index and rebuilt.protein_count == 3
    assert json.loads((rebuilt.path / "metadata.json").read_text())["genome_ids"] == [
        "VTEST001.1", "VTEST002.1"]


def test_transcript_blocks_map_intervals_onto_both_strands():
    plus = ((11, 20, "+"), (31, 40, "+"))
    assert transcript_blocks(plus, 0, 20) == ((11, 20, "+", 0), (31, 40, "+", 10))
    assert transcript_blocks(plus, 5, 10) == ((16, 20, "+", 0), (31, 35, "+", 5))
    minus = ((31, 40, "-"), (11, 20, "-"))
    assert transcript_blocks(minus, 0, 20) == ((31, 40, "-", 0), (11, 20, "-", 10))
    assert transcript_blocks(minus, 5, 10) == ((31, 35, "-", 0), (16, 20, "-", 5))
    assert transcript_blocks(plus, 0, 0) == ()


def test_viral_genome_opens_for_downstream_tools(viral_panel, tmp_path):
    pytest.importorskip("pyensembl")
    bundle = _install(tmp_path)
    genome = bundle.pyensembl_genome(cache_root=str(tmp_path))
    # The fourth record annotates no CDS, so the annotation places nothing on it.
    assert sorted(genome.contigs()) == ["VTEST001.1", "VTEST002.1", "VTEST003.1"]

    transcript = genome.transcript_by_id("VTEST001.1_cds1")
    assert transcript.contig == "VTEST001.1" and transcript.strand == "+"
    assert transcript.complete
    assert transcript.sequence == "ATG" + "GCT" * 8 + "TAA"
    assert str(transcript.coding_sequence) == "ATG" + "GCT" * 8 + "TAA"
    assert transcript.protein_sequence == "MAAAAAAAA"
    assert genome.transcript_by_id("VTEST002.1_cds1").strand == "-"
    assert [exon.start for exon in genome.transcript_by_id("VTEST002.1_cds1").exons] == [61, 11]


def test_viral_genome_annotates_variants_for_downstream_tools(viral_panel, tmp_path):
    pytest.importorskip("varcode")
    from varcode import Variant

    genome = _install(tmp_path).pyensembl_genome(cache_root=str(tmp_path))
    # VTEST001.1 codes ATG GCT... from base 1 on the plus strand.
    effect = Variant(contig="VTEST001.1", start=4, ref="G", alt="A",
                     genome=genome).effects().top_priority_effect()
    assert effect.short_description == "p.A2T"
    assert effect.transcript_id == "VTEST001.1_cds1"


def test_summary_reads_the_installed_layers_without_writing(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    whole = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert whole["collection"] == "gdc-viruses" and whole["selection"]["selected"] is False
    assert (whole["genome_count"], whole["transcript_count"]) == (4, 4)
    assert whole["complete_transcript_count"] == 3  # One CDS is 5'-partial.
    # A genome whose source record annotates no CDS is reported, not hidden.
    assert whole["genomes_without_transcripts"] == ["TV4"]
    assert [layer["name"] for layer in whole["layers"]] == list(VIRAL_LAYERS)

    selected = viral_reference_summary(subsets=None, contexts=[CONTEXT], species="human",
                                       cache_root=str(tmp_path))
    assert selected["genome_count"] == 2 and selected["layers"][0]["records"] == 2
    assert selected["selection"]["contexts"] == [CONTEXT]
    # Summarizing never materializes a derived index.
    assert not (bundle.path.parent.parent.parent / "indexes").exists()
