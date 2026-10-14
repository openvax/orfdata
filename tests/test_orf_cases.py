"""Real Ensembl 115 ORF cases: non-AUG starts, recoding, frameshifts, bicistronic and circular RNAs, products.

Transcript slices come from the pinned Ensembl 115 files in tests/data/references/ensembl115.json;
isoform sequences and literature facts from tests/data/references/uniprot_literature.json.
"""

import json
from pathlib import Path

import pytest
from pyensembl import Genome

from orfdata import (
    IRES,
    ORFEvidence,
    RNAFeature,
    SourceRNA,
    TranscriptORF,
    annotated_orfs,
    extra_orfs,
    products_from_uniprot,
    rfam_ires,
)

DATA = Path(__file__).parent / "data" / "references"
ENSEMBL = json.loads((DATA / "ensembl115.json").read_text())
SOURCES = json.loads((DATA / "uniprot_literature.json").read_text())
ISOFORMS = SOURCES["isoforms"]


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """A real PyEnsembl Genome over the Ensembl 115 rows of the fixture transcripts."""
    path = tmp_path_factory.mktemp("ensembl115")
    rows, genes, cdna, pep = [], {}, [], []
    for t in ENSEMBL["transcripts"]:
        gene = (f'gene_id "{t["gene_id"]}"; gene_version "{t["gene_version"]}"; gene_name "{t["gene_name"]}"; '
                f'gene_biotype "{t["gene_biotype"]}";')
        attrs = gene + f' transcript_id "{t["id"]}"; transcript_version "{t["version"]}"; transcript_biotype "{t["biotype"]}";'
        exons = t["exons"]
        start, end = exons[0][0], exons[-1][1]
        span = genes.get(t["gene_id"], (t["contig"], t["strand"], start, end, gene))
        genes[t["gene_id"]] = (t["contig"], t["strand"], min(start, span[2]), max(end, span[3]), gene)
        line = f'{t["contig"]}\tensembl\t{{}}\t{{}}\t{{}}\t.\t{t["strand"]}\t.\t{attrs}{{}}\n'
        rows.append(line.format("transcript", start, end, ""))
        rows += [line.format("exon", a, b, f' exon_number "{i + 1}"; exon_id "{t["id"]}-{i}";') for i, (a, b) in enumerate(exons)]
        protein = f' protein_id "{t["protein_id"]}"; protein_version "{t["protein_version"]}";' if "protein_id" in t else ""
        rows += [line.format("CDS", a, b, protein) for a, b in t["cds"]]
        rows += [line.format(feature, a, b, "") for feature in ("start_codon", "stop_codon") for a, b in t[feature]]
        cdna.append(f'>{t["id"]}.{t["version"]}\n{t["sequence"]}\n')
        if "protein_id" in t:
            pep.append(f'>{t["protein_id"]}.{t["protein_version"]}\n{t["protein"]}\n')
    rows += [f"{c}\tensembl\tgene\t{a}\t{b}\t.\t{s}\t.\t{attrs}\n" for c, s, a, b, attrs in genes.values()]
    (path / "slice.gtf").write_text("".join(rows))
    (path / "cdna.fa").write_text("".join(cdna))
    (path / "pep.fa").write_text("".join(pep))
    genome = Genome("GRCh38", "ensembl", "115", gtf_path_or_url=str(path / "slice.gtf"),
                    transcript_fasta_paths_or_urls=[str(path / "cdna.fa")],
                    protein_fasta_paths_or_urls=[str(path / "pep.fa")], cache_directory_path=str(path))
    genome.species = "homo_sapiens"
    genome.index()
    yield genome
    genome.close()


def test_selenocysteine_is_a_decoded_stop_codon(genome):
    gpx1 = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000419783"))
    assert (gpx1.start, gpx1.end, gpx1.start_codon, gpx1.stop_codon) == (75, 687, "ATG", "TAG")
    assert gpx1.decoded_stop_codons == ((48, "TGA", "U"),)
    assert gpx1.transcript.sequence[219:222] == "TGA"


def test_programmed_frameshifts_are_subcodon_introns(genome):
    # Ensembl skips one base for OAZ1's +1 frameshift and two for PEG10's -1 slippage.
    oaz1 = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000602676"))
    assert (oaz1.start, oaz1.end, len(oaz1.protein_sequence)) == (79, 766, 228)
    assert oaz1.subcodon_introns == ((283 - 79, 2271443, 2271443, 1),)
    peg10 = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000612748"))
    assert peg10.subcodon_introns == ((1396 - 217, 94664508, 94664509, -1),)
    assert len(peg10.coding_sequence) == 3 * (len(peg10.protein_sequence) + 1)


@pytest.mark.parametrize("isoform,transcript_id,codon", (
    ("P78344-1", "ENST00000339995", "GTG"), ("P28347-1", "ENST00000527636", "ATT"), ("P01106-2", "ENST00000621592", "CTG")))
def test_non_aug_starts_are_complete_and_named(genome, isoform, transcript_id, codon):
    transcript = genome.transcript_by_id(transcript_id)
    annotated = TranscriptORF.from_annotation(transcript)
    # Ensembl has no start_codon row for GTG or ATT starts; the CDS still starts inside the RNA.
    assert annotated.start_complete and annotated.start_codon == codon
    assert annotated.protein_sequence == ISOFORMS[isoform]["sequence"]
    mapping, = [m for m in ISOFORMS[isoform]["mappings"] if m["transcript"].startswith(transcript_id)]
    assert (annotated.start, annotated.end) == tuple(mapping["slice"])
    source = TranscriptORF.at(transcript, annotated.start, orf_id=isoform, source="UniProt 2026_03",
                              start_codon=codon, protein=ISOFORMS[isoform]["sequence"])
    assert (source.end, source.cds_relation, source.start_codon) == (annotated.end, "annotated_cds", codon)
    if codon != "CTG":  # Outside the standard code's initiation codons: the source must declare it.
        with pytest.raises(ValueError, match="declare start_codon"):
            TranscriptORF.at(transcript, annotated.start, orf_id=isoform, source="UniProt 2026_03")


def test_cug_initiated_isoforms_are_n_terminal_extensions(genome):
    fgf2 = genome.transcript_by_id("ENST00000644866")
    annotated = TranscriptORF.from_annotation(fgf2)
    assert (annotated.start, annotated.protein_sequence) == (344, ISOFORMS["P09038-2"]["sequence"])
    for isoform, start in (("P09038-1", 179), ("P09038-3", 221)):
        orf = TranscriptORF.at(fgf2, start, orf_id=isoform, source="UniProt 2026_03",
                               protein=ISOFORMS[isoform]["sequence"], orf_type="alternative_initiation")
        assert (orf.start_codon, orf.end, orf.cds_relation) == ("CTG", annotated.end, "N_terminal_extension")
        assert orf.protein_sequence.endswith(annotated.protein_sequence)
    with pytest.raises(ValueError, match="does not translate"):
        TranscriptORF.at(fgf2, 179, orf_id="wrong", source="test", protein=ISOFORMS["P09038-3"]["sequence"])


def test_readthrough_extension_from_a_declared_recoded_stop(genome):
    vegfa = genome.transcript_by_id("ENST00000672860")
    annotated = TranscriptORF.from_annotation(vegfa)
    paper = SOURCES["vegf_ax"]["paper"]
    # Eswarappa et al. 2014: the UGA stop is decoded as serine, adding a 22-residue extension.
    ax = TranscriptORF.at(vegfa, annotated.start, orf_id="VEGF-Ax", source=paper["citation"],
                          recoded={len(annotated.protein_sequence): "S"})
    extension = ax.protein_sequence[len(annotated.protein_sequence):]
    assert extension == SOURCES["vegf_ax"]["extension"] and len(extension) == paper["extension_length_aa"]
    assert ax.end - annotated.end == paper["inter_stop_nt"] + 3 and extension.endswith("SLTRKD")
    assert ax.decoded_stop_codons == ((395, "TGA", "S"),) and ax.cds_relation == "C_terminal_extension"
    # Ensembl's own readthrough model shows the read-through stop as X.
    model = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000425836"))
    assert (371, "TGA", "X") in model.decoded_stop_codons
    with pytest.raises(ValueError, match="not a stop codon"):
        TranscriptORF.at(vegfa, annotated.start, orf_id="bad", source="test", recoded={3: "S"})


def test_other_genes_annotated_cds_on_bicistronic_rnas(genome):
    mief1 = annotated_orfs(genome.transcript_by_id("ENST00000325301"))
    assert [(o.orf_id, o.source_gene_name, o.start, o.end, o.cds_relation) for o in mief1] == [
        ("ENSP00000490747.1", "MIURF", 176, 389, "upstream"), ("ENSP00000327124.2", "MIEF1", 487, 1879, "annotated_cds")]
    snrpn = annotated_orfs(genome.transcript_by_id("ENST00000390687"))
    assert [(o.source_gene_name, o.start, o.end, o.cds_relation, o.on_source_transcript) for o in snrpn] == [
        ("SNURF", 62, 278, "upstream", False), ("SNRPN", 466, 1189, "annotated_cds", True)]


def test_catalogue_categories_and_computed_relations_agree_on_mief1(genome):
    found = extra_orfs(genome)  # The bundled human catalogue, on these real transcripts.
    mief1 = found.for_transcript("ENST00000325301")
    assert [(o.orf_id, o.orf_type, o.cds_relation) for o in mief1] == [
        ("c22riboseqorf42", "uORF", "upstream"), ("c22norep69", "intORF", "internal")]


def test_circular_rna_orf_crosses_its_back_splice_junction(genome):
    facts = SOURCES["circZNF609"]
    linear = genome.transcript_by_id("ENST00000326648")
    circ = SourceRNA.from_transcript(linear, *facts["circ_transcript_slice"], id="circZNF609",
                                     source=facts["paper"]["citation"], topology="circular")
    assert (len(circ.sequence), circ.exons, circ.related_transcript_ids) == (874, ((64499293, 64500166),), ("ENST00000326648.5",))
    orf = TranscriptORF.at(circ, facts["orf"]["start"], orf_id="circZNF609-ORF", source="Legnini et al. 2017",
                           protein=facts["orf"]["protein"], evidence=ORFEvidence())
    assert orf.wraps_junction and orf.segments == ((127, 874), (0, 6))
    assert len(orf.coding_sequence) == 753 and len(orf.protein_sequence) == 250 and orf.stop_codon == "TGA"
    assert orf.coding_blocks == ((64499293, 64499298), (64499420, 64500166))
    assert orf.cds_relation == "no_cds" and orf.source_transcript_id == "circZNF609"
    assert json.loads(json.dumps(orf.to_dict()))["segments"] == [[127, 874], [0, 6]]
    linear_protein = TranscriptORF.from_annotation(linear).protein_sequence
    assert orf.protein_sequence[:249] == linear_protein[:249]
    stopless = SourceRNA("loop", "ATGAAA", circ.reference, "synthetic", topology="circular")
    with pytest.raises(ValueError, match="one traversal"):
        TranscriptORF.at(stopless, 0, orf_id="loop", source="synthetic")
    runoff = TranscriptORF.at(SourceRNA("runoff", "ATGAAACCC", circ.reference, "synthetic"), 0,
                              orf_id="runoff", source="synthetic")
    assert (runoff.stop_complete, runoff.stop_codon, runoff.protein_sequence) == (False, None, "MKP")


def test_precursor_products_map_to_codons_and_genome(genome):
    pomc = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000395826"))
    entry = SOURCES["pomc_entry"]
    assert pomc.protein_sequence == entry["sequence"]["value"] and (pomc.start, pomc.end) == (103, 907)
    products = products_from_uniprot(entry)
    assert len(products) == 12 and products[0].kind == "signal_peptide" and (products[0].start, products[0].end) == (0, 26)
    acth = next(p for p in products if (p.start, p.end) == (137, 176))
    assert acth.coding_slice(pomc) == (514, 631)
    assert acth.sequence(pomc) == pomc.protein_sequence[137:176] and acth.sequence(pomc).startswith("SYSMEHFRWG")
    blocks = acth.genomic_blocks(pomc)
    assert sum(b - a + 1 for a, b in blocks) == 3 * (176 - 137)
    assert all(any(a <= x and y <= b for a, b in pomc.coding_blocks) for x, y in blocks)


def test_rfam_ires_elements_located_exactly_on_transcripts(genome):
    myc = genome.transcript_by_id("ENST00000621592")
    element, = rfam_ires(myc)
    assert (element.feature_id.split(":")[0], element.kind, element.so_term) == ("RF00216", *IRES)
    assert (element.start, element.end) == (13, 315) and element.sequence == myc.sequence[13:315]
    assert "Activity in this transcript or condition is not established" in element.evidence
    assert element.overlaps("chr8", *element.genomic_blocks[0]) and not element.overlaps("8", 1, 10)
    # VEGFA's IRES lies inside the CTG-initiated CDS, spanning the AUG of the shorter isoform.
    vegfa = genome.transcript_by_id("ENST00000672860")
    element, = rfam_ires(vegfa)
    cds = TranscriptORF.from_annotation(vegfa)
    assert (element.start, element.end) == (742, 1051) and cds.start < element.start < element.end < cds.end
    assert vegfa.sequence[1036:1039] == "ATG"
    assert json.loads(json.dumps(element.to_dict()))["start"] == 742
    assert rfam_ires(genome.transcript_by_id("ENST00000419783")) == ()  # GPX1: none.


def test_rna_features_are_placed_only_at_one_exact_occurrence(genome):
    myc = genome.transcript_by_id("ENST00000621592")
    found = RNAFeature.find(myc, myc.sequence[13:315].replace("T", "U"), feature_id="x", kind=IRES[0],
                            so_term=IRES[1], source="test")
    assert (found.start, found.end) == (13, 315)
    with pytest.raises(ValueError, match="does not occur"):
        RNAFeature.find(myc, "ACGT" * 30, feature_id="absent", kind=IRES[0], so_term=IRES[1], source="test")
    with pytest.raises(ValueError, match="more than once"):
        RNAFeature.find(myc, "GG", feature_id="ambiguous", kind=IRES[0], so_term=IRES[1], source="test")


def test_catalogue_models_can_declare_a_non_aug_start(genome):
    from dataclasses import replace

    from orfdata import ORFModel, ORFRecord

    tead1 = genome.transcript_by_id("ENST00000527636")
    cds = TranscriptORF.from_annotation(tead1)
    model = ORFModel("P28347-1", tead1.id, tead1.gene_id, "TEAD1", tead1.contig, tead1.strand, cds.coding_blocks,
                     cds.protein_sequence, "annotated_CDS", "GRCh38", 9606, "homo_sapiens", source="UniProt",
                     start_codon="ATT")
    bound = ORFRecord.from_model(model).bind(tead1)
    assert (bound.start, bound.end, bound.start_codon) == (cds.start, cds.end, "ATT")
    assert ORFRecord.from_model(model).to_dict()["start_codon"] == "ATT"
    with pytest.raises(ValueError, match="initiation codon ATT"):
        ORFRecord.from_model(replace(model, start_codon=None)).bind(tead1)
    with pytest.raises(ValueError, match="declared start codon GTG"):
        ORFRecord.from_model(replace(model, start_codon="GTG")).bind(tead1)


def test_fragments_without_a_known_start(genome):
    pomc = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000395826"))
    # Reading from inside the CDS without claiming a start: no forced methionine.
    fragment = TranscriptORF.at(pomc.transcript, pomc.start + 3 * 137, orf_id="ACTH-fragment", source="test",
                                start_complete=False)
    assert fragment.protein_sequence == pomc.protein_sequence[137:]
    assert (fragment.start_complete, fragment.start_codon, fragment.end) == (False, None, pomc.end)
    assert fragment.decoded_stop_codons == ()  # The reading frame is the source's assertion; no recoding claimed.
