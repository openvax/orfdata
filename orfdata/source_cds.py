"""A CDS a source annotates on a viral genome, placed as a SourceRNA and read as an ORF."""

from .model import TranscriptORF
from .rna import SourceRNA

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def source_cds_orf(genome, blocks, strand, *, context, label, accession, protein_id, gene, source,
                   protein, translation_table, start_complete, start_codon=None, phase=0):
    """The ORF of one source CDS on ``genome`` (named ``label``), exactly as the source gives it.

    ``blocks`` are the CDS's 1-based inclusive genomic blocks in reading order,
    which wrap the origin of a circular genome or share a base across a
    frameshift where the source says so. The ORF is read from the first base to
    the first stop, and must translate to the source ``protein``; otherwise this
    raises ``ValueError`` and the caller keeps the CDS out, with the reason.
    """
    coding = "".join(genome[a - 1:b] if strand == "+" else genome[a - 1:b].translate(_COMPLEMENT)[::-1]
                     for a, b in blocks)
    rna = SourceRNA(
        id=f"{accession}:CDS:{protein_id}", sequence=coding, reference=context, source=source,
        contig=label, strand=strand, exons=blocks, gene_name=gene, biotype="source_CDS_fragment",
        contig_length=len(genome))
    orf = TranscriptORF.at(
        rna, 0, orf_id=protein_id, source=source, protein=protein, translation_table=translation_table,
        orf_type="reference_CDS", start_codon=start_codon, start_complete=start_complete, phase=phase)
    if orf.coding_sequence != coding:
        raise ValueError("coding boundaries do not translate exactly")
    return orf
