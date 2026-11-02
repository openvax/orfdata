"""A CDS a source annotates is placed as the source gives it, wherever on the genome it lies."""

import pytest

from orfdata.reference import ReferenceContext
from orfdata.rna import SourceRNA
from orfdata.source_cds import source_cds_orf

CONTEXT = ReferenceContext(taxon_id=1, species="test_virus", genome_build="G", annotation_name="t",
                           annotation_version="1", sequence_sha256="0" * 64)
FLIP = str.maketrans("ACGT", "TGCA")


def orf(genome, blocks, strand, protein, **options):
    return source_cds_orf(genome, blocks, strand, context=CONTEXT, label="G", accession="A.1",
                          protein_id="P.1", gene="g", source="test", protein=protein,
                          translation_table=1, start_complete=options.pop("start_complete", True),
                          **options)


def test_a_plain_cds_on_either_strand():
    coding = "ATGGCTAAATAA"
    genome = "CC" + coding + "CC"
    assert orf(genome, ((3, 14),), "+", "MAK").coding_blocks == ((3, 14),)
    flipped = "CC" + coding[::-1].translate(FLIP) + "CC"
    assert orf(flipped, ((3, 14),), "-", "MAK").coding_sequence == coding


def test_a_cds_continues_across_the_origin_of_a_circular_genome():
    coding = "ATGGCTAAATAA"
    genome = coding[6:] + "CCCCC" + coding[:6]   # the CDS runs off the end and restarts at base 1
    last = len(genome)
    found = orf(genome, ((last - 5, last), (1, 6)), "+", "MAK")
    assert found.coding_sequence == coding and sorted(found.coding_blocks) == [(1, 6), (last - 5, last)]
    # On the minus strand the wrap runs from the first base to the last.
    reverse = coding[::-1].translate(FLIP)
    genome = reverse[6:] + "CCCCC" + reverse[:6]
    assert orf(genome, ((1, 6), (last - 5, last)), "-", "MAK").coding_sequence == coding


def test_a_frameshift_that_skips_bases_is_recorded_as_a_subcodon_intron():
    # Bases 1..12 then 14..16: base 13 is skipped, as OAZ1's +1 frameshift skips one.
    genome = "ATGGCTAAAGCT" + "C" + "TAA" + "CC"
    found = orf(genome, ((1, 12), (14, 16)), "+", "MAKA")
    assert found.coding_sequence == "ATGGCTAAAGCTTAA" and found.protein_sequence == "MAKA"
    assert found.subcodon_introns == ((12, 13, 13, 1),)


def test_a_base_read_twice_is_kept_out_until_the_orf_can_record_it():
    # Bases 1..12 then 12..14 read base 12 twice, which the ORF's merged blocks would lose.
    genome = "ATGGCTAAAGCT" + "AA" + "CCCC"
    with pytest.raises(ValueError, match="read twice"):
        orf(genome, ((1, 12), (12, 14)), "+", "MAKA")


def test_a_cds_without_a_start_or_a_stop_keeps_what_the_source_established():
    five_prime_partial = orf("GCTAAATAA", ((1, 9),), "+", "AK*"[:2], start_complete=False)
    assert not five_prime_partial.start_complete and five_prime_partial.stop_complete
    stopless = orf("ATGGCTAAA", ((1, 9),), "+", "MAK")
    assert stopless.start_complete and not stopless.stop_complete


@pytest.mark.parametrize("blocks, strand, contig_length", [
    (((10, 30), (1, 9)), "+", None),         # Backwards, and no circular contig to wrap.
    (((10, 30), (1, 9)), "+", 40),           # Backwards, but not from the contig's last base.
    (((1, 20), (10, 30)), "+", None),        # Overlaps by more than one base.
    (((1, 10), (11, 30)), "-", None),        # Forwards on the minus strand.
    (((1, 10), (12, 40)), "+", 30),          # Beyond the end of the contig.
])
def test_blocks_that_do_not_follow_the_strand_are_refused(blocks, strand, contig_length):
    length = sum(b - a + 1 for a, b in blocks)
    with pytest.raises(ValueError, match="Exons"):
        SourceRNA("r", "A" * length, CONTEXT, "test", contig="G", strand=strand, exons=blocks,
                  contig_length=contig_length)
