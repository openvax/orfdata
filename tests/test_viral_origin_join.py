"""A gene split by a circular genome's origin is one gene, when the join proves it."""

import pytest

from orfdata.viral_annotation import coding_feature, extract, join_across_origin, parse_genbank

# A 60-base circular genome whose only gene starts at 49 and ends at 9, across the origin.
ORF = "ATG" + "GCT" * 5 + "TAA"
GENOME = list("T" * 60)
for index, base in enumerate(ORF[:12]):
    GENOME[48 + index] = base
for index, base in enumerate(ORF[12:]):
    GENOME[index] = base
GENOME = "".join(GENOME)
PROTEIN = "MAAAAA"


def _half(ordinal, exons, *, partial, translation):
    return coding_feature(ordinal, exons, extract(GENOME, exons),
                          {"gene": "P", "product": "polymerase", "translation": translation},
                          partial=partial)


def _halves(*, tail_translation=PROTEIN[:4], head_translation=PROTEIN[4:]):
    return (_half(1, ((49, 60, "+"),), partial=(False, True), translation=tail_translation),
            _half(2, ((1, 9, "+"),), partial=(True, False), translation=head_translation))


def test_the_two_halves_become_one_complete_gene():
    joined = join_across_origin(_halves(), GENOME)
    assert len(joined) == 1
    gene, = joined
    assert gene["location"] == "join(49..60,1..9)"
    assert gene["exons"] == ((49, 60, "+"), (1, 9, "+"))
    assert gene["nucleotides"] == ORF and gene["translation"] == PROTEIN
    assert gene["partial"] == "" and gene["complete"] is True
    assert gene["annotation_source"] == "genbank:origin_join"
    assert gene["ordinal"] == 1


def test_halves_that_do_not_read_as_one_gene_are_left_alone():
    # The published translations must meet end to end.
    untouched = join_across_origin(_halves(head_translation="WWW"), GENOME)
    assert len(untouched) == 2 and all(f["partial"] for f in untouched)

    # A partial pair that does not reach the record's edges is not an origin split.
    inner = (_half(1, ((40, 51, "+"),), partial=(False, True), translation=PROTEIN[:4]),
             _half(2, ((2, 10, "+"),), partial=(True, False), translation=PROTEIN[4:]))
    assert len(join_across_origin(inner, GENOME)) == 2

    # One half alone stays as the source wrote it.
    assert len(join_across_origin(_halves()[:1], GENOME)) == 1


def test_a_linear_record_is_never_joined():
    record, = parse_genbank(_record(circular=False).encode("utf-8"))
    assert not record.circular
    assert [f["partial"] for f in record.coding_features] == ["3_prime", "5_prime"]
    assert {f["annotation_source"] for f in record.coding_features} == {"genbank"}


def test_a_circular_record_is_joined_when_it_is_read():
    record, = parse_genbank(_record(circular=True).encode("utf-8"))
    assert record.circular
    gene, = record.coding_features
    assert gene["annotation_source"] == "genbank:origin_join"
    assert gene["translation"] == PROTEIN and gene["complete"] is True


def _record(*, circular):
    topology = "circular" if circular else "linear"
    lines = [f"LOCUS       TESTCIRC                  60 bp    DNA     {topology} VRL 01-JAN-2020",
             "VERSION     TESTCIRC.1",
             "FEATURES             Location/Qualifiers",
             "     CDS             49..>60",
             '                     /gene="P"',
             f'                     /translation="{PROTEIN[:4]}"',
             "     CDS             <1..9",
             '                     /gene="P"',
             f'                     /translation="{PROTEIN[4:]}"',
             "ORIGIN"]
    for offset in range(0, len(GENOME), 60):
        block = GENOME[offset:offset + 60].lower()
        lines.append(f"{offset + 1:>9} " + " ".join(block[i:i + 10] for i in range(0, len(block), 10)))
    lines.append("//")
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("translation", [PROTEIN[:4], "M"])
def test_the_join_needs_a_start_codon_and_a_clean_frame(translation):
    # Shifting the tail so the joined sequence no longer begins at a start codon.
    shifted = (_half(1, ((50, 60, "+"),), partial=(False, True), translation=translation),
               _half(2, ((1, 9, "+"),), partial=(True, False), translation=PROTEIN[4:]))
    assert len(join_across_origin(shifted, GENOME)) == 2


def test_a_record_that_names_nothing_still_describes_its_product():
    # HBV's polymerase has neither /gene nor /product, only a /note.
    from orfdata.viral_annotation import coding_feature, name_transcripts

    feature = coding_feature(1, ((49, 60, "+"), (1, 9, "+")), ORF,
                             {"note": "P gene product, put. DNA polymerase",
                              "translation": PROTEIN})
    assert feature["gene"] == ""
    assert feature["product"] == "P gene product, put. DNA polymerase"
    # The note describes the product; it is never used as a gene label.
    transcript, = name_transcripts("X.1", (feature,), unique_protein_ids=set())
    assert transcript["gene_name"] == "cds1"
    assert transcript["gene_id"] == "X.1_cds1_plus"
