"""The GenBank reader keeps what the file says, including what it says it does not know."""

import pytest

from orfdata.genbank import Complement, Join, Order, Range, parse, parse_location

RECORD = """LOCUS       TEST0001                 120 bp    DNA     circular VRL 01-JAN-2020
DEFINITION  Synthetic test virus, complete genome.
ACCESSION   TEST0001
VERSION     TEST0001.1
KEYWORDS    .
SOURCE      Synthetic test virus
  ORGANISM  Synthetic test virus
            Viruses; Synthetic.
REFERENCE   1  (bases 1 to 120)
  AUTHORS   Nobody.
FEATURES             Location/Qualifiers
     source          1..120
                     /organism="Synthetic test virus"
                     /mol_type="genomic DNA"
     CDS             complement(join(11..40,
                     61..90))
                     /gene="rev"
                     /product="reverse strand protein,
                     spliced"
                     /note="says \\"maybe\\" here"
                     /codon_start=1
                     /transl_table=11
                     /translation="MKTLLSAIEQ
                     RSTVW"
     CDS             <1..>60
                     /gene="partial"
     CDS             join(100..120,1..30)
                     /gene="wraps"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
""".replace('\\"', '""')


def test_locations_keep_their_structure_order_and_markers():
    assert parse_location("1..30") == Range(1, 30)
    assert parse_location("<1..>60") == Range(1, 60, before=True, after=True)
    assert parse_location("complement(11..40)") == Complement(Range(11, 40))
    # Order matters: a join is read in the order written, not sorted.
    assert parse_location("join(100..120,1..30)") == Join([Range(100, 120), Range(1, 30)])
    assert parse_location("complement(join(11..40,61..90))") == Complement(
        Join([Range(11, 40), Range(61, 90)]))
    assert parse_location("order(>104408,104493..105098)") == Order(
        [Range(104408, 104408, after=True), Range(104493, 105098)])
    assert parse_location("join(1..10,complement(20..30))") == Join(
        [Range(1, 10), Complement(Range(20, 30))])
    # Whitespace from a wrapped location line is not part of it.
    assert parse_location("join(11..40,\n                     61..90)") == Join(
        [Range(11, 40), Range(61, 90)])


@pytest.mark.parametrize("text", [
    "123^124",          # a site between bases
    "(1.5)..100",       # an uncertain position
    "J00194.1:10..20",  # another record's coordinates
    "join(1..10",       # unbalanced
    "join(1..10,)",     # empty part
    ">1..100",          # a marker ORFdata will not guess at
])
def test_locations_outside_the_supported_subset_are_refused(text):
    with pytest.raises(ValueError, match="Unsupported|Unbalanced"):
        parse_location(text)


def test_a_record_reads_back_whole():
    record = parse(RECORD)
    assert record.accession_version == "TEST0001.1"
    assert record.circular is True
    assert record.sequence == "ACGT" * 30
    assert [feature.kind for feature in record.features] == ["source", "CDS", "CDS", "CDS"]

    spliced, partial, wrapping = record.features[1:]
    assert spliced.location == Complement(Join([Range(11, 40), Range(61, 90)]))
    assert partial.location == Range(1, 60, before=True, after=True)
    assert wrapping.location == Join([Range(100, 120), Range(1, 30)])


def test_qualifiers_keep_their_order_values_and_wrapping():
    qualifiers = dict(parse(RECORD).features[1].qualifiers)
    # A wrapped sequence closes up; wrapped prose keeps the break for the reader.
    assert qualifiers["translation"] == "MKTLLSAIEQRSTVW"
    assert qualifiers["product"] == "reverse strand protein,\nspliced"
    assert qualifiers["note"] == 'says "maybe" here'
    assert qualifiers["codon_start"] == "1" and qualifiers["transl_table"] == "11"
    assert dict(parse(RECORD).features[0].qualifiers)["organism"] == "Synthetic test virus"


def test_a_linear_record_without_features_still_reads():
    record = parse("LOCUS       X 4 bp DNA linear VRL 01-JAN-2020\n"
                   "VERSION     X.1\nORIGIN\n        1 acgt\n//\n")
    assert (record.circular, record.sequence, record.features) == (False, "ACGT", ())


def test_a_record_without_a_version_is_an_error():
    with pytest.raises(ValueError, match="no VERSION line"):
        parse("LOCUS       X 4 bp DNA linear VRL 01-JAN-2020\nORIGIN\n        1 acgt\n//\n")


WRAPPED = '''LOCUS       TEST0002                  12 bp    DNA     linear   VRL 01-JAN-2020
VERSION     TEST0002.1
FEATURES             Location/Qualifiers
     CDS             1..12
                     /note="dosing is 5 units
                     /day in the first week"
                     /gene="A"
                     /db_xref="GO:0003677"
                     /pseudo
ORIGIN
        1 acgtacgtacgt
//
'''


def test_a_wrapped_value_keeps_a_line_that_starts_like_a_qualifier():
    qualifiers = parse(WRAPPED).features[0].qualifiers
    values = dict(qualifiers)
    # The "/day" line continues the note; it does not start a qualifier.
    assert values["note"] == "dosing is 5 units\n/day in the first week"
    assert [name for name, _ in qualifiers] == ["note", "gene", "db_xref", "pseudo"]
    assert values["gene"] == "A" and values["db_xref"] == "GO:0003677"
    assert values["pseudo"] is None


def test_a_truncated_record_is_caught_by_what_it_says_about_itself():
    # The LOCUS line promises 120 bases; a cut-short ORIGIN is a damaged record.
    truncated = RECORD.replace(
        "       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n", "")
    with pytest.raises(ValueError, match="declares 120 bases but carries 60"):
        parse(truncated)
    # A record with no ORIGIN at all, such as a CONTIG entry, is not a truncation.
    assert parse(RECORD.split("ORIGIN")[0] + "//\n").sequence == ""


def test_a_feature_whose_location_starts_on_the_next_line_still_reads():
    text = RECORD.replace("     CDS             complement(join(11..40,\n                     61..90))",
                          "     CDS\n                     complement(join(11..40,61..90))")
    feature = parse(text).features[1]
    assert feature.kind == "CDS"
    assert feature.location == Complement(Join([Range(11, 40), Range(61, 90)]))


def test_the_locus_name_and_topology_are_kept():
    record = parse(RECORD)
    assert record.name == "TEST0001" and record.circular is True
    assert parse(RECORD.replace(" circular ", " linear   ")).circular is False


def test_a_feature_outside_the_sequence_is_refused():
    # A record whose coordinates do not fit its own bases is damaged, not unusual.
    text = RECORD.replace("     CDS             <1..>60", "     CDS             <1..>600")
    with pytest.raises(ValueError, match="places its CDS feature at 600, outside its 120 bases"):
        parse(text)
    # Without a sequence there is nothing to check against.
    assert parse(text.split("ORIGIN")[0] + "//\n").features[2].location == Range(
        1, 600, before=True, after=True)
