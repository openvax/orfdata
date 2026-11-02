"""The readers return exactly what the writers wrote."""

from orfdata.fasta import read_fasta, write_fasta
from orfdata.gtf import read_exons, transcript_lines


def test_read_fasta_returns_what_write_fasta_wrote(tmp_path):
    path = tmp_path / "x.fasta"
    long = "ACGT" * 50  # Wraps over several lines.
    with path.open("w") as stream:
        write_fasta(stream, "a.1", long, description="first record")
        write_fasta(stream, "b.1", "MKV")
    assert read_fasta(path) == {"a.1": long, "b.1": "MKV"}
    assert list(read_fasta(path)) == ["a.1", "b.1"]


def test_read_exons_gives_transcript_order_for_both_strands(tmp_path):
    path = tmp_path / "x.gtf"
    plus = ((10, 20, "+"), (30, 40, "+"))
    minus = ((30, 40, "-"), (10, 20, "-"))  # Transcript order runs high to low.
    wrapped = ((90, 100, "+"), (1, 5, "+"))  # Joined across a circular origin.
    rows = []
    for name, exons in (("p", plus), ("m", minus), ("w", wrapped)):
        rows += transcript_lines("g", exons, gene_id=name, gene_name=name,
                                 transcript_id=f"t_{name}", protein_id=f"p_{name}")
    path.write_text("# header\n" + "".join(rows))
    assert read_exons(path) == {
        "t_p": (((10, 20), (30, 40)), "+"),
        "t_m": (((30, 40), (10, 20)), "-"),
        "t_w": (((90, 100), (1, 5)), "+"),
    }


def test_unreadable_files_are_errors_not_silent_losses(tmp_path):
    import pytest

    path = tmp_path / "x.fasta"
    path.write_text(">a\nAC\n>a\nGT\n")
    with pytest.raises(ValueError, match="repeated"):
        read_fasta(path)
    path.write_text(">\nAC\n")
    with pytest.raises(ValueError, match="without an identifier"):
        read_fasta(path)
    gtf = tmp_path / "x.gtf"
    gtf.write_text('g\torfdata\texon\t1\t5\t.\t+\t.\ttranscript_id "t";\n')
    with pytest.raises(ValueError, match="exon_number"):
        read_exons(gtf)


def test_inconsistent_input_is_refused_with_its_location(tmp_path):
    import pytest

    path = tmp_path / "x.fasta"
    path.write_text("AC\n>a\nGT\n")
    with pytest.raises(ValueError, match="before the first header"):
        read_fasta(path)
    gtf = tmp_path / "x.gtf"
    row = 'g\torfdata\texon\t{}\t{}\t.\t{}\t.\ttranscript_id "t"; exon_number "{}";\n'
    gtf.write_text(row.format(1, 5, "+", 1) + row.format(9, 12, "+", 1))
    with pytest.raises(ValueError, match="repeats an exon number"):
        read_exons(gtf)
    gtf.write_text(row.format(1, "NA", "+", 1))
    with pytest.raises(ValueError, match="x.gtf:1"):
        read_exons(gtf)


def test_a_row_that_lost_its_tabs_is_refused_not_skipped(tmp_path):
    import pytest

    gtf = tmp_path / "x.gtf"
    gtf.write_text("# header\n\ng orfdata exon 1 5 . + . transcript_id \"t\"; exon_number \"1\";\n")
    with pytest.raises(ValueError, match="x.gtf:3"):
        read_exons(gtf)


def test_upper_case_is_applied_while_reading(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(">a\nacgt\nnn\n")
    assert read_fasta(path, upper=True) == {"a": "ACGTNN"}
    assert read_fasta(path) == {"a": "acgtnn"}


def test_attribute_values_with_quotes_and_semicolons_survive_a_round_trip(tmp_path):
    from orfdata.gtf import attributes, parse_attributes

    column = attributes([("gene_name", 'a "b"; c\td'), ("empty", ""), ("gene_id", "g1")])
    assert parse_attributes(column) == {"gene_name": "a 'b'; c d", "gene_id": "g1"}
    # The row stays readable as a whole, and a PyEnsembl-style parser agrees.
    path = tmp_path / "x.gtf"
    path.write_text(
        "".join(transcript_lines("g", ((1, 9, "+"),), gene_id="g1", gene_name='say "hi";',
                                 transcript_id="t", protein_id="p")))
    assert read_exons(path) == {"t": (((1, 9),), "+")}
    import gtfparse

    frame = gtfparse.read_gtf(str(path))
    assert set(frame["gene_name"].to_list()) == {"say 'hi';"}
