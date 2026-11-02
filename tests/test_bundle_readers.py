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
