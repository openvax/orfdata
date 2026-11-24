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


def test_gtf_attribute_parser_round_trips_and_handles_quotes_and_strict_mode():
    from orfdata.gtf import attributes, parse_attribute_pairs, parse_attributes

    pairs = (
        ("gene_id", "ENSG000001"),
        ("gene_name", "HLA-A; alpha chain"),
        ("note", "O'Brien's \"quoted\" text"),
        ("transcript_id", "ENST000001"),
    )
    rendered = attributes(pairs)
    parsed = parse_attributes(rendered)
    assert parsed == {
        "gene_id": "ENSG000001",
        "gene_name": "HLA-A; alpha chain",
        "note": "O'Brien's 'quoted' text",
        "transcript_id": "ENST000001",
    }
    # Semicolon inside quotes, missing trailing semicolon, no space after semicolon, unquoted int:
    tricky = 'gene_id "G1";gene_name "A; B"; exon_number 2; transcript_id "T1"'
    assert parse_attribute_pairs(tricky) == (
        ("gene_id", "G1"),
        ("gene_name", "A; B"),
        ("exon_number", "2"),
        ("transcript_id", "T1"),
    )
    # Repeated tag allowed in strict mode, conflicting non-repeated key rejected:
    strict_ok = 'gene_id "G1"; tag "basic"; tag "CCDS";'
    assert parse_attributes(strict_ok, strict=True, repeated_keys=("tag",)) == {
        "gene_id": "G1",
        "tag": "CCDS",
    }
    with pytest.raises(ValueError, match="Conflicting GTF attribute: gene_name"):
        parse_attributes('gene_name "A"; gene_name "B";', strict=True, repeated_keys=("tag",))


def test_all_gtf_callers_handle_quoted_attributes_and_malformed_rows_identically(tmp_path):
    """Every GTF consumer parses the same quoted attribute value and rejects the same malformed row."""
    import gzip
    import hashlib
    import json
    from dataclasses import replace
    from pathlib import Path

    from orfdata import DatasetManifest, ReferenceContext
    from orfdata.build import load_transcript_reference
    from orfdata.ensembl import _annotation as ensembl_annotation
    from orfdata.gtf import (
        read_exons,
        read_gtf_rows,
        read_transcript_exons,
        read_transcript_loci,
    )
    from orfdata.sources.riborf_native import _annotations as riborf_native_annotations
    from orfdata.sources.wang import _annotation as wang_annotation

    quoted_name = "HLA-A; alpha chain"
    # Note: quoted attribute with space & semicolon AND final attribute without trailing ';'.
    attrs_tx = (
        f'gene_id "ENSG00000169800"; gene_version "1"; gene_name "{quoted_name}"; '
        'transcript_id "ENST00000303766"; transcript_version "1"; '
        'transcript_biotype "protein_coding"'
    )
    attrs_exon = f'{attrs_tx}; exon_number "1"'
    attrs_cds = f'{attrs_tx}; protein_id "ENSP00000307155"; protein_version "1"'
    valid_text = (
        "#!genome-build GRCh38.p5\n"
        "#!genome-version GRCh38\n"
        "#!genome-build-accession GCA_000001405.20\n"
        f"1\tensembl\ttranscript\t10\t30\t.\t+\t.\t{attrs_tx}\n"
        f"1\tensembl\texon\t10\t30\t.\t+\t.\t{attrs_exon}\n"
        f"1\tensembl\tCDS\t10\t27\t.\t+\t0\t{attrs_cds}\n"
    )
    malformed_text = valid_text + "1\tensembl\tgene\t10\n"

    valid_plain = tmp_path / "valid.gtf"
    valid_gz = tmp_path / "valid.gtf.gz"
    bad_plain = tmp_path / "malformed.gtf"
    bad_gz = tmp_path / "malformed.gtf.gz"
    valid_plain.write_text(valid_text)
    valid_gz.write_bytes(gzip.compress(valid_text.encode(), mtime=0))
    bad_plain.write_text(malformed_text)
    bad_gz.write_bytes(gzip.compress(malformed_text.encode(), mtime=0))

    # 1. Shared reader: plain & gzipped agree on headers and quoted attribute value.
    for path in (valid_plain, valid_gz):
        headers = {}
        rows = list(read_gtf_rows(path, headers=headers))
        assert len(rows) == 3
        assert headers["genome-build"] == "GRCh38.p5"
        assert all(row[5]["gene_name"] == quoted_name for row in rows)
        assert all(row[5]["transcript_biotype"] == "protein_coding" for row in rows)
        assert read_exons(path) == {"ENST00000303766": (((10, 30),), "+")}

        tx_exons = read_transcript_exons(path, transcript_ids={"ENST00000303766"})
        assert tx_exons["ENST00000303766"]["gene_name"] == quoted_name
        assert tx_exons["ENST00000303766"]["exons"] == ((10, 30),)

        loci = read_transcript_loci(path, genome_build="GRCh38")
        assert loci["ENST00000303766.1"][0]["gene_name"] == quoted_name
        assert riborf_native_annotations(path, "GRCh38") == loci

    # 2. Ensembl annotation reader parses the quoted gene_name identically.
    ens_ref = ReferenceContext(
        9606, "human", "GRCh38", "Ensembl", "115", assembly_accession="GCA_000001405.20"
    )
    ens_tx, ens_products = ensembl_annotation(
        valid_gz, hashlib.sha256(valid_gz.read_bytes()).hexdigest(), ens_ref
    )
    (tx_entry,) = ens_tx["ENST00000303766.1"]
    assert tx_entry.name == quoted_name
    assert tx_entry.biotype == "protein_coding"
    assert "ENSP00000307155.1" in ens_products

    # 3. build.load_transcript_reference reads the same GTF (plain or gz) cleanly.
    cdna_gz = tmp_path / "cdna.fa.gz"
    cdna_gz.write_bytes(gzip.compress(b">ENST00000303766.1\n" + b"A" * 21 + b"\n", mtime=0))
    built_ref = load_transcript_reference(valid_gz, cdna_gz, {"ENST00000303766"})
    assert built_ref["ENST00000303766"]["transcript_id"] == "ENST00000303766.1"
    assert built_ref["ENST00000303766"]["exons"] == [(10, 30)]

    # 4. Wang annotation reader parses the same GTF cleanly.
    wang_root = Path(__file__).parent / "data" / "wang"
    wang_manifest = DatasetManifest.from_dict(json.loads((wang_root / "manifest.json").read_text()))
    valid_sha = hashlib.sha256(valid_gz.read_bytes()).hexdigest()
    wang_sources = tuple(
        replace(s, sha256=valid_sha, size=valid_gz.stat().st_size)
        if s.asset_id == "annotation.gtf.gz"
        else s
        for s in wang_manifest.sources
    )
    wang_ref = replace(wang_manifest.references[0], source_sha256=valid_sha)
    wang_cov = dict(wang_manifest.coverage)
    wang_cov["annotation_counts"] = {"genes": 1, "transcripts": 1, "proteins": 1}
    wang_m = replace(wang_manifest, references=(wang_ref,), sources=wang_sources, coverage=wang_cov)
    wang_inputs = {s.asset_id: wang_root / s.asset_id for s in wang_manifest.sources}
    wang_inputs["annotation.gtf.gz"] = valid_gz
    _, genes, transcripts, proteins, _, _ = wang_annotation(wang_inputs, wang_m)
    assert genes == {"ENSG00000169800"}
    assert transcripts == {"ENST00000303766": "ENSG00000169800"}
    assert proteins == {"ENSP00000307155": ("ENST00000303766", "ENSG00000169800")}

    # 5. Every single caller raises ValueError naming line 7 on the malformed 4-field row.
    for path in (bad_plain, bad_gz):
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            list(read_gtf_rows(path))
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            read_exons(path)
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            read_transcript_exons(path)
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            read_transcript_loci(path, genome_build="GRCh38")
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            riborf_native_annotations(path, "GRCh38")
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            load_transcript_reference(path, cdna_gz, {"ENST00000303766"})
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            ensembl_annotation(path, hashlib.sha256(path.read_bytes()).hexdigest(), ens_ref)
        bad_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        bad_sources = tuple(
            replace(s, sha256=bad_sha, size=path.stat().st_size)
            if s.asset_id == "annotation.gtf.gz"
            else s
            for s in wang_manifest.sources
        )
        bad_m = replace(
            wang_m,
            references=(replace(wang_manifest.references[0], source_sha256=bad_sha),),
            sources=bad_sources,
        )
        bad_inputs = dict(wang_inputs, **{"annotation.gtf.gz": path})
        with pytest.raises(ValueError, match=r"nine tab-separated fields.*:7"):
            wang_annotation(bad_inputs, bad_m)


def test_promoted_model_measurements_and_gtf_helpers_are_public_and_documented():
    import orfdata.gtf as gtf_mod
    import orfdata.measurements as measurements_mod
    import orfdata.model as model_mod
    import orfdata.sources.riborf as riborf_mod

    assert not hasattr(model_mod, "_species")
    for name in (
        "normalize_contig",
        "exon_layout",
        "rna_layout",
        "layout_offset",
        "layout_blocks",
        "merge_genomic_ranges",
        "same_assembly",
        "match_transcript_id",
        "transcript_id_forms",
        "versioned_id",
    ):
        obj = getattr(model_mod, name)
        assert callable(obj) and obj.__doc__
    for name in (
        "Record",
        "require_identifier",
        "freeze_metadata",
        "require_finite_number",
        "to_plain_data",
    ):
        obj = getattr(measurements_mod, name)
        assert obj.__doc__
    for name in (
        "parse_attribute_pairs",
        "parse_attributes",
        "open_gtf",
        "read_gtf_rows",
        "read_exons",
        "read_transcript_exons",
        "read_transcript_loci",
        "transcript_interval_blocks",
    ):
        obj = getattr(gtf_mod, name)
        assert callable(obj) and obj.__doc__
    for name in ("read_riborf_rows", "parse_riborf_rows"):
        obj = getattr(riborf_mod, name)
        assert callable(obj) and obj.__doc__
