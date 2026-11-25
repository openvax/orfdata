"""Reference evidence must survive exports without becoming mutant evidence."""

import gzip
import hashlib
import io
import json
import tarfile
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib.resources import files
from pathlib import Path

import pytest

from orfdata import UORFCatalog, UORFPeptide, UORFRecord, load_uorf_catalog
from orfdata import build as builder
from orfdata.build import intervals
from orfdata.cli import run as _run


def run(argv=None, prog=None):
    return _run(argv, prog, catalog_loader=load_uorf_catalog)


@pytest.fixture(scope="module")
def catalog():
    return load_uorf_catalog()


def test_native_tpst1_observed_tail_does_not_support_fusion(catalog):
    native = catalog.get("c7riboseqorf56")
    assert native.protein_sequence == "MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV"
    assert native.transcript_id == "ENST00000304842"
    assert native.coding_blocks == ((66205483, 66205522), (66240325, 66240404))
    assert len(native.ribosome_studies) == 6
    (peptide,) = native.peptides
    assert (peptide.method, peptide.sequence, peptide.protein_intervals) == (
        "immunopeptidomics",
        "KIILSKISV",
        ((31, 39),),
    )
    assert peptide.spectrum_count == 28 and peptide.unique_in_source
    assert peptide.quality == "reported"  # No unimported supplementary tier claim.
    fusion = "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"
    assert fusion[:13] == native.protein_sequence[:13]
    assert intervals(fusion, peptide.sequence) == []
    start = native.map_genomic_position(66205483, genome_build="GRCh38")
    assert start["role"] == "start_codon" and start["amino_acid_position"] == 1
    assert not start["reference_ms_support"] and start["reference_peptides"] == []
    tail = native.map_genomic_position(66240380, genome_build="GRCh38")
    assert tail["amino_acid_position"] == 32 and tail["reference_amino_acid"] == "I"
    assert tail["reference_ms_support"]
    assert tail["mutant_translation_evidence"] == "not_assessed"


def test_proteomics_and_mixed_assay_sources_are_not_confused(catalog):
    native = catalog.get("c1riboseqorf94")
    assert native.gene_name == "CDKN2C"
    assert any(
        p.method == "proteomics"
        and p.sequence == "FPGLSPSPEAPTKPGK"
        and p.datasets == ("Zhu_et_al_2018",)
        for p in native.peptides
    )
    assert native.has_evidence("proteomics")
    assert not native.has_evidence("proteomics", require_reviewed=True)
    assert not native.has_evidence("proteomics", require_unique=True)
    kat6a = catalog.get("c8riboseqorf42")
    assert kat6a.has_evidence("ms_unclassified")
    assert not kat6a.has_evidence("proteomics")  # Actual proteomics reviews are poor.
    assert any(
        p.method == "proteomics" and p.quality == "reviewed_low_quality" for p in kat6a.peptides
    )
    assert any(p.method == "ms_unclassified" and p.quality == "reported" for p in kat6a.peptides)


def test_review_scores_mapping_ambiguity_and_counts_keep_their_scopes(catalog):
    bmpr2 = catalog.get("c2riboseqorf155")
    assert bmpr2.has_evidence("immunopeptidomics", require_reviewed=True)
    assert bmpr2.has_evidence("immunopeptidomics", require_unique=True)
    # Distinct observations cannot lend their review or mapping count to each other.
    assert not bmpr2.has_evidence("immunopeptidomics", require_reviewed=True, require_unique=True)
    reviewed = [p for p in bmpr2.peptides if p.quality == "reviewed_support"]
    assert reviewed and all(p.attribution == "sequence_compatible" for p in reviewed)
    assert all(p.mapping_count is None and p.spectrum_ids for p in reviewed)
    assert all(all(score >= 4 for _, score in p.review_ratings) for p in reviewed)
    assert any(p.quality == "reviewed_mixed" and not p.is_supporting for p in bmpr2.peptides)
    assert (
        "not RNA abundance"
        in next(s for s in catalog.metadata["sources"] if s["id"] == "author_hla")["count_scope"]
    )


def test_all_source_review_scores_and_peptide_positions_validate(catalog):
    for record in catalog.records:
        for peptide in record.peptides:
            if peptide.quality == "reviewed_support":
                assert all(score >= 4 for _, score in peptide.review_ratings)
            elif peptide.quality == "reviewed_low_quality":
                assert all(score < 4 for _, score in peptide.review_ratings)
            elif peptide.quality == "reviewed_mixed":
                assert (
                    min(s for _, s in peptide.review_ratings)
                    < 4
                    <= max(s for _, s in peptide.review_ratings)
                )
            assert list(map(list, peptide.protein_intervals)) == intervals(
                record.protein_sequence, peptide.sequence
            )


def test_queries_do_not_transfer_gene_evidence_to_another_isoform(catalog):
    tpst1 = catalog.get("c7riboseqorf56")
    assert catalog.query(gene="ENSG00000169902", transcript_id="ENST00000304842") == (tpst1,)
    assert catalog.query(transcript_id="ENST00000304842.5") == ()
    assert catalog.query(gene="TPST1", evidence="proteomics") == ()
    assert catalog.query(gene="TPST1", evidence="immunopeptidomics") == (tpst1,)
    assert catalog.query(contig="chr7", start=66210000, end=66210000, genome_build="GRCh38") == ()
    assert catalog.query(
        gene="TPST1", contig="7", start=66205483, end=66205483, genome_build="GRCh38"
    ) == (tpst1,)
    assert tpst1.map_genomic_position(66210000, genome_build="GRCh38") is None
    assert catalog.query(gene="not_a_gene") == ()
    with pytest.raises(KeyError):
        catalog.get("not_an_orf")


def test_minus_strand_spliced_codon_and_stop_mapping():
    # Transcript order: 302,301,300,105,104,103,102,101,100 -> MK*.
    record = UORFRecord(
        "minus", "gene", "GENE", "tx", "1", "-", ((100, 105), (300, 302)), "MK", "uORF", (), ()
    )
    catalog = UORFCatalog([record], {"genome_build": "GRCh38"})
    assert catalog.query(contig="chr1", start=200, end=201, genome_build="GRCh38") == ()
    for position, offset, aa, residue, role in [
        (302, 0, 1, "M", "start_codon"),
        (300, 2, 1, "M", "start_codon"),
        (105, 3, 2, "K", "body"),
        (103, 5, 2, "K", "body"),
        (100, 8, 3, "*", "stop_codon"),
    ]:
        value = record.map_genomic_position(position, genome_build="GRCh38")
        assert (
            value["coding_nucleotide_offset"],
            value["amino_acid_position"],
            value["reference_amino_acid"],
            value["role"],
        ) == (offset, aa, residue, role)
    # A codon can itself straddle two exons, on either strand.
    record = replace(record, coding_blocks=((100, 106), (300, 301)))
    assert record.map_genomic_position(106, genome_build="GRCh38")["coding_nucleotide_offset"] == 2
    assert record.map_genomic_position(106, genome_build="GRCh38")["role"] == "start_codon"
    plus = replace(record, strand="+")
    assert plus.map_genomic_position(300, genome_build="GRCh38")["coding_nucleotide_offset"] == 7
    assert plus.map_genomic_position(300, genome_build="GRCh38")["role"] == "stop_codon"


@pytest.mark.parametrize(
    "options",
    [
        {"contig": "7", "start": 1, "end": 2},
        {"genome_build": "GRCh37"},
        {"contig": "7", "start": 0, "end": 2, "genome_build": "GRCh38"},
        {"contig": "7", "start": True, "end": 2, "genome_build": "GRCh38"},
        {"contig": "7", "start": 3, "end": 2, "genome_build": "GRCh38"},
        {"evidence": "tandem"},
        {"biotype": "dORF"},
        {"require_unique": True},
        {"evidence": "ribosome", "require_reviewed": True},
        {"min_protein_length": 0},
        {"min_protein_length": True},
        {"min_protein_length": 1.5},
    ],
)
def test_ambiguous_or_invalid_queries_are_errors(catalog, options):
    with pytest.raises(ValueError):
        catalog.query(**options)


def test_position_mapping_requires_correct_build(catalog):
    native = catalog.get("c7riboseqorf56")
    with pytest.raises(TypeError):
        native.map_genomic_position(66205483)
    with pytest.raises(ValueError, match="GRCh38"):
        native.map_genomic_position(66205483, genome_build="GRCh37")


def test_reference_records_and_metadata_do_not_mutate(catalog):
    record = catalog.get("c7riboseqorf56")
    with pytest.raises(FrozenInstanceError):
        record.protein_sequence = "M"
    metadata = catalog.metadata
    metadata["sources"].clear()
    data = record.to_dict()
    data["peptides"][0]["sequence"] = "wrong"
    assert catalog.metadata["sources"] and record.peptides[0].sequence == "KIILSKISV"


def test_multiple_I_L_equivalent_intervals_and_invalid_assignments():
    peptide = UORFPeptide(
        "immunopeptidomics", "IL", ((2, 3), (4, 5)), "source", "pep", sequence_match="I_L_equivalent"
    )
    record = UORFRecord(
        "test", "gene", "GENE", "tx", "1", "+", ((10, 27),), "MLILL", "uORF", (), (peptide,)
    )
    UORFCatalog([record], {"genome_build": "GRCh38"})
    assert intervals(record.protein_sequence, peptide.sequence) == [[2, 3], [3, 4], [4, 5]]
    # Catalogue builders must preserve all matches; malformed intervals are rejected.
    with pytest.raises(ValueError, match="sequence does not match"):
        UORFCatalog(
            [replace(record, peptides=(replace(peptide, sequence="AK"),))],
            {"genome_build": "GRCh38"},
        )
    with pytest.raises(ValueError, match="coding blocks"):
        UORFCatalog([replace(record, coding_blocks=((10, 26),))], {"genome_build": "GRCh38"})
    with pytest.raises(ValueError, match="count"):
        UORFCatalog(
            [replace(record, peptides=(replace(peptide, mapping_count=True),))],
            {"genome_build": "GRCh38"},
        )


def test_packaged_dataset_manifest_and_local_replay(catalog, tmp_path):
    from orfdata import load_catalog

    assert not files("orfdata").joinpath("data/uorf-evidence").is_dir()
    root = files("orfdata").joinpath("data/orf-evidence")
    raw = root.joinpath("catalog.json.gz").read_bytes()
    manifest = json.loads(root.joinpath("manifest.json").read_text())
    full = load_catalog()
    assert hashlib.sha256(raw).hexdigest() == manifest["catalog_sha256"]
    assert manifest["statistics"]["records"] == len(full.records) == 7264
    assert len(catalog.records) == 3771
    assert len(catalog.query(evidence="ms")) == 1219
    assert len(catalog.query(evidence="ms", require_reviewed=True)) == 16
    assert catalog.records == tuple(r for r in full.records if r.biotype in ("uORF", "uoORF"))
    local = tmp_path / "copied.json.gz"
    local.write_bytes(raw)
    assert (
        load_uorf_catalog(local, expected_sha256=manifest["catalog_sha256"]).records
        == catalog.records
    )
    with pytest.raises(ValueError, match="checksum"):
        load_uorf_catalog(local, expected_sha256="0" * 64)
    invalid = json.loads(gzip.decompress(raw))
    invalid["metadata"]["coordinates"] = "0-based"
    local.write_text(json.dumps(invalid))
    with pytest.raises(ValueError, match="coordinate"):
        load_uorf_catalog(local)
    invalid = deepcopy(invalid)
    invalid["schema"] = "future_schema"
    local.write_text(json.dumps(invalid))
    with pytest.raises(ValueError, match="schema"):
        load_uorf_catalog(local)


@pytest.mark.parametrize("fmt", ["json", "tsv", "fasta"])
def test_cli_exports_native_sequence_and_peptide_evidence(fmt, capsys):
    run(["--orf-id", "c7riboseqorf56", "--format", fmt])
    output = capsys.readouterr().out
    assert "MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV" in output
    if fmt == "json":
        data = json.loads(output)
        assert data["records"][0]["peptides"][0]["protein_intervals"] == [[31, 39]]
    elif fmt == "tsv":
        assert "coding_blocks_json\t" in output and "KIILSKISV" in output
    else:
        assert output.startswith(">c7riboseqorf56 gene=TPST1")


def test_cli_position_annotation_and_review_filter(capsys):
    run(["--position", "chr7:66240380", "--genome-build", "GRCh38", "--json"])
    output = json.loads(capsys.readouterr().out)
    assert output["position_annotations"][0]["amino_acid_position"] == 32
    assert output["position_annotations"][0]["mutant_translation_evidence"] == "not_assessed"
    run(["--evidence", "immunopeptidomics", "--require-reviewed", "--json"])
    output = json.loads(capsys.readouterr().out)
    assert len(output["records"]) == 16


@pytest.mark.parametrize(
    "args, status",
    [
        (["--position", "chr7:66205483"], 2),
        (["--region", "bad"], 2),
        (["--region", "bad", "--genome-build", "GRCh38"], 2),
        (["--position", "7:1", "--format", "tsv"], 2),
        (["--orf-id", "typo"], 1),  # A well-formed ID the catalogue doesn't have.
        (["--position", "7:1", "--region", "7:1-2"], 2),
        (["--require-reviewed"], 2),
    ],
)
def test_cli_invalid_inputs_fail_explicitly(args, status, capsys):
    with pytest.raises(SystemExit) as error:
        run(args)
    assert error.value.code == status
    assert ("--help' for usage" in capsys.readouterr().err) == (status == 2)


def test_source_preparation_preserves_an_existing_Rdata_file(tmp_path, monkeypatch):
    # A user's existing unpacked source is not the builder's scratch file.
    existing = tmp_path / "filtered_peptides.Rdata"
    existing.write_bytes(b"keep this existing source")
    pins = dict(builder.PINS)
    for name in (
        "ncorf_list.xlsx",
        "peptide_sample_msrun_counts.tar.bz2",
        "author-hla-counts.json",
        "media-4_with_dot_products.csv",
        *builder.ENSEMBL_INPUTS,
        *builder.NCRNA_INPUTS,
    ):
        path = tmp_path / name
        path.write_bytes(b"cached")
        pins[name] = hashlib.sha256(b"cached").hexdigest()
    tar_path = tmp_path / "filtered_peptides.tar.bz2"
    with tarfile.open(tar_path, "w:bz2") as archive:
        member = tarfile.TarInfo("filtered_peptides.Rdata")
        member.size = 5
        archive.addfile(member, io.BytesIO(b"input"))
    pins["filtered_peptides.tar.bz2"] = hashlib.sha256(tar_path.read_bytes()).hexdigest()
    exported = b"exported mappings\n"
    pins["author-hla-mappings.tsv"] = hashlib.sha256(exported).hexdigest()
    monkeypatch.setattr(builder, "PINS", pins)

    def run_R(args, **kwargs):
        scratch = Path(args[-2])
        assert scratch != existing and scratch.read_bytes() == b"input"
        Path(args[-1]).write_bytes(exported)

    monkeypatch.setattr(builder.subprocess, "run", run_R)
    builder.prepare_sources(tmp_path)
    assert existing.read_bytes() == b"keep this existing source"
    assert (tmp_path / "author-hla-mappings.tsv").read_bytes() == exported
    assert not list(tmp_path.glob("uorf-rdata-*"))
    assert not list(tmp_path.glob(".author-hla-mappings.tsv-*"))


def test_interrupted_prepare_sources_leaves_no_truncated_files_and_retries_cleanly(tmp_path, monkeypatch):
    """Interrupted downloads and derived writes leave no file at the final path; retry needs no cleanup."""
    import re
    import zipfile

    import requests

    def tar_bz2_bytes(member_name, payload):
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:bz2") as archive:
            info = tarfile.TarInfo(member_name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
        return buf.getvalue()

    mappings_bytes = b"peptide\torf_name\tpeptide_accession\tprotein\tn_pep\nMPEPTIDE\torf1\tp1\tMPEPTIDE\t1\n"
    counts_tsv = b"MPEPTIDE\ts1\t3\tPXD0001\ttag\trun1\n"
    counts_json = (
        json.dumps(
            {"MPEPTIDE": {"datasets": ["PXD0001"], "run_count": 1, "spectrum_count": 3}},
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode()
    benchmark_bytes = b"pmid,Rating..1.5.,USI,curator,peptide,HLA,orf\n"
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(builder.BENCHMARK_MEMBER, benchmark_bytes)
    zip_bytes = zip_buf.getvalue()

    payloads = {
        "ensembl101.gtf.gz": gzip.compress(b"#!genome-build GRCh38\n", mtime=0),
        "ensembl101.cdna.fa.gz": gzip.compress(b">ENST1.1\nATG\n", mtime=0),
        "ensembl101.ncrna.fa.gz": gzip.compress(b">ENST2.1\nATG\n", mtime=0),
        "ncorf_list.xlsx": b"fake-xlsx-content-for-checksum-pin",
        "filtered_peptides.tar.bz2": tar_bz2_bytes("filtered_peptides.Rdata", b"rdata-bytes"),
        "peptide_sample_msrun_counts.tar.bz2": tar_bz2_bytes(
            "dir/peptide_sample_msrun_counts.tsv", counts_tsv
        ),
        "author-hla-mappings.tsv": mappings_bytes,
        "author-hla-counts.json": counts_json,
        "media-4_with_dot_products.csv": benchmark_bytes,
    }
    pins = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    monkeypatch.setattr(builder, "PINS", pins)

    url_map = {
        **{url: payloads[name] for name, url in {**builder.ENSEMBL_INPUTS, **builder.NCRNA_INPUTS}.items()},
        **{
            builder.AUTHOR + name: payloads[name]
            for name in (
                "ncorf_list.xlsx",
                "filtered_peptides.tar.bz2",
                "peptide_sample_msrun_counts.tar.bz2",
            )
        },
        builder.BENCHMARK_URL: zip_bytes,
    }

    fail_state = {
        "truncate_url": builder.AUTHOR + "ncorf_list.xlsx",
        "fail_rscript": True,
        "fail_zip": True,
    }

    class FakeRaw(io.BytesIO):
        def __init__(self, data, *, cut_at=None):
            super().__init__(data)
            self._cut_at = cut_at

        def read(self, amt=-1, decode_content=False):
            if self._cut_at is not None and self.tell() >= self._cut_at:
                raise requests.ConnectionError("simulated mid-transfer interruption")
            if self._cut_at is not None and (amt < 0 or self.tell() + amt > self._cut_at):
                amt = self._cut_at - self.tell()
            return super().read(amt)

    def fake_send(self, request, **kwargs):
        url = request.url
        body = url_map[url]
        response = requests.Response()
        response.url = url
        range_header = request.headers.get("Range")
        if range_header:
            if fail_state["fail_zip"] and url == builder.BENCHMARK_URL:
                raise requests.ConnectionError("simulated remote zip interruption")
            match = re.fullmatch(r"bytes=(\d+)-(\d+)", range_header)
            start, end = int(match.group(1)), int(match.group(2))
            chunk = body[start : end + 1]
            response.status_code = 206
            response.headers["Content-Range"] = f"bytes {start}-{start + len(chunk) - 1}/{len(body)}"
            response.headers["Content-Length"] = str(len(chunk))
            response.headers["ETag"] = '"benchmark-zip-v1"'
            response.raw = FakeRaw(chunk)
        else:
            response.status_code = 200
            response.headers["Content-Length"] = str(len(body))
            cut_at = len(body) // 2 if url == fail_state["truncate_url"] else None
            response.raw = FakeRaw(body, cut_at=cut_at)
        return response

    monkeypatch.setattr(requests.sessions.Session, "send", fake_send)

    def fake_rscript(args, **kwargs):
        out_path = Path(args[-1])
        if fail_state["fail_rscript"]:
            out_path.write_bytes(mappings_bytes[:10])
            raise RuntimeError("simulated Rscript interruption")
        out_path.write_bytes(mappings_bytes)

    monkeypatch.setattr(builder.subprocess, "run", fake_rscript)

    # 1. Interrupted download must not leave a truncated file under ncorf_list.xlsx.
    with pytest.raises(Exception):
        builder.prepare_sources(tmp_path)
    assert not (tmp_path / "ncorf_list.xlsx").exists()

    # 2. Retry downloads cleanly, then Rscript fails partway through writing mappings.
    fail_state["truncate_url"] = None
    with pytest.raises(RuntimeError, match="simulated Rscript interruption"):
        builder.prepare_sources(tmp_path)
    assert (tmp_path / "ncorf_list.xlsx").read_bytes() == payloads["ncorf_list.xlsx"]
    assert not (tmp_path / "author-hla-mappings.tsv").exists()

    # 3. Retry completes mappings and counts, then remote ZIP fails partway through.
    fail_state["fail_rscript"] = False
    with pytest.raises(requests.ConnectionError, match="simulated remote zip interruption"):
        builder.prepare_sources(tmp_path)
    assert (tmp_path / "author-hla-mappings.tsv").read_bytes() == mappings_bytes
    assert (tmp_path / "author-hla-counts.json").read_bytes() == counts_json
    assert not (tmp_path / "media-4_with_dot_products.csv").exists()

    # 4. Final retry completes with no manual cleanup and every file verified.
    fail_state["fail_zip"] = False
    builder.prepare_sources(tmp_path)
    for name, expected in payloads.items():
        assert (tmp_path / name).read_bytes() == expected
    assert not list(tmp_path.glob(".*"))
