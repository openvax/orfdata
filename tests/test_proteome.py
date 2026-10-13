"""Layered proteomes keep annotated, extra-ORF and viral sources separate.

The genome, catalogue and virus below are synthetic API fixtures, not observations.
"""

import io
import json

import pytest
from pyensembl import Genome

from orfdata import (
    ORFCatalog,
    ORFEvidence,
    ORFModel,
    ORFRecord,
    PeptideEvidence,
    Proteome,
    ViralProteome,
    extra_orfs,
    load_catalog,
)
from orfdata.cli import run

LNCRNA = "CCCATGAAACCCGGGTAAATGAAATAA"  # Start codon spans the exon junction.
CODING = "ATGAAACCCTAA"


def write_genome(tmp_path, reference_name="GRCm39", species="mus_musculus"):
    rows = []

    def add(contig, transcript_id, gene_id, biotype, exons, cds=None, protein_id=None):
        attrs = (f'gene_id "{gene_id}"; gene_name "{gene_id.upper()}"; transcript_id "{transcript_id}"; '
                 f'transcript_version "1"; gene_biotype "{biotype}"; transcript_biotype "{biotype}";')
        start, end = exons[0][0], exons[-1][1]
        rows.extend(f"{contig}\tfixture\t{kind}\t{start}\t{end}\t.\t+\t.\t{attrs}\n" for kind in ("gene", "transcript"))
        rows.extend(f'{contig}\tfixture\texon\t{a}\t{b}\t.\t+\t.\t{attrs} exon_number "{i + 1}"; exon_id "{transcript_id}-{i}";\n'
                    for i, (a, b) in enumerate(exons))
        if cds:
            rows.append(f'{contig}\tfixture\tCDS\t{cds[0]}\t{cds[1]}\t.\t+\t0\t{attrs} exon_number "1"; '
                        f'protein_id "{protein_id}"; protein_version "1";\n')

    add("1", "test_rna", "test_gene", "lncRNA", ((100, 104), (205, 226)))
    add("1", "coding_rna", "coding_gene", "protein_coding", ((1000, 1011),), (1000, 1008), "coding_protein")
    add("JH584304.1", "alt_rna", "alt_gene", "protein_coding", ((50, 61),), (50, 58), "alt_protein")
    gtf = tmp_path / "genome.gtf"
    gtf.write_text("".join(rows))
    cdna = tmp_path / "cdna.fa"
    cdna.write_text(f">test_rna.1\n{LNCRNA}\n>coding_rna.1\n{CODING}\n>alt_rna.1\n{CODING}\n")
    pep = tmp_path / "pep.fa"
    pep.write_text(">coding_protein.1\nMKP\n>alt_protein.1\nMKP\n")
    genome = Genome(reference_name, "synthetic", "1", gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(cdna)], protein_fasta_paths_or_urls=[str(pep)],
                    cache_directory_path=str(tmp_path))
    genome.species = species
    genome.index()
    return genome


@pytest.fixture
def genome(tmp_path):
    genome = write_genome(tmp_path)
    yield genome
    genome.close()


def orf(orf_id, transcript_id="test_rna", blocks=((103, 104), (205, 217)), protein="MKPG", **evidence):
    model = ORFModel(orf_id, transcript_id, "test_gene", "TEST", "1", "+", blocks, protein, "lncRNA",
                     "GRCm39", 10090, "mus_musculus", source="synthetic")
    return ORFRecord.from_model(model, ORFEvidence(**evidence))


@pytest.fixture
def catalog():
    peptide = (PeptideEvidence("shotgun", "KPG", ((2, 4),), "synthetic", "fixture"),)
    records = [
        orf("ms-orf", peptides=peptide),
        orf("ribo-orf", blocks=((218, 226),), protein="MK", ribosome_studies=("synthetic-study",)),
        orf("absent-orf", transcript_id="missing_rna", peptides=peptide),
        orf("changed-orf", protein="MKPA", peptides=(PeptideEvidence("shotgun", "MKP", ((1, 3),), "synthetic", "fixture"),)),
    ]
    return ORFCatalog(records, dict(genome_build="GRCm39", taxon_id=10090, species="mus_musculus",
                                    coordinates="1-based-inclusive", dataset_version="fixture-v1",
                                    sources=[dict(id="synthetic", scope="test fixture")]))


def virus(hosts=(10090,)):
    return ViralProteome("synthetic-virus", 999999, hosts, (("virus-protein", "MSTVK"),),
                         "synthetic fixture", "test-declared host range")


def test_annotated_proteome_matches_the_genome_annotation(genome):
    proteome = Proteome.from_genome(genome)
    assert [(e.accession, e.sequence, e.kind, e.contig) for e in proteome] == [
        ("alt_protein.1", "MKP", "annotated", "JH584304.1"), ("coding_protein.1", "MKP", "annotated", "1")]
    entry = proteome.entries[1]
    assert (entry.transcript_id, entry.gene_id, entry.biotype, entry.taxon_id) == ("coding_rna.1", "coding_gene", "protein_coding", 10090)
    assert proteome.summary()["distinct_sequences"] == 1
    assert [dict(layer)["kind"] for layer in proteome.layers] == ["annotated"]
    primary = proteome.filter(lambda e: "." not in e.contig, description="primary contigs")
    assert [e.accession for e in primary] == ["coding_protein.1"]
    assert dict(primary.layers[-1]) == dict(kind="filter", description="primary contigs", removed=1, added=0)


def test_extra_orfs_are_selected_by_evidence_and_unbound_ones_are_reported(genome, catalog):
    base = Proteome.from_genome(genome)
    proteome = base.with_extra_orfs(catalog=catalog)  # MS support by default
    extra = [e for e in proteome if e.kind == "extra_orf"]
    assert [(e.accession, e.sequence, e.evidence, e.biotype) for e in extra] == [("ms-orf", "MKPG", ("shotgun",), "lncRNA")]
    assert extra[0].orf.transcript.biotype == "lncRNA" and extra[0].transcript_id == "test_rna.1"
    reasons = {e.orf_id: e.reason for e in proteome.excluded}
    assert reasons == {"absent-orf": "Transcript absent from this annotation", "changed-orf": "ORF reference protein/stop mismatch"}
    layer = dict(proteome.layers[-1])
    assert (layer["evidence"], layer["selected"], layer["added"], layer["excluded"]) == ("ms", 3, 1, 2)
    assert layer["source"] == "orfdata-catalog:fixture-v1"
    everything = base.with_extra_orfs(None, catalog=catalog)
    assert {e.accession for e in everything if e.kind == "extra_orf"} == {"ms-orf", "ribo-orf"}
    ribosome = base.with_extra_orfs("ribosome", catalog=catalog)
    assert [e.accession for e in ribosome if e.kind == "extra_orf"] == ["ribo-orf"]
    assert base.with_extra_orfs("ms", catalog=catalog, require_resolved=True).layers[-1]["added"] == 1
    with pytest.raises(ValueError, match="Duplicate proteome accessions"):
        proteome.with_extra_orfs(catalog=catalog)


def test_extra_orfs_lookup_by_transcript_and_gene(genome, catalog):
    found = extra_orfs(genome, catalog=catalog)
    assert [o.model.orf_id for o in found] == ["ms-orf", "ribo-orf"]
    transcript = genome.transcript_by_id("test_rna")
    assert found.for_transcript(transcript) == found.for_transcript("test_rna.1") == found.orfs
    assert found.for_gene("TEST") == found.orfs and found.for_gene("coding_gene") == ()
    assert {e.orf_id for e in found.excluded} == {"absent-orf", "changed-orf"}


def test_uninstalled_annotation_is_never_downloaded_or_indexed(tmp_path):
    genome = write_genome(tmp_path)
    genome.close()
    unindexed = Genome("GRCm39", "synthetic", "1", gtf_path_or_url=genome.gtf_path,
                       transcript_fasta_paths_or_urls=genome.transcript_fasta_paths,
                       cache_directory_path=str(tmp_path / "empty-cache"))
    unindexed.species = "mus_musculus"
    for build in (Proteome.from_genome, extra_orfs):
        with pytest.raises(ValueError, match="annotation is not indexed; install it first"):
            build(unindexed)
    assert not (tmp_path / "empty-cache").exists()


def test_species_without_a_catalogue_or_mismatched_catalogue_fails(genome):
    base = Proteome.from_genome(genome)
    with pytest.raises(ValueError, match="No extra-ORF catalogue is available for mus_musculus GRCm39"):
        base.with_extra_orfs()
    with pytest.raises(ValueError, match="cannot bind"):
        base.with_extra_orfs(catalog=load_catalog())


def test_bundled_human_catalogue_is_the_default_for_grch38(tmp_path):
    human = write_genome(tmp_path, "GRCh38", "homo_sapiens")
    try:
        proteome = Proteome.from_genome(human).with_extra_orfs()
    finally:
        human.close()
    layer = dict(proteome.layers[-1])
    # No bundled ORF's transcript is in this slice: all are accounted for, none dropped silently.
    assert (layer["selected"], layer["added"], layer["excluded"]) == (2186, 0, 2186)
    assert layer["source"] == "orfdata-catalog:" + load_catalog().metadata["dataset_version"]


def test_viral_proteins_require_a_recorded_host(genome):
    base = Proteome.from_genome(genome)
    with pytest.raises(ValueError, match="no recorded host matching mus_musculus"):
        base.with_viral(virus(hosts=(9606,)))
    proteome = base.with_viral(virus())
    entry = proteome.entries[-1]
    assert (entry.accession, entry.kind, entry.taxon_id, entry.source) == ("virus-protein", "viral", 999999, "viral:synthetic-virus")
    layer = dict(proteome.layers[-1])
    assert layer["host_taxon_ids"] == [10090] and layer["host_source"] == "test-declared host range"
    for bad in (dict(hosts=()), dict(hosts=("mouse",))):
        with pytest.raises(ValueError):
            virus(**bad)


def test_viral_fasta_is_pinned_and_parsed(tmp_path):
    path = tmp_path / "virus.fa"
    path.write_text(">p1 description\nMSTV\nK*\n>p2\nmaa\n")
    parsed = ViralProteome.from_fasta(path, name="v", taxon_id=1, host_taxon_ids={9606}, host_source="test")
    assert parsed.proteins == (("p1", "MSTVK"), ("p2", "MAA")) and len(parsed.source_sha256) == 64


def test_fasta_export_round_trips_sequences_and_provenance(genome, catalog):
    proteome = Proteome.from_genome(genome).with_extra_orfs(catalog=catalog).with_viral(virus())
    handle = io.StringIO()
    proteome.write_fasta(handle)
    records, header = {}, None
    for line in handle.getvalue().splitlines():
        if line.startswith(">"):
            header = line
            records[header] = ""
        else:
            records[header] += line
    assert records == {e.fasta_header(): e.sequence for e in proteome}
    assert ">ms-orf kind=extra_orf taxon=10090 source=orfdata-catalog:fixture-v1 gene=test_gene gene_name=TEST transcript=test_rna.1 contig=1 biotype=lncRNA evidence=shotgun" in records


def test_cli_writes_fasta_and_reports_layers(genome, catalog, tmp_path, capsys, monkeypatch):
    from orfdata import proteome_cli

    monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda release, species: genome)
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog.to_dict()))
    output = tmp_path / "proteome.fa"
    run(["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path), "--output", str(output)])
    summary = json.loads(capsys.readouterr().err)
    assert summary["kinds"] == {"annotated": 2, "extra_orf": 1}
    assert {e["orf_id"] for e in summary["excluded_orfs"]} == {"absent-orf", "changed-orf"}
    assert output.read_text().count(">") == 3
    run(["proteome", "--species", "mouse", "--release", "1", "--extra-orf-evidence", "none"])
    captured = capsys.readouterr()
    assert captured.out.count(">") == 2 and json.loads(captured.err)["kinds"] == {"annotated": 2}
    with pytest.raises(SystemExit):
        run(["proteome", "--species", "mouse", "--release", "1"])  # No mouse catalogue yet.
    assert "No extra-ORF catalogue" in capsys.readouterr().err
