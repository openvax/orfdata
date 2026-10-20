"""Layered proteomes keep annotated, extra-ORF and viral sources separate.

The genome, catalogue and virus below are synthetic API fixtures, not observations.
"""

import gzip
import io
import json
from dataclasses import replace
from types import SimpleNamespace

import pytest
from pyensembl import Genome

from orfdata import (
    MS,
    ExtraORFs,
    InitiationContext,
    ORFCatalog,
    ORFEvidence,
    ORFModel,
    ORFRecord,
    PeptideEvidence,
    Proteome,
    ReferenceContext,
    ReferenceORF,
    ReferenceTranscriptome,
    RiboSeq,
    TranscriptORF,
    ViralProteome,
    annotated_orfs,
    extra_orfs,
    load_catalog,
    viral_reference,
)
from orfdata.cli import run
from orfdata.proteome import _sha256

LNCRNA = "CCCATGAAACCCGGGTAAATGAAATAA"  # Start codon spans the exon junction.
CODING = "ATGAAACCCTAA"
MIXED = "ATGAAATAACCCATGCCCTAA"  # uORF MK, then the annotated CDS MP.
ENS_RNA, ENS_GENE = "ENSMUST00000000001", "ENSMUSG00000000001"  # Synthetic Ensembl-style IDs.


def write_genome(path, reference_name="GRCm39", species="mus_musculus", *, lnc_id="test_rna",
                 lnc_gene="test_gene", version="1", twin=False, fragment=False):
    path.mkdir(parents=True, exist_ok=True)
    rows, genes = [], {}

    def add(contig, transcript_id, gene_id, biotype, exons, cds=None, protein_id=None, codons=False):
        attrs = (f'gene_id "{gene_id}"; gene_version "{version}"; gene_name "{gene_id.upper()}"; '
                 f'transcript_id "{transcript_id}"; transcript_version "{version}"; '
                 f'gene_biotype "{biotype}"; transcript_biotype "{biotype}";')
        start, end = exons[0][0], exons[-1][1]
        gene_attrs = (f'gene_id "{gene_id}"; gene_version "{version}"; gene_name "{gene_id.upper()}"; '
                      f'gene_biotype "{biotype}";')
        span = genes.get(gene_id, (contig, start, end, gene_attrs))
        genes[gene_id] = (contig, min(start, span[1]), max(end, span[2]), gene_attrs)
        rows.append(f"{contig}\tfixture\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}\n")
        rows.extend(f'{contig}\tfixture\texon\t{a}\t{b}\t.\t+\t.\t{attrs} exon_number "{i + 1}"; exon_id "{transcript_id}-{i}";\n'
                    for i, (a, b) in enumerate(exons))
        if cds:
            protein = f' exon_number "1"; protein_id "{protein_id}"; protein_version "1";'
            rows.append(f'{contig}\tfixture\tCDS\t{cds[0]}\t{cds[1]}\t.\t+\t0\t{attrs}{protein}\n')
            if codons:  # Ensembl marks the stop codon separately from the CDS.
                rows.append(f'{contig}\tfixture\tstart_codon\t{cds[0]}\t{cds[0] + 2}\t.\t+\t0\t{attrs}\n')
                rows.append(f'{contig}\tfixture\tstop_codon\t{cds[1] + 1}\t{cds[1] + 3}\t.\t+\t0\t{attrs}\n')

    add("1", lnc_id, lnc_gene, "lncRNA", ((100, 104), (205, 226)))
    # A longer isoform of the same gene carrying the identical spliced ORF paths.
    add("1", lnc_id + "_long", lnc_gene, "lncRNA", ((95, 104), (205, 226)))
    add("1", "mixed_rna", "mixed_gene", "protein_coding", ((2000, 2020),), (2012, 2017), "mixed_protein", True)
    # Another gene's RNA whose sequence contains mixed_rna's CDS, as MIEF1's mRNA contains MIURF's.
    add("1", "host_rna", "host_gene", "lncRNA", ((1995, 2025),))
    add("1", "coding_rna", "coding_gene", "protein_coding", ((1000, 1011),), (1000, 1008), "coding_protein")
    if twin:  # A second transcript and protein ID with the identical CDS.
        add("1", "twin_rna", "coding_gene", "protein_coding", ((1000, 1011),), (1000, 1008), "twin_protein")
    if fragment:  # A 3'-truncated RNA whose CDS has no stop: coding_rna reads on past it.
        add("1", "fragment_rna", "coding_gene", "protein_coding", ((1000, 1005),), (1000, 1005), "fragment_protein")
    add("JH584304.1", "alt_rna", "alt_gene", "protein_coding", ((50, 61),), (50, 58), "alt_protein")
    rows += [f"{contig}\tfixture\tgene\t{a}\t{b}\t.\t+\t.\t{attrs}\n" for contig, a, b, attrs in genes.values()]
    gtf = path / "genome.gtf"
    gtf.write_text("".join(rows))
    cdna = path / "cdna.fa"
    cdna.write_text(f">{lnc_id}.{version}\n{LNCRNA}\n>{lnc_id}_long.{version}\nGGGGG{LNCRNA}\n"
                    f">coding_rna.{version}\n{CODING}\n>alt_rna.{version}\n{CODING}\n>mixed_rna.{version}\n{MIXED}\n"
                    f">host_rna.{version}\nCCCCC{MIXED}GGGGG\n")
    with cdna.open("a") as handle:
        if twin:
            handle.write(f">twin_rna.{version}\n{CODING}\n")
        if fragment:
            handle.write(f">fragment_rna.{version}\n{CODING[:6]}\n")
    pep = path / "pep.fa"
    pep.write_text(">coding_protein.1\nMKP\n>alt_protein.1\nMKP\n>mixed_protein.1\nMP\n>twin_protein.1\nMKP\n"
                   ">fragment_protein.1\nMK\n")
    genome = Genome(reference_name, "synthetic", "1", gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(cdna)], protein_fasta_paths_or_urls=[str(pep)],
                    cache_directory_path=str(path))
    genome.species = species
    genome.index()
    return genome


@pytest.fixture
def genome(tmp_path):
    genome = write_genome(tmp_path / "genome")
    yield genome
    genome.close()


def orf(orf_id, transcript_id="test_rna", blocks=((103, 104), (205, 217)), protein="MKPG", **evidence):
    model = ORFModel(orf_id, transcript_id, "test_gene", "TEST", "1", "+", blocks, protein, "lncRNA",
                     "GRCm39", 10090, "mus_musculus", source="synthetic")
    return ORFRecord.from_model(model, ORFEvidence(**evidence))


def make_catalog(transcript_id="test_rna"):
    peptide = (PeptideEvidence("proteomics", "KPG", ((2, 4),), "synthetic", "fixture"),)
    # Spliced start flanks: CCC | ATG (across the junction) | AAACCC.
    flanks = InitiationContext("CCC", "ATG", "AAACCC", (100, 101, 102, 103, 104, 205, 206, 207, 208, 209, 210, 211))
    records = [
        replace(orf("ms-orf", transcript_id, peptides=peptide), initiation_context=flanks),
        orf("ribo-orf", transcript_id, blocks=((218, 226),), protein="MK", ribosome_studies=("synthetic-study",)),
        orf("absent-orf", transcript_id="missing_rna", peptides=peptide),
        orf("changed-orf", transcript_id, protein="MKPA",
            peptides=(PeptideEvidence("proteomics", "MKP", ((1, 3),), "synthetic", "fixture"),)),
    ]
    return ORFCatalog(records, dict(genome_build="GRCm39", taxon_id=10090, species="mus_musculus",
                                    coordinates="1-based-inclusive", dataset_version="fixture-v1",
                                    sources=[dict(id="synthetic", scope="test fixture")]))


@pytest.fixture
def catalog():
    return make_catalog()


@pytest.fixture
def ensembl_style(tmp_path):
    """Build genomes that share Ensembl-style IDs but differ in version or reference."""
    genomes = []

    def make(name, reference_name="GRCm39", species="mus_musculus", version="6"):
        genomes.append(write_genome(tmp_path / name, reference_name, species, lnc_id=ENS_RNA,
                                    lnc_gene=ENS_GENE, version=version))
        return genomes[-1]

    yield make
    for genome in genomes:
        genome.close()


def virus(hosts=(10090,)):
    return ViralProteome("synthetic-virus", 999999, hosts, (("virus-protein", "MSTVK"),),
                         "synthetic fixture", "test-declared host range")


def test_annotated_proteome_matches_the_genome_annotation(genome):
    proteome = Proteome.from_genome(genome)
    assert [(e.accession, e.sequence, e.kind, e.contig) for e in proteome] == [
        ("alt_protein.1", "MKP", "annotated", "JH584304.1"), ("coding_protein.1", "MKP", "annotated", "1"),
        ("mixed_protein.1", "MP", "annotated", "1")]
    coding, mixed = proteome.entries[1], proteome.entries[2]
    assert (coding.transcript_id, coding.gene_id, coding.transcript_biotype, coding.orf_type, coding.taxon_id) == (
        "coding_rna.1", "coding_gene", "protein_coding", None, 10090)
    # Coordinates slice the transcript sequence; this CDS has no stop codon feature.
    occurrence, = mixed.occurrences
    assert (occurrence.start, occurrence.end, occurrence.coding_blocks) == (12, 21, ((2012, 2020),))
    assert MIXED[occurrence.start:occurrence.end] == "ATGCCCTAA" and occurrence.on_source_transcript
    assert (coding.occurrences[0].start, coding.occurrences[0].end) == (0, 9)
    # Identical proteins from different genes stay separate entries, grouped on request.
    assert [e.gene_id for e in proteome.by_sequence()["MKP"]] == ["alt_gene", "coding_gene"]
    assert proteome.summary()["distinct_sequences"] == 2
    assert [layer["role"] for layer in proteome.layers] == ["host"]
    primary = proteome.filter(lambda e: "." not in e.contig, description="primary contigs")
    assert [e.accession for e in primary] == ["coding_protein.1", "mixed_protein.1"]
    assert dict(primary.layers[-1]) == dict(kind="filter", role="filter", description="primary contigs",
                                            removed=1, added=0)


def test_extra_orfs_sit_on_every_transcript_carrying_their_coding_path(genome, catalog):
    found = extra_orfs(genome, catalog=catalog)
    assert [(o.orf_id, o.transcript.id, o.on_source_transcript) for o in found] == [
        ("ms-orf", "test_rna", True), ("ms-orf", "test_rna_long", False),
        ("ribo-orf", "test_rna", True), ("ribo-orf", "test_rna_long", False),
        ("absent-orf", "test_rna", False), ("absent-orf", "test_rna_long", False)]
    long_form = found.orfs[1]
    assert (long_form.start, long_form.end) == (8, 23)
    assert long_form.transcript.sequence[long_form.start:long_form.end] == long_form.coding_sequence
    assert long_form.protein_sequence == "MKPG" and long_form.source_transcript_id == "test_rna"
    # Source-binding history is kept even when another transcript carries the ORF.
    assert {f.orf_id: f.reason for f in found.source_failures} == {
        "absent-orf": "Transcript absent from this annotation", "changed-orf": "ORF reference protein/stop mismatch"}
    assert [f.orf_id for f in found.excluded] == ["changed-orf"]
    assert [o.transcript.id for o in found.for_orf("absent-orf")] == ["test_rna", "test_rna_long"]
    assert [o.on_source_transcript for o in found.for_orf("ms-orf")] == [True, False]
    assert found.for_orf("changed-orf") == ()
    with pytest.raises(ValueError, match="identity mismatch"):
        catalog.get("ms-orf").bind(genome.transcript_by_id("test_rna_long"))  # bind() stays strict.
    bound = found.orfs[0]
    assert (bound.orf_id, bound.id, bound.orf_type, bound.transcript_id) == ("ms-orf", "ms-orf", "lncRNA", "test_rna")
    # Gene fields describe this genome; the source's own labels are kept alongside.
    assert (bound.gene_id, bound.gene_name, bound.source_gene_name) == ("test_gene", "TEST_GENE", "TEST")
    assert bound.gene.id == "test_gene" and bound.coding_blocks == ((103, 104), (205, 217))
    exported = json.loads(json.dumps(bound.to_dict()))
    assert (exported["start"], exported["end"], exported["on_source_transcript"]) == (3, 18, True)


def test_proteome_entries_keep_every_occurrence_and_source_binding(genome, catalog):
    base = Proteome.from_genome(genome)
    proteome = base.with_extra_orfs(catalog=catalog)  # MS support by default
    extra = {e.accession: e for e in proteome if e.kind == "extra_orf"}
    assert set(extra) == {"ms-orf", "absent-orf"}
    ms = extra["ms-orf"]
    # The host transcript's biotype and the source ORF category are separate fields.
    assert (ms.sequence, ms.evidence, ms.transcript_biotype, ms.orf_type, ms.source_binding) == (
        "MKPG", ("proteomics",), "lncRNA", "lncRNA", None)
    assert [(o.transcript_id, o.on_source_transcript, o.start, o.end) for o in ms.occurrences] == [
        ("test_rna.1", True, 3, 18), ("test_rna_long.1", False, 8, 23)]
    assert ms.evidence_detail["peptides"][0]["sequence"] == "KPG"
    rescued = extra["absent-orf"]
    assert rescued.source_binding == "Transcript absent from this annotation"
    assert [o.on_source_transcript for o in rescued.occurrences] == [False, False]
    assert [f.orf_id for f in proteome.excluded] == ["changed-orf"]
    layer = dict(proteome.layers[-1])
    assert (layer["role"], layer["evidence"], layer["selected"], layer["added"], layer["excluded"],
            layer["not_on_source_transcript"]) == ("host", MS().to_dict(), 3, 2, 1, 1)
    assert layer["source"] == "orfdata-catalog:fixture-v1" and len(layer["catalog"]["sha256"]) == 64
    everything = base.with_extra_orfs(None, catalog=catalog)
    assert {e.accession for e in everything if e.kind == "extra_orf"} == {"ms-orf", "ribo-orf", "absent-orf"}
    ribosome = base.with_extra_orfs("ribosome", catalog=catalog)
    assert [e.accession for e in ribosome if e.kind == "extra_orf"] == ["ribo-orf"]
    assert base.with_extra_orfs("ms", catalog=catalog, require_resolved=True).layers[-1]["added"] == 2
    with pytest.raises(ValueError, match="Duplicate proteome accessions"):
        proteome.with_extra_orfs(catalog=catalog)


def test_extra_orfs_lookup_by_transcript_and_gene(genome, catalog):
    found = extra_orfs(genome, catalog=catalog)
    transcript = genome.transcript_by_id("test_rna")
    on_test_rna = found.for_transcript(transcript)
    assert [o.orf_id for o in on_test_rna] == ["absent-orf", "ms-orf", "ribo-orf"]  # RNA order.
    assert found.for_transcript("test_rna.1") == on_test_rna
    assert found.for_transcript("test_rna.2") == ()  # Not an Ensembl version: no stripping.
    # The genome names this gene TEST_GENE; the source catalogue calls it TEST.
    assert found.for_gene("TEST_GENE") == found.for_gene("TEST", source=True) == found.orfs
    assert found.for_gene("TEST") == found.for_gene("coding_gene") == found.for_gene("") == ()


def mixed_catalog():
    uorf = ORFModel("mixed-uorf", "mixed_rna", "mixed_gene", "MIXED", "1", "+", ((2000, 2008),), "MK", "uORF",
                    "GRCm39", 10090, "mus_musculus", source="synthetic")
    record = ORFRecord.from_model(uorf, ORFEvidence(ribosome_studies=("synthetic-study",)))
    return ORFCatalog([record], dict(genome_build="GRCm39", taxon_id=10090, species="mus_musculus",
                                     coordinates="1-based-inclusive", dataset_version="mixed-v1",
                                     sources=[dict(id="synthetic", scope="test fixture")]))


def test_annotated_cds_is_a_transcript_orf_taken_as_annotated(genome):
    assert ReferenceORF is TranscriptORF  # The original name still works for Varcode and Isovar.
    complete = TranscriptORF.from_annotation(genome.transcript_by_id("mixed_rna"))
    assert (complete.kind, complete.orf_id, complete.id, complete.on_source_transcript) == (
        "annotated", "mixed_protein.1", "mixed_protein.1", True)
    assert (complete.coding_sequence, complete.protein_sequence, complete.start, complete.end) == (
        "ATGCCCTAA", "MP", 12, 21)
    assert complete.start_complete and complete.stop_complete and complete.coding_blocks == ((2012, 2020),)
    assert complete.model is None and complete.orf_type is None and complete.genetic_code is None
    assert not complete.evidence.has_evidence("ribosome") and complete.gene_id == "mixed_gene"
    # No start or stop codon features: the CDS is reported as incomplete, not completed.
    partial = TranscriptORF.from_annotation(genome.transcript_by_id("coding_rna"))
    assert (partial.coding_sequence, partial.start_complete, partial.stop_complete) == ("ATGAAACCC", False, False)
    exported = json.loads(json.dumps(complete.to_dict()))
    assert (exported["kind"], exported["orf_id"], exported["model"], exported["stop_complete"]) == (
        "annotated", "mixed_protein.1", None, True)
    with pytest.raises(ValueError, match="no annotated CDS"):
        TranscriptORF.from_annotation(genome.transcript_by_id("test_rna"))


def test_another_genes_annotated_cds_on_a_carrying_transcript(genome):
    # Like MIURF's CDS in MIEF1's mRNA: same coding path, a different gene's RNA.
    host = genome.transcript_by_id("host_rna")
    carried = TranscriptORF.from_annotation(genome.transcript_by_id("mixed_rna"), on=host)
    assert (carried.orf_id, carried.transcript_id, carried.start, carried.end) == ("mixed_protein.1", "host_rna", 17, 26)
    assert (carried.gene_id, carried.source_gene_id, carried.source_transcript_id) == ("host_gene", "mixed_gene", "mixed_rna.1")
    assert not carried.on_source_transcript and host.sequence[17:26] == carried.coding_sequence
    assert annotated_orfs(host) == (carried,) and carried.equivalent_ids == ()
    with pytest.raises(ValueError, match="not one contiguous part"):
        TranscriptORF.from_annotation(genome.transcript_by_id("mixed_rna"), on=genome.transcript_by_id("test_rna"))


def test_identical_annotated_placements_are_reported_once_with_their_ids(tmp_path):
    genome = write_genome(tmp_path / "twins", twin=True, fragment=True)
    try:
        coding = genome.transcript_by_id("coding_rna")
        # fragment_rna's stopless CDS is not an ORF on coding_rna, which reads on past it.
        placement, = annotated_orfs(coding)
        assert (placement.orf_id, placement.on_source_transcript, placement.equivalent_ids) == (
            "coding_protein.1", True, ("twin_protein.1",))
        twin, = annotated_orfs(genome.transcript_by_id("twin_rna"))
        assert (twin.orf_id, twin.equivalent_ids) == ("twin_protein.1", ("coding_protein.1",))
        # The proteome keeps both proteins: they are distinct annotated identities.
        assert {"coding_protein.1", "twin_protein.1"} <= {e.accession for e in Proteome.from_genome(genome)}
    finally:
        genome.close()


def test_dotted_non_ensembl_ids_are_compared_whole(tmp_path):
    # C. elegans transcript IDs contain dots; Y74C9A.3.1 is not version 1 of Y74C9A.3.
    rows, rna = [], "CCCATGAAACCCGGGTAAGG"
    for transcript_id in ("Y74C9A.3.1", "Y74C9A.3.2"):
        attrs = (f'gene_id "WBGene1"; gene_name "gen-1"; transcript_id "{transcript_id}"; '
                 'gene_biotype "lncRNA"; transcript_biotype "lncRNA";')
        rows += [f"I\tfixture\ttranscript\t100\t119\t.\t+\t.\t{attrs}\n",
                 f'I\tfixture\texon\t100\t119\t.\t+\t.\t{attrs} exon_number "1"; exon_id "{transcript_id}-e";\n']
    rows.append('I\tfixture\tgene\t100\t119\t.\t+\t.\tgene_id "WBGene1"; gene_name "gen-1"; gene_biotype "lncRNA";\n')
    (tmp_path / "g.gtf").write_text("".join(rows))
    (tmp_path / "c.fa").write_text(f">Y74C9A.3.1\n{rna}\n>Y74C9A.3.2\n{rna}\n")
    worm = Genome("WBcel235", "synthetic", "1", gtf_path_or_url=str(tmp_path / "g.gtf"),
                  transcript_fasta_paths_or_urls=[str(tmp_path / "c.fa")], cache_directory_path=str(tmp_path))
    worm.species = "caenorhabditis_elegans"
    worm.index()
    try:
        model = ORFModel("worm-orf", "Y74C9A.3.1", "WBGene1", "gen-1", "I", "+", ((103, 117),), "MKPG", "lncRNA",
                         "WBcel235", 6239, "caenorhabditis_elegans", source="synthetic")
        catalog = ORFCatalog([ORFRecord.from_model(model, ORFEvidence())], dict(
            genome_build="WBcel235", taxon_id=6239, species="caenorhabditis_elegans",
            coordinates="1-based-inclusive", dataset_version="fixture", sources=[dict(id="s", scope="test")]))
        context = ReferenceContext(6239, "caenorhabditis_elegans", "WBcel235", "synthetic", "1")
        found = extra_orfs(worm, catalog=catalog, context=context)
        assert [(o.transcript.id, o.on_source_transcript) for o in found] == [
            ("Y74C9A.3.1", True), ("Y74C9A.3.2", False)]
        assert found.source_failures == () and found.excluded == ()
    finally:
        worm.close()


def test_one_transcript_lists_its_annotated_cds_with_extra_orfs_in_rna_order(genome):
    found = extra_orfs(genome, catalog=mixed_catalog())
    transcript = genome.transcript_by_id("mixed_rna")
    assert [o.orf_id for o in found.for_transcript(transcript)] == ["mixed-uorf"]
    every = found.for_transcript(transcript, include_annotated=True)
    assert [(o.kind, o.orf_id, o.start, o.end) for o in every] == [
        ("extra_orf", "mixed-uorf", 0, 9), ("annotated", "mixed_protein.1", 12, 21)]
    assert found.for_transcript("mixed_rna.1", include_annotated=True) == every
    # The other gene's RNA carries both the uORF and the annotated CDS.
    on_host = found.for_transcript("host_rna", include_annotated=True)
    assert [(o.kind, o.orf_id, o.start, o.on_source_transcript) for o in on_host] == [
        ("extra_orf", "mixed-uorf", 5, False), ("annotated", "mixed_protein.1", 17, False)]
    assert found.for_transcript("coding_rna", include_annotated=True)[0].orf_id == "coding_protein.1"
    assert found.for_transcript("test_rna", include_annotated=True) == ()  # Noncoding, no catalogue ORFs.
    detached = ExtraORFs(found.orfs, (), found.catalog, found.reference)
    with pytest.raises(ValueError, match="Pass a transcript object"):
        detached.for_transcript("coding_rna", include_annotated=True)


def test_overlapping_uses_coding_paths_and_optional_initiation_flanks(genome, catalog):
    found = extra_orfs(genome, catalog=catalog)

    def at(start, end, **kw):
        return sorted({o.orf_id for o in found.overlapping("chr1", start, end, genome_build="GRCm39", **kw)})

    assert at(105, 204) == []  # Intron between the ORF's coding blocks.
    assert at(100, 103) == ["absent-orf", "ms-orf"] and at(218, 218) == ["ribo-orf"]
    # absent-orf has no source context; on test_rna its flanks are read from that RNA.
    assert at(101, 101) == [] and at(101, 101, include_initiation_context=True) == ["absent-orf", "ms-orf"]
    for orf in (o for o in found if o.orf_id == "ms-orf"):
        context, rna = orf.initiation_context, orf.transcript_sequence
        assert context.upstream == rna[orf.start - len(context.upstream) : orf.start]
        assert context.source_transcript_id == (orf.record.initiation_context.source_transcript_id
                                                if orf.on_source_transcript else orf.transcript.versioned_id)
    assert found.overlapping("2", 1, 10**6, genome_build="GRCm39") == ()
    with pytest.raises(ValueError, match="Expected genome build GRCm39"):
        found.overlapping("1", 100, 103, genome_build="GRCm38")
    with pytest.raises(ValueError, match="end >= start"):
        found.overlapping("1", 103, 100, genome_build="GRCm39")


def test_snapshot_round_trips_offline_and_rejects_corruption(genome, catalog, tmp_path):
    proteome = Proteome.from_genome(genome).with_extra_orfs(catalog=catalog).with_viral(virus())
    path = tmp_path / "proteome.json.gz"
    identity = proteome.save(path)
    loaded = Proteome.load(path)
    assert identity == proteome.identity == loaded.identity and loaded.genome is None
    assert [e.to_dict() for e in loaded] == [e.to_dict() for e in proteome]
    assert loaded.excluded == proteome.excluded and loaded.layers == proteome.layers
    assert [e.kind for e in loaded.host()] == ["annotated"] * 3 + ["extra_orf"] * 2
    # Identity is content-based, and saving is deterministic: re-saving gives the same bytes.
    assert loaded.save(tmp_path / "again.json.gz") == identity
    assert (tmp_path / "again.json.gz").read_bytes() == path.read_bytes()
    with pytest.raises(ValueError, match="bind to a genome"):
        loaded.with_extra_orfs(catalog=catalog)
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    snapshot["content"]["entries"][0]["sequence"] = "MKA"
    corrupt = tmp_path / "corrupt.json"
    corrupt.write_text(json.dumps(snapshot))
    with pytest.raises(ValueError, match="does not match its identity"):
        Proteome.load(corrupt)
    # Hash-consistent but non-canonical content does not load with a mismatched identity.
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    snapshot["content"]["entries"][0]["evidence"] = "ab"
    snapshot["identity"] = _sha256(snapshot["content"])
    corrupt.write_text(json.dumps(snapshot))
    with pytest.raises(ValueError, match="does not round-trip"):
        Proteome.load(corrupt)


def test_snapshot_identity_follows_content_selection_and_catalogue(genome, catalog):
    base = Proteome.from_genome(genome)
    ms = base.with_extra_orfs(catalog=catalog)
    assert ms.identity == base.with_extra_orfs(catalog=make_catalog()).identity
    assert ms.identity != base.with_extra_orfs("ribosome", catalog=catalog).identity
    assert ms.identity != ms.filter(lambda e: e.kind == "annotated").identity
    # The same selected ORFs from a different catalogue release are a different reference.
    changed = make_catalog()
    changed = ORFCatalog(changed.records, dict(changed.metadata, dataset_version="fixture-v2"))
    assert ms.identity != base.with_extra_orfs(catalog=changed).identity


def test_snapshot_rejects_wrong_species_and_reference_declarations(genome, catalog, tmp_path):
    proteome = Proteome.from_genome(genome).with_extra_orfs(catalog=catalog).with_viral(virus())
    content = proteome._content()

    def load(content):
        path = tmp_path / "edited.json"
        path.write_text(json.dumps(dict(schema="orfdata.proteome.v1", identity=_sha256(content), content=content)))
        return Proteome.load(path)

    human = json.loads(json.dumps(content))
    human["entries"][3]["taxon_id"] = 9606
    with pytest.raises(ValueError, match="from another species"):
        load(human)
    orphan = json.loads(json.dumps(content))
    orphan["layers"] = [layer for layer in orphan["layers"] if layer["role"] != "viral"]
    with pytest.raises(ValueError, match="no host layer"):
        load(orphan)
    elsewhere = json.loads(json.dumps(content))
    elsewhere["layers"][0]["reference"] = "ref:" + "0" * 64
    with pytest.raises(ValueError, match="different reference"):
        load(elsewhere)


def test_transcript_lookup_keeps_version_and_reference_identity(ensembl_style):
    # Regression for #62: lookups once stripped versions and ignored the query's genome.
    genome = ensembl_style("main")
    found = extra_orfs(genome, catalog=make_catalog(ENS_RNA))
    on_rna = found.for_transcript(ENS_RNA)
    assert [o.orf_id for o in on_rna] == ["absent-orf", "ms-orf", "ribo-orf"]
    assert found.for_transcript(ENS_RNA + ".6") == found.for_transcript(genome.transcript_by_id(ENS_RNA)) == on_rna
    with pytest.raises(ValueError, match=r"ENSMUST00000000001\.999 requested; this reference has ENSMUST00000000001\.6"):
        found.for_transcript(ENS_RNA + ".999")
    for name, changes in (("old-assembly", dict(reference_name="GRCm38")),
                          ("same-ids-in-rat", dict(species="rattus_norvegicus"))):
        other = ensembl_style(name, **changes)
        with pytest.raises(ValueError, match="not from this collection's reference"):
            found.for_transcript(other.transcript_by_id(ENS_RNA))
    newer = ensembl_style("newer", version="7")
    with pytest.raises(ValueError, match=r"ENSMUST00000000001\.7 differs from this reference's ENSMUST00000000001\.6"):
        found.for_transcript(newer.transcript_by_id(ENS_RNA))


def test_gene_lookup_separates_this_annotation_from_source_labels(ensembl_style):
    genome = ensembl_style("main")
    found = extra_orfs(genome, catalog=make_catalog(ENS_RNA))
    assert found.for_gene(ENS_GENE) == found.for_gene(ENS_GENE + ".6") == found.orfs
    assert found.for_gene(genome.gene_by_id(ENS_GENE)) == found.orfs
    with pytest.raises(ValueError, match=r"ENSMUSG00000000001\.5 requested"):
        found.for_gene(ENS_GENE + ".5")
    # The source catalogue calls this gene test_gene / TEST; that is not its ID here.
    assert found.for_gene("test_gene") == found.for_gene("TEST") == ()
    assert found.for_gene("test_gene", source=True) == found.for_gene("TEST", source=True) == found.orfs
    assert (found.orfs[0].gene_id, found.orfs[0].model.gene_id) == (ENS_GENE, "test_gene")
    with pytest.raises(TypeError):
        found.for_gene(genome.gene_by_id(ENS_GENE), source=True)
    other = ensembl_style("old-assembly", reference_name="GRCm38")
    with pytest.raises(ValueError, match="not from this collection's reference"):
        found.for_gene(other.gene_by_id(ENS_GENE))
    newer = ensembl_style("newer", version="7")
    with pytest.raises(ValueError, match="differs from this reference"):
        found.for_gene(newer.gene_by_id(ENS_GENE))


def test_empty_collection_keeps_its_reference(ensembl_style):
    genome = ensembl_style("main")
    nowhere = ORFRecord.from_model(ORFModel("nowhere", "ENSMUST00000099999", "g", "G", "1", "+", ((5000, 5008),), "MK",
                                            "lncRNA", "GRCm39", 10090, "mus_musculus", source="synthetic"))
    catalog = ORFCatalog([nowhere], make_catalog().metadata)
    empty = extra_orfs(genome, catalog=catalog)  # No transcript carries this path.
    assert len(empty) == 0 and len(empty.excluded) == 1 and empty.reference.genome_build == "GRCm39"
    other = ensembl_style("old-assembly", reference_name="GRCm38")
    with pytest.raises(ValueError, match="not from this collection's reference"):
        empty.for_transcript(other.transcript_by_id(ENS_RNA))
    assert empty.for_transcript(ENS_RNA) == ()
    with pytest.raises(TypeError, match="ReferenceContext"):
        ExtraORFs((), (), None, None)
    with pytest.raises(ValueError, match="bound to its reference"):
        ExtraORFs(extra_orfs(genome, catalog=make_catalog(ENS_RNA)).orfs, (), None, other_reference(empty))


def other_reference(collection):
    from dataclasses import replace

    return replace(collection.reference, genome_build="GRCm38")


def test_uninstalled_annotation_is_never_downloaded_or_indexed(tmp_path):
    genome = write_genome(tmp_path / "genome")
    genome.close()
    unindexed = Genome("GRCm39", "synthetic", "1", gtf_path_or_url=genome.gtf_path,
                       transcript_fasta_paths_or_urls=genome.transcript_fasta_paths,
                       cache_directory_path=str(tmp_path / "empty-cache"))
    unindexed.species = "mus_musculus"
    for build in (Proteome.from_genome, extra_orfs):
        with pytest.raises(ValueError, match="annotation is not indexed; install it first"):
            build(unindexed)
    assert not (tmp_path / "empty-cache").exists()


def test_extra_orfs_follow_separate_ribo_seq_and_ms_rules(genome, catalog):
    base = Proteome.from_genome(genome)

    def added(rule):
        return sorted(e.accession for e in base.with_extra_orfs(rule, catalog=catalog) if e.kind == "extra_orf")

    assert added(RiboSeq()) == ["ribo-orf"] and added(RiboSeq(min_studies=2)) == []
    assert added(MS(methods="proteomics")) == ["absent-orf", "ms-orf"]
    assert added(RiboSeq() | MS()) == ["absent-orf", "ms-orf", "ribo-orf"] and added(RiboSeq() & MS()) == []
    layer = base.with_extra_orfs(RiboSeq() | MS(), catalog=catalog).layers[-1]
    assert layer["evidence"] == {"any": [MS().to_dict(), RiboSeq().to_dict()]}
    assert extra_orfs(genome, RiboSeq() | MS(), catalog=catalog).evidence == RiboSeq() | MS()
    with pytest.raises(ValueError, match="on the rule itself"):
        base.with_extra_orfs(MS(), catalog=catalog, require_unique=True)


def test_species_without_a_catalogue_or_mismatched_catalogue_fails(genome):
    base = Proteome.from_genome(genome)
    with pytest.raises(ValueError, match=r"openprot-mouse 2.2-1 catalogue, which is not installed: build it"):
        base.with_extra_orfs()
    with pytest.raises(ValueError, match="cannot bind"):
        base.with_extra_orfs(catalog=load_catalog())


def test_a_registered_catalogue_is_the_default_for_its_species(genome, monkeypatch):
    class OneCatalogue:
        def list(self):
            return (SimpleNamespace(catalog_id="mouse-orfs", version="1", genome_build="GRCm39",
                                    references=(SimpleNamespace(species="mus_musculus"),)),)

        def info(self, name, version=None):
            return self.list()[0]

        def published(self, name, version=None):
            return True

        def inspect(self, name, version=None, verify_files=True):
            return SimpleNamespace(status="available")

        def open(self, name, version=None):
            return make_catalog()

    monkeypatch.setattr("orfdata.catalogs.CatalogRegistry", OneCatalogue)
    base = Proteome.from_genome(genome)
    assert base.with_extra_orfs(None).identity == base.with_extra_orfs(None, catalog=make_catalog()).identity


def test_bundled_human_catalogue_is_the_default_for_grch38(tmp_path):
    human = write_genome(tmp_path / "human", "GRCh38", "homo_sapiens")
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


def test_reference_bound_ncorf_only_proteome_skips_annotated_proteins(genome, catalog, tmp_path, monkeypatch):
    transcriptome = ReferenceTranscriptome.from_genome(genome)
    assert transcriptome.reference == ReferenceContext.from_genome(genome)
    with monkeypatch.context() as patch:
        patch.setattr(genome, "protein_sequence", lambda accession: 1 / 0)
        only = transcriptome.proteome(host_sources="nc_orfs", extra_orf_evidence="all", catalog=catalog)
    assert only.summary()["kinds"] == {"extra_orf": 3}
    assert only.summary()["protein_sources"] == {"nc_orf": 3}
    assert [layer["kind"] for layer in only.layers] == ["extra_orf"]
    assert {entry.accession for entry in only} == {"ms-orf", "ribo-orf", "absent-orf"}
    assert [entry.to_dict() for entry in only] == [entry.to_dict() for entry in
                                                Proteome.build(genome, extra_orf_evidence="all",
                                                               catalog=catalog) if entry.kind == "extra_orf"]
    assert [failure.orf_id for failure in only.excluded] == ["changed-orf"]
    assert only.identity == Proteome.build(genome, extra_orf_evidence="all", annotated=False,
                                           catalog=catalog).identity
    assert only.identity == Proteome.build(transcriptome, host_sources=("nc_orfs",),
                                           extra_orf_evidence="all", catalog=catalog).identity
    assert {entry.protein_source for entry in only} == {"nc_orf"}
    fasta = tmp_path / "ncorfs.fa"
    only.write_fasta(fasta)
    assert fasta.read_text().count(">") == 3
    cached = Proteome.build(genome, extra_orf_evidence="all", annotated=False,
                            catalog=catalog, cache=True, cache_root=tmp_path / "cache")
    assert cached.identity == only.identity
    assert Proteome.build(genome, extra_orf_evidence="all", catalog=catalog,
                          cache=True, cache_root=tmp_path / "cache").identity != only.identity


def test_reference_transcriptome_selects_host_protein_sources(genome, catalog):
    transcriptome = ReferenceTranscriptome.from_genome(genome)
    cds = transcriptome.proteome(host_sources=("reference_cds",))
    assert cds.summary()["protein_sources"] == {"reference_cds": 3}
    assert cds.identity == Proteome.from_genome(genome).identity
    assert "protein_source=reference_cds" in cds.entries[0].fasta_header()
    viral_only = transcriptome.proteome(host_sources=(), viral=[virus()])
    assert viral_only.summary()["protein_sources"] == {"viral": 1}
    assert viral_only.entries[0].protein_source == "viral"
    with pytest.raises(ValueError, match="Unknown host protein sources"):
        transcriptome.proteome(host_sources=("annotated",))
    with pytest.raises(ValueError, match="host_sources or annotated/extra_orfs"):
        Proteome.build(genome, host_sources=("nc_orfs",), annotated=False, catalog=catalog)
    with pytest.raises(ValueError, match="context differs"):
        Proteome.build(transcriptome, host_sources=("reference_cds",),
                       context=replace(transcriptome.reference, annotation_version="other"))


def test_multiple_synthetic_viral_sources_share_one_proteome(genome):
    second = ViralProteome("synthetic-virus-two", 999998, {10090}, (("toy-protein-two", "MAA"),),
                           "synthetic fixture", "test-declared host range")
    combined = Proteome.from_genome(genome).with_viral([virus(), second])
    assert combined.summary()["kinds"] == {"annotated": 3, "viral": 2}
    assert {entry.source for entry in combined if entry.kind == "viral"} == {
        "viral:synthetic-virus", "viral:synthetic-virus-two"}
    assert combined.by_accession()["toy-protein-two"].sequence == "MAA"


def test_bundled_virus_builds_into_one_proteome_without_manual_conversion(tmp_path, catalog):
    human = write_genome(tmp_path / "human", "GRCh38", "homo_sapiens")
    try:
        ref = viral_reference("NC_001526.4")
        human_catalog = ORFCatalog(
            [replace(record, genome_build="GRCh38", taxon_id=9606, species="homo_sapiens")
             for record in catalog.records],
            dict(catalog.metadata, genome_build="GRCh38", taxon_id=9606, species="homo_sapiens"))
        options = dict(evidence=None, catalog=human_catalog)
        by_id = Proteome.build(human, viral=ref.id, **options)
        by_type = Proteome.build(human, viral="HPV16", **options)
        by_object = Proteome.build(human, viral=ref, **options)
        all_bundled = Proteome.build(human, viral=True, **options)
        assert by_id.identity == by_type.identity == by_object.identity == all_bundled.identity
        assert by_id.identity == Proteome.build(human, viral="all", **options).identity
        assert by_id.identity == Proteome.build(human, viral="all", cache=True,
                                               cache_root=tmp_path / "cache", **options).identity
        assert by_id.identity == Proteome.from_genome(human).with_extra_orfs(
            None, catalog=human_catalog).with_viral("HPV16").identity
        assert by_id.summary()["kinds"] == {"annotated": 3, "extra_orf": 3, "viral": 9}
        assert len(by_id.by_accession()) == len(by_id) == 15
        viral = [entry for entry in by_id if entry.kind == "viral"]
        assert {entry.accession for entry in viral} == {product.protein_id for product in ref.products}
        assert all(entry.taxon_id == ref.taxon_id and entry.source == "viral:HPV16@NC_001526.4"
                   for entry in viral)
        assert by_id.layers[-1]["provenance"] == ref.source_accession
        assert by_id.layers[-1]["source_sha256"] == ref.source_sha256
        fasta = tmp_path / "combined.fa"
        by_id.write_fasta(fasta)
        assert fasta.read_text().count(">") == 15
        product = next(p for p in ref.products if p.protein_id == "NP_041327.2")
        assert (f">{product.protein_id} kind=viral protein_source=viral taxon=333760 "
                f"source=viral:HPV16@NC_001526.4 gene_name={product.gene} "
                f"product={product.product.replace(' ', '_')}\n") in fasta.read_text()
        assert by_id.by_accession()[product.protein_id].gene_name == product.gene
        assert by_id.by_accession()[product.protein_id].viral_product == product.product
        snapshot = tmp_path / "combined.json.gz"
        by_id.save(snapshot)
        assert Proteome.load(snapshot).by_accession()[product.protein_id].viral_product == product.product
        assert ViralProteome.from_reference(ref).proteins == tuple(
            (orf.orf_id, orf.protein_sequence) for orf in ref.orfs)
        gdc = Proteome.build(human, viral="GRCh38.d1.vd1:HPV16", **options)
        assert gdc.identity != by_id.identity
        assert all(entry.source == "viral:HPV16@GRCh38.d1.vd1:HPV16"
                   for entry in gdc if entry.kind == "viral")
        with pytest.raises(ValueError, match="No bundled viral reference"):
            Proteome.build(human, viral="HPV18", **options)
        with pytest.raises(ValueError, match="Choose viral='all'"):
            Proteome.build(human, viral=["all", ref.id], **options)
        alias = Proteome.build(human, extra_orf_evidence=None, catalog=human_catalog, viral="HPV16")
        assert alias.identity == by_id.identity
        with pytest.raises(ValueError, match="extra_orf_evidence or evidence"):
            Proteome.build(human, evidence="ribosome", extra_orf_evidence="all", catalog=human_catalog)
        with pytest.raises(ValueError, match="extra_orf_evidence or evidence"):
            Proteome.build(human, evidence="ms", extra_orf_evidence="all", catalog=human_catalog)
    finally:
        human.close()


def test_bundled_viral_reference_rejects_an_unrecorded_host(genome):
    with pytest.raises(ValueError, match="no recorded host matching mus_musculus"):
        Proteome.build(genome, extra_orfs=False, viral="NC_001526.4")
    with pytest.raises(ValueError, match="No bundled viral references for host taxon 10090"):
        Proteome.build(genome, extra_orfs=False, viral="all")


def test_cli_adds_bundled_viral_proteins_to_fasta(tmp_path, monkeypatch, capsys):
    from orfdata import proteome_cli

    human = write_genome(tmp_path / "human", "GRCh38", "homo_sapiens")
    try:
        monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda release, species: human)
        fasta = tmp_path / "combined.fa"
        run(["proteome", "--species", "human", "--release", "1", "--extra-orf-evidence", "none",
             "--viral", "all", "--output", str(fasta)])
        assert json.loads(capsys.readouterr().err)["kinds"] == {"annotated": 3, "viral": 9}
        assert fasta.read_text().count(">") == 12
        assert "source=viral:HPV16@NC_001526.4" in fasta.read_text()
    finally:
        human.close()


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
    assert (">ms-orf kind=extra_orf protein_source=nc_orf taxon=10090 source=orfdata-catalog:fixture-v1 gene=test_gene gene_name=TEST_GENE "
            "transcript=test_rna.1 contig=1 transcript_biotype=lncRNA orf_type=lncRNA evidence=proteomics transcripts=2") in records
    # An ORF whose source transcript is gone names it, beside the transcript that carries it.
    assert any(h.startswith(">absent-orf ") and "transcript=test_rna.1 " in h and "source_transcript=missing_rna" in h
               for h in records)


def test_cli_writes_fasta_and_reports_layers(genome, catalog, tmp_path, capsys, monkeypatch):
    from orfdata import proteome_cli

    monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda release, species: genome)
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog.to_dict()))
    output = tmp_path / "proteome.fa"
    run(["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path), "--output", str(output)])
    summary = json.loads(capsys.readouterr().err)
    assert summary["kinds"] == {"annotated": 3, "extra_orf": 2}
    assert [e["orf_id"] for e in summary["excluded_orfs"]] == ["changed-orf"]
    assert output.read_text().count(">") == 5
    nc_only = tmp_path / "ncorfs.fa"
    run(["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path),
         "--extra-orf-evidence", "all", "--host-sources", "nc-orfs", "--output", str(nc_only)])
    assert json.loads(capsys.readouterr().err)["kinds"] == {"extra_orf": 3}
    assert nc_only.read_text().count(">") == 3
    run(["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path),
         "--extra-orf-evidence", "all", "--no-annotated", "--output", str(nc_only)])
    assert json.loads(capsys.readouterr().err)["protein_sources"] == {"nc_orf": 3}
    snapshot = tmp_path / "proteome.json.gz"
    run(["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path), "--output", str(output),
         "--snapshot", str(snapshot)])
    assert Proteome.load(snapshot).identity == json.loads(capsys.readouterr().err)["identity"]
    run(["proteome", "--species", "mouse", "--release", "1", "--extra-orf-evidence", "none"])
    captured = capsys.readouterr()
    assert captured.out.count(">") == 3 and json.loads(captured.err)["kinds"] == {"annotated": 3}
    with pytest.raises(SystemExit):
        run(["proteome", "--species", "mouse", "--release", "1"])  # The mouse catalogue isn't built.
    assert "openprot-mouse 2.2-1 catalogue, which is not installed" in capsys.readouterr().err
    rules = ["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path)]
    for options, kinds in ((["--min-ribo-seq-studies", "1", "--min-ms-peptides", "1", "--either"], 3),
                           (["--min-ribo-seq-studies", "1", "--min-ms-peptides", "1"], 0),
                           (["--min-ms-peptides", "1", "--ms-methods", "proteomics", "--require-resolved"], 2)):
        run(rules + options)
        summary = json.loads(capsys.readouterr().err)
        assert summary["kinds"].get("extra_orf", 0) == kinds
    for options, message in ((["--either", "--min-ms-peptides", "1"], "--either needs both"),
                             (["--extra-orf-evidence", "ms", "--min-ms-peptides", "2"], "not both"),
                             (["--min-ribo-seq-studies", "2", "--require-unique"], "need an MS requirement"),
                             (["--min-ms-peptides", "0"], "min_peptides must be a positive integer"),
                             (["--ms-methods", "ribosome"], "MS methods must be drawn from"),
                             (["--ms-methods", ""], "MS methods must be drawn from")):
        with pytest.raises(SystemExit):
            run(rules + options)
        assert message in capsys.readouterr().err
    for count in (["--min-ms-peptides", "0"], ["--min-ribo-seq-studies", "0"]):  # Zero is still a requirement.
        with pytest.raises(SystemExit):
            run(["proteome", "--species", "mouse", "--release", "1", "--extra-orf-evidence", "none"] + count)
        assert "Extra-ORF options require" in capsys.readouterr().err


def test_build_reuses_an_identical_cached_proteome(genome, catalog, tmp_path, monkeypatch, capsys):
    from orfdata import proteome_cli

    built = Proteome.build(genome, "ms", catalog=catalog, cache=True, cache_root=tmp_path)
    assert built.identity == Proteome.from_genome(genome).with_extra_orfs("ms", catalog=catalog).identity
    calls = []
    original = Proteome.from_genome.__func__
    monkeypatch.setattr(Proteome, "from_genome", classmethod(lambda cls, *a, **k: calls.append(1) or original(cls, *a, **k)))
    again = Proteome.build(genome, "ms", catalog=catalog, cache=True, cache_root=tmp_path)
    assert calls == [] and again.identity == built.identity and again.genome is genome
    assert [e.accession for e in again] == [e.accession for e in built]
    assert again.with_extra_orfs("ribosome", catalog=catalog).layers[-1]["added"] == 1  # Still bound to the genome.
    Proteome.build(genome, MS(min_peptides=2), catalog=catalog, cache=True, cache_root=tmp_path)
    assert calls == [1]  # Another rule is another build.
    assert len(Proteome.build(genome, extra_orfs=False).entries) == 3
    monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda release, species: genome)
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog.to_dict()))
    command = ["proteome", "--species", "mouse", "--release", "1", "--catalog", str(path), "--cache",
               "--cache-root", str(tmp_path)]
    run(command)
    first = json.loads(capsys.readouterr().err)["identity"]
    run(command)
    assert json.loads(capsys.readouterr().err)["identity"] == first == built.identity and calls == [1, 1]


def test_cache_keys_viral_layers_and_annotation_files(genome, catalog, tmp_path, monkeypatch):
    with_virus = Proteome.build(genome, catalog=catalog, viral=[virus()], cache=True, cache_root=tmp_path)
    assert with_virus.identity == Proteome.build(genome, catalog=catalog, viral=[virus()]).identity
    # Same names, different annotation files: never the same cache entry.
    twins = write_genome(tmp_path / "twins", twin=True)
    try:
        assert "twin_protein.1" in {e.accession for e in Proteome.build(twins, extra_orfs=False, cache=True,
                                                                         cache_root=tmp_path)}
        plain = Proteome.build(genome, extra_orfs=False, cache=True, cache_root=tmp_path)
        assert "twin_protein.1" not in {e.accession for e in plain}
    finally:
        twins.close()
    with pytest.raises(ValueError, match="needs extra_orfs=True"):
        Proteome.build(genome, catalog=catalog, extra_orfs=False)


def test_cache_races_and_unreadable_stores(genome, catalog, tmp_path, monkeypatch):
    from orfdata import proteome as module

    calls = []

    def make():
        calls.append(1)
        return Proteome.from_genome(genome)

    key = dict(test="race")
    # Another program installs the identical build while this one is building.
    racing = module._cached_snapshot(key, lambda: (module._cached_snapshot(key, make, tmp_path), make())[1], tmp_path)
    assert Proteome.load(racing).identity == Proteome.from_genome(genome).identity and len(calls) == 2
    import datacache

    monkeypatch.setattr(datacache, "inspect_bundle", lambda store: SimpleNamespace(status="inaccessible"))
    with pytest.raises(PermissionError, match="inaccessible"):
        module._cached_snapshot(dict(test="locked"), make, tmp_path)
    assert len(calls) == 2  # Nothing was built for an unreadable store.


def test_cache_root_is_where_catalogues_are_found_and_killed_builds_are_swept(genome, tmp_path, monkeypatch):
    import os

    from orfdata import proteome as module

    seen = []
    monkeypatch.setattr(module, "load_catalog", lambda **kw: seen.append(kw["registry"]) or make_catalog())
    Proteome.build(genome, cache_root=tmp_path / "root")  # Without cache=True too: lookups use it.
    assert seen[0]._cache_root == tmp_path / "root"
    proteomes = tmp_path / "root" / "orfdata" / "proteomes"
    proteomes.mkdir(parents=True)
    stale, fresh = proteomes / ".partial-old", proteomes / ".partial-new"
    for path in (stale, fresh):
        path.mkdir()
    os.utime(stale, (0, 0))
    Proteome.build(genome, cache=True, cache_root=tmp_path / "root")
    assert not stale.exists() and fresh.exists()


def test_a_species_without_a_catalogue_says_how_to_build_annotated_proteins_only():
    with pytest.raises(ValueError, match=r"pass extra_orfs=False \(orfdata proteome --extra-orf-evidence none\)"):
        load_catalog(species="dog")


def test_a_start_codon_that_reads_differently_is_named_as_such():
    from orfdata.proteome import _binding_reason

    error = ValueError("Unsupported ORF initiation codon ATA for NCBI table 1; an explicit initiation model is required")
    assert _binding_reason(error) == "Start codon reads ATA (not a start in NCBI table 1)"
    assert _binding_reason(ValueError("ORF protein mismatch")) == "ORF protein mismatch"
