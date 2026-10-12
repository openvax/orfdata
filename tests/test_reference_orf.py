"""Synthetic mouse data exercise the API; they are not experimental evidence."""

import json
from dataclasses import replace

import pytest
from pyensembl import Genome

from orfdata import ORFCatalog, ORFEvidence, ORFModel, ORFRecord, load_catalog, load_uorf_catalog


@pytest.fixture
def mouse_transcript(tmp_path):
    sequence = "CCCATGAAACCCGGGTAAATGAAATAA"
    attrs = 'gene_id "test_gene"; gene_name "TEST"; transcript_id "test_rna"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";'
    gtf = tmp_path / "mouse.gtf"
    rows = []
    for feature, start, end, extra in (
        ("gene", 100, 226, ""),
        ("transcript", 100, 226, ""),
        ("exon", 100, 104, ' exon_number "1"; exon_id "exon1";'),
        ("exon", 205, 226, ' exon_number "2"; exon_id "exon2";'),
    ):
        rows.append(f"1\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}{extra}\n")
    # Sequence is 27 bases; the start codon spans the exon junction.
    assert len(sequence) == 27
    gtf.write_text("".join(rows))
    fasta = tmp_path / "mouse.fa"
    fasta.write_text(">test_rna\n" + sequence + "\n")
    genome = Genome(
        reference_name="GRCm39",
        annotation_name="synthetic",
        gtf_path_or_url=str(gtf),
        transcript_fasta_paths_or_urls=[str(fasta)],
        cache_directory_path=str(tmp_path),
    )
    genome.species = "mus_musculus"
    genome.index()
    yield genome.transcript_by_id("test_rna")
    genome.close()


def model(**kwargs):
    values = dict(
        orf_id="test-orf",
        transcript_id="test_rna",
        gene_id="test_gene",
        gene_name="TEST",
        contig="1",
        strand="+",
        coding_blocks=((103, 104), (205, 217)),
        protein_sequence="MKPG",
        orf_type="lncRNA",
        genome_build="GRCm39",
        taxon_id=10090,
        species="mus_musculus",
        translation_table=1,
        source="synthetic",
    )
    values.update(kwargs)
    return ORFModel(**values)


def catalog(record):
    return ORFCatalog(
        [record],
        dict(
            genome_build="GRCm39",
            taxon_id=10090,
            species="mus_musculus",
            coordinates="1-based-inclusive",
            dataset_version="fixture-v1",
            sources=[dict(id="synthetic", scope="test fixture; no biological observations")],
        ),
    )


def test_real_pyensembl_noncoding_transcript_with_multiple_orfs(mouse_transcript):
    first = ORFRecord.from_model(model())
    second = ORFRecord.from_model(
        model(orf_id="second", coding_blocks=((218, 226),), protein_sequence="MK")
    )
    first_catalog = catalog(first)
    combined = ORFCatalog([first, second], first_catalog.metadata)
    bound, other = combined.for_transcript(mouse_transcript)
    assert other.protein_sequence == "MK" and other.cdna_interval == (18, 27)
    assert bound.transcript is mouse_transcript
    assert not mouse_transcript.is_protein_coding and mouse_transcript.protein is None
    assert bound.cdna_interval == (3, 18)
    assert bound.coding_sequence == "ATGAAACCCGGGTAA"
    assert bound.protein_sequence == "MKPG"
    assert bound.model.orf_type == "lncRNA"
    assert not bound.evidence.has_evidence("ms")
    exported = json.loads(json.dumps(bound.to_dict()))
    assert exported["resolved_transcript"]["biotype"] == "lncRNA"
    assert exported["model"]["taxon_id"] == 10090


def test_species_and_assembly_are_not_inferred_from_shared_ids(mouse_transcript):
    with pytest.raises(ValueError, match="species"):
        ORFRecord.from_model(model(species="homo_sapiens", taxon_id=9606)).bind(mouse_transcript)
    with pytest.raises(ValueError, match="assembly"):
        ORFRecord.from_model(model(genome_build="GRCm38")).bind(mouse_transcript)
    with pytest.raises(ValueError, match="translation table"):
        ORFRecord.from_model(model(translation_table=3)).bind(mouse_transcript)


def test_custom_species_catalogue_roundtrip_and_missing_identity(tmp_path):
    original = catalog(ORFRecord.from_model(model(), ORFEvidence()))
    path = tmp_path / "mouse.json"
    data = original.to_dict()
    path.write_text(json.dumps(data))
    restored = load_catalog(path)
    assert restored.records == original.records
    assert restored.taxon_id == 10090
    del data["records"][0]["taxon_id"]
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="explicit species"):
        load_catalog(path)
    with pytest.raises(ValueError, match="species"):
        catalog(replace(original.records[0], taxon_id=9606, species="homo_sapiens"))


def test_all_primary_source_orf_classes_and_legacy_subset():
    all_orfs = load_catalog()
    upstream = load_uorf_catalog()
    assert len(all_orfs.records) == 7264 and len(upstream.records) == 3771
    for kind in ("lncRNA", "processed_transcript", "intORF", "dORF", "doORF"):
        assert all_orfs.query(biotype=kind)
    assert all_orfs.query(biotype="lncRNA", evidence="ms")
    assert all_orfs.get("c7riboseqorf56").evidence == upstream.get("c7riboseqorf56").evidence


def test_nonhuman_mhc_evidence_is_scoped_and_preserved(tmp_path):
    from orfdata import PeptideEvidence

    observation = PeptideEvidence("mhc", "MKP", ((1, 3),), "synthetic", "fixture")
    record = ORFRecord.from_model(model(), ORFEvidence(peptides=(observation,)))
    data = catalog(record)
    assert data.query(evidence="mhc") == (record,)
    assert data.query(evidence="hla") == ()
    path = tmp_path / "mouse-ms.json"
    path.write_text(json.dumps(data.to_dict()))
    assert load_catalog(path).get(record.orf_id).evidence == record.evidence
    with pytest.raises(ValueError, match="non-human"):
        catalog(replace(record, peptides=(replace(observation, method="hla"),)))
