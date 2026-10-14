"""Declared NCBI genetic codes bind native mitochondrial ORFs and fail closed otherwise."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from pyensembl import Genome
from pyensembl.species import find_species_by_name

from orfdata import EvidenceSummary, ORFEvidence, ORFModel, ORFRecord, ReferenceContext
from orfdata.sequence import CODONS, GENETIC_CODES, genetic_code

FIXTURES = json.loads((Path(__file__).parent / "data/references/mitochondrial.json").read_text())["examples"]


def record(example, **changes):
    fields = dict(example["model"], **changes)
    fields["reference_context"] = ReferenceContext(**fields["reference_context"])
    summaries = tuple(EvidenceSummary.from_dict(s) for s in example.get("summaries", ()))
    return ORFRecord.from_model(ORFModel(**fields), ORFEvidence(summaries=summaries))


@pytest.fixture
def native(tmp_path, request):
    """Build an actual PyEnsembl Genome from a transcript/exon slice of the pinned GTF."""
    example = request.param
    orf = record(example)
    rna = example["transcript"]
    context = orf.reference_context
    attrs = (f'gene_id "{orf.gene_id}"; gene_name "{orf.gene_name}"; transcript_id "{rna["id"]}"; '
             f'gene_biotype "{rna["biotype"]}"; transcript_biotype "{rna["biotype"]}";')
    if rna["version"] is not None:
        attrs += f' transcript_version "{rna["version"]}";'
    rows = [("gene", rna["start"], rna["end"], ""), ("transcript", rna["start"], rna["end"], "")]
    rows += [("exon", a, b, f' exon_id "e{i}"; exon_number "{i + 1}";') for i, (a, b) in enumerate(rna["exons"])]
    gtf = tmp_path / "slice.gtf"
    gtf.write_text("".join(f"{orf.contig}\tfixture\t{kind}\t{a}\t{b}\t.\t{orf.strand}\t.\t{attrs}{extra}\n"
                           for kind, a, b, extra in rows))
    fasta = tmp_path / "slice.fa"
    versioned = rna["id"] + ("" if rna["version"] is None else "." + rna["version"])
    fasta.write_text(f">{versioned}\n{rna['sequence']}\n")
    genome = Genome(context.genome_build, context.annotation_name, context.annotation_version,
                    gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name(context.species)
    genome.assembly_accession = context.assembly_accession
    genome.strain = context.strain
    genome.strain_taxon_id = context.strain_taxon_id
    genome.index()
    try:
        yield example, orf, genome.transcript_by_id(rna["id"])
    finally:
        genome.close()


def test_tables_match_ncbi_definitions():
    assert set(GENETIC_CODES) == {1, 2, 3}
    assert dict(genetic_code(1).codons) == CODONS and len(CODONS) == 64
    assert genetic_code(1).start_codons == {"ATG", "CTG", "TTG"}
    assert genetic_code(2).start_codons == {"ATT", "ATC", "ATA", "ATG", "GTG"}
    assert genetic_code(3).start_codons == {"ATA", "ATG", "GTG"}
    standard, vertebrate, yeast = (genetic_code(t).codons for t in (1, 2, 3))
    changed = {c for c in standard if standard[c] != vertebrate[c]}
    assert changed == {"AGA", "AGG", "ATA", "TGA"}  # NCBI table 2 differences
    assert (vertebrate["AGA"], vertebrate["AGG"], vertebrate["ATA"], vertebrate["TGA"]) == ("*", "*", "M", "W")
    changed = {c for c in standard if standard[c] != yeast[c]}
    assert changed == {"ATA", "CTT", "CTC", "CTA", "CTG", "TGA"}  # NCBI table 3 differences
    assert {yeast[c] for c in ("CTT", "CTC", "CTA", "CTG")} == {"T"}


@pytest.mark.parametrize("table", (0, 4, 11, "2", None))
def test_unsupported_codes_never_fall_back_to_the_standard_code(table):
    with pytest.raises(ValueError, match="translation table"):
        genetic_code(table)


@pytest.mark.parametrize("native", FIXTURES, ids=[e["name"] for e in FIXTURES], indirect=True)
def test_native_mitochondrial_orfs_bind_only_under_their_declared_code(native):
    example, orf, transcript = native
    bound = orf.bind(transcript)
    assert bound.protein_sequence == example["model"]["protein_sequence"]
    assert bound.genetic_code.table_id == example["model"]["translation_table"]
    assert bound.to_dict()["model"]["translation_table"] == example["model"]["translation_table"]
    assert bound.transcript is transcript and transcript.biotype == example["transcript"]["biotype"]
    for table in {1, 2, 3} - {orf.translation_table}:
        with pytest.raises(ValueError, match="mismatch|initiation codon"):
            replace(orf, translation_table=table).bind(transcript)
    with pytest.raises(ValueError, match="translation table"):
        replace(orf, translation_table=4).bind(transcript)


@pytest.mark.parametrize("native", FIXTURES[:1], ids=[FIXTURES[0]["name"]], indirect=True)
def test_openprot_mitochondrial_row_keeps_native_coordinates_and_counts(native):
    example, orf, transcript = native
    row = example["native"]
    assert (row["protein accession"], row["transcript accession"]) == (orf.orf_id, orf.transcript_id)
    assert example["source_line"] == 423554
    assert orf.coding_blocks == ((int(row["start genomic coordinates"]) + 1, int(row["stop genomic coordinates"]) + 1),)
    bound = orf.bind(transcript)
    assert bound.cdna_interval == (int(row["start transcript coordinates"]) - 1, int(row["stop transcript coordinates"]) - 1)
    assert bound.coding_sequence[:3] == "ATA" and bound.coding_sequence[-3:] == "AGG"
    standard = "".join(CODONS[bound.coding_sequence[i:i + 3]] for i in range(0, len(bound.coding_sequence), 3))
    assert standard[0] == "I" and "*" not in standard  # Table 1 would read ATA as I and AGG as R.
    assert {s.metric: s.value for s in bound.evidence.summaries} == {
        "MS score": int(row["MS score"]), "TE Translation Event score": int(row["TE Translation Event score"])}
    assert bound.evidence.has_evidence("ms") and not bound.evidence.has_evidence("ms", require_resolved=True)
    assert not bound.evidence.peptides  # Aggregate counts never become peptide observations.


@pytest.mark.parametrize("table,expected", ((1, "L"), (3, "T")))
def test_codon_annotation_uses_the_declared_code(table, expected):
    from orfdata import InitiationContext

    # Synthetic codon-table check, not a biological observation.
    model = ORFModel("synthetic", "rna", "gene", "GENE", "1", "+", ((100, 108),), "MT" if table == 3 else "ML",
                     "uORF", "R64-1-1", 4932, "saccharomyces_cerevisiae", translation_table=table, source="synthetic")
    context = InitiationContext("", "ATG", "CTGTAA", tuple(range(100, 109)))
    orf = ORFRecord.from_model(model, initiation_context=context)
    codon = orf.annotate_initiation_snv(103, "C", "A", genome_build="R64-1-1")["codon_consequence"]
    assert codon["reference_codon"] == "CTG" and codon["reference_codon_translation"] == expected
    assert codon["alternate_codon"] == "ATG" and codon["alternate_codon_translation"] == "M"
    with pytest.raises(ValueError, match="translation table"):
        replace(orf, translation_table=4).annotate_initiation_snv(103, "C", "A", genome_build="R64-1-1")
