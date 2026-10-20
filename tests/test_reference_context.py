"""Reference identities must survive aliases without crossing species or releases."""

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from pyensembl import EnsemblRelease, Genome
from pyensembl.species import find_species_by_name

from orfdata import ORFModel, ORFRecord, ReferenceContext, ScopedTarget, load_catalog
from orfdata.model import _species
from orfdata.reference import canonical_species, species_name_for_taxon, species_taxon

SPECIES = (
    (9606, "homo_sapiens", "GRCh38"),
    (10090, "mus_musculus", "GRCm39"),
    (9615, "canis_lupus_familiaris", "ROS_Cfam_1.0"),
    (4932, "saccharomyces_cerevisiae", "R64-1-1"),
)


@pytest.mark.parametrize("taxon,aliases", [
    (9606, ("human", "Homo sapiens")),
    (10090, ("mouse", "Mus musculus")),
    (9615, ("dog", "Canis familiaris", "Canis lupus familiaris")),
    (4932, ("yeast", "budding_yeast", "Saccharomyces cerevisiae")),
    (9685, ("cat", "Felis catus")),
    (9913, ("cattle", "cow", "Bos taurus")),
    (10116, ("rat", "Rattus norvegicus")),
])
def test_taxonomy_names_round_trip_without_annotation(taxon, aliases):
    for alias in aliases:
        assert species_taxon(alias) == taxon
        assert species_name_for_taxon(taxon) == canonical_species(alias)


@pytest.mark.parametrize("taxon", [999999999, 559292])
def test_reverse_taxonomy_does_not_guess_unknown_or_strain_species(taxon):
    assert species_name_for_taxon(taxon) is None


@pytest.mark.parametrize("taxon", [None, True, False, 0, -1, 9913.0, "9913"])
def test_reverse_taxonomy_requires_a_positive_integer(taxon):
    with pytest.raises(ValueError, match="taxon_id must be a positive NCBI taxonomy ID"):
        species_name_for_taxon(taxon)


def context(taxon=9615, species="dog", genome_build="ROS_Cfam_1.0", **kwargs):
    return ReferenceContext(taxon, species, genome_build, "ensembl", "110", **kwargs)


def transcript(reference, strand="+"):
    # Deliberately reuse every local identifier and sequence across taxa.
    return SimpleNamespace(
        id="same_rna", versioned_id="same_rna.1", contig="1", strand=strand,
        exons=[SimpleNamespace(start=1, end=9)], sequence="ATGAAATAA",
        genome=SimpleNamespace(
            reference_name=reference.genome_build, species=reference.species,
            annotation_name=reference.annotation_name, annotation_version=reference.annotation_version,
        ),
    )


def model(reference, strand="+"):
    return ORFModel(
        "same_orf", "same_rna.1", "same_gene", "SAME", "1", strand, ((1, 9),),
        "MK", "ncRNA", reference.genome_build, reference.taxon_id, reference.species,
        reference_context=reference, source="synthetic:1",
    )


@pytest.mark.parametrize("taxon,species,assembly", SPECIES)
def test_pyensembl_species_reference_resolution_is_offline(taxon, species, assembly):
    genome = EnsemblRelease(110, species=species)
    resolved = ReferenceContext.from_genome(genome)
    assert resolved.taxon_id == taxon
    assert resolved.genome_build == assembly
    assert resolved.annotation_version == "110"
    assert ReferenceContext(**json.loads(json.dumps(resolved.to_dict()))) == resolved


def test_dog_alias_regression_and_no_breed_assembly_collapse():
    expected = context()
    assert _species(find_species_by_name("canis_lupus_familiaris")) == _species("canis_lupus_familiaris")
    keys = {replace(expected, species=name).key
            for name in ("dog", "canis_familiaris", "Canis lupus familiaris")}
    assert keys == {expected.key}
    assert replace(expected, genome_build="Dog10K_Boxer_Tasha").key != expected.key
    assert replace(expected, strain="boxer").key != expected.key
    bound = ORFRecord.from_model(model(expected)).bind(transcript(replace(expected, species="canis_familiaris")))
    assert bound.context.key == expected.key
    assert bound.model.species == "dog"  # Keep the source's label.


@pytest.mark.parametrize("strand", ("+", "-"))
def test_colliding_local_ids_never_bind_across_species(strand):
    contexts = [context(t, s, "same_assembly") for t, s, _ in SPECIES]
    targets = [ScopedTarget(c, "orf", "same_orf", "same_dataset:1") for c in contexts]
    assert len({t.key for t in targets}) == 4
    for first in contexts:
        record = ORFRecord.from_model(model(first, strand))
        for second in contexts:
            if first == second:
                bound = record.bind(transcript(second, strand))
                assert bound.context.key == first.key
                assert bound.genomic_positions == (tuple(range(1, 10)) if strand == "+" else tuple(range(9, 0, -1)))
            else:
                with pytest.raises(ValueError, match="species|taxon"):
                    record.bind(transcript(second, strand))


@pytest.mark.parametrize("change", (
    dict(annotation_version="111"), dict(annotation_name="refseq"),
    dict(genome_build="CanFam3.1"), dict(strain="beagle"),
))
def test_pinned_annotation_and_reference_mismatches_fail(change):
    source = context()
    record = ORFRecord.from_model(model(source))
    altered = replace(source, **change)
    t = transcript(altered)
    t.genome.strain = altered.strain
    with pytest.raises(ValueError, match="assembly|reference context"):
        record.bind(t)


def test_missing_custom_genome_species_requires_explicit_context():
    c = context()
    t = transcript(c)
    del t.genome.species
    record = ORFRecord.from_model(model(c))
    with pytest.raises(ValueError, match="explicit ReferenceContext"):
        record.bind(t)
    assert record.bind(t, context=c).context == c
    with pytest.raises(ValueError, match="reference_name"):
        record.bind(t, context=replace(c, genome_build="wrong"))


def test_custom_taxon_and_strain_are_not_forced_into_registered_species():
    c = context(123456789, "synthetic_species", "custom.1")
    t = transcript(c)
    with pytest.raises(ValueError, match="explicit ReferenceContext"):
        ORFRecord.from_model(model(c)).bind(t)
    assert ORFRecord.from_model(model(c)).bind(t, context=c).context == c
    yeast = context(4932, "yeast", "R64-1-1")
    strain = replace(yeast, strain="S288C", strain_taxon_id=559292)
    assert yeast.key != strain.key
    with pytest.raises(ValueError, match="Strain"):
        replace(context(), strain_taxon_id=559292)
    with pytest.raises(ValueError, match="separate strain"):
        context(559292, "saccharomyces_cerevisiae_s288c", "R64-1-1")


def test_target_namespace_type_and_reference_are_part_of_identity():
    first = ScopedTarget(context(), "protein", "shared", "experiment:1")
    for changed in (replace(first, namespace="experiment:2"), replace(first, kind="gene"),
                    replace(first, reference=replace(context(), annotation_version="111"))):
        assert changed.key != first.key


def test_genomic_coordinates_are_annotation_independent_but_species_scoped():
    c = context()
    other = transcript(replace(c, annotation_version="111")).genome
    with pytest.raises(ValueError, match="annotation_version"):
        c.validate_genome(other)
    c.validate_genome(other, include_annotation=False)
    other.species = "mouse"
    with pytest.raises(ValueError, match="taxon"):
        c.validate_genome(other, include_annotation=False)


def test_catalogue_preserves_context_roundtrip(tmp_path):
    from orfdata import ORFCatalog

    c = context()
    record = ORFRecord.from_model(model(c))
    catalogue = ORFCatalog([record], dict(
        taxon_id=c.taxon_id, species="canis_familiaris", genome_build=c.genome_build,
        coordinates="1-based-inclusive", dataset_version="test:1", sources=[{"id": "synthetic"}],
    ))
    path = tmp_path / "dog.json"
    path.write_text(json.dumps(catalogue.to_dict()))
    restored = load_catalog(path)
    assert restored.records == catalogue.records
    assert restored.for_transcript(transcript(c))[0].context.key == c.key


def test_cli_keeps_pinned_reference_in_json_tsv_and_fasta(tmp_path, capsys):
    import csv
    import io

    from orfdata import ORFCatalog
    from orfdata.cli import run

    c = context()
    catalogue = ORFCatalog([ORFRecord.from_model(model(c))], dict(
        taxon_id=c.taxon_id, species=c.species, genome_build=c.genome_build,
        coordinates="1-based-inclusive", dataset_version="test:1", sources=[{"id": "synthetic"}],
    ))
    path = tmp_path / "dog.json"
    path.write_text(json.dumps(catalogue.to_dict()))
    run(["--catalog", str(path), "--format", "json"])
    assert json.loads(capsys.readouterr().out)["records"][0]["reference_context"]["taxon_id"] == 9615
    run(["--catalog", str(path), "--format", "tsv"])
    row = next(csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t"))
    assert json.loads(row["reference_context_json"])["annotation_version"] == "110"
    run(["--catalog", str(path), "--format", "fasta"])
    fasta = capsys.readouterr().out
    assert f"taxon=9615 reference={c.key}" in fasta
    assert f"kind=extra_orf source={catalogue.source}" in fasta


def test_taxonomy_conflicts_and_unsupported_alias_versions_fail():
    with pytest.raises(ValueError, match="disagree"):
        context(9606, "dog")
    with pytest.raises(ValueError, match="taxonomy alias version"):
        context(taxonomy_version="future")


@pytest.mark.parametrize("species", [s for _, s, _ in SPECIES])
def test_real_reference_slices_bind_to_pyensembl(tmp_path, species):
    data = json.loads((Path(__file__).parent / "data" / "references" / (species + ".json")).read_text())
    fields = data["model"]
    fields["reference_context"] = c = ReferenceContext(**fields["reference_context"])
    record = ORFRecord.from_model(ORFModel(**fields))
    rna = data["transcript"]
    attrs = f'gene_id "{record.gene_id}"; gene_name "{record.gene_name}"; transcript_id "{rna["id"]}"; gene_biotype "{rna["biotype"]}"; transcript_biotype "{rna["biotype"]}";'
    if rna["version"] is not None:
        attrs += f' transcript_version "{rna["version"]}";'
    rows = [("gene", rna["start"], rna["end"], ""), ("transcript", rna["start"], rna["end"], "")]
    rows += [("exon", a, b, f' exon_id "e{i}"; exon_number "{i+1}";')
             for i, (a, b) in enumerate(rna["exons"])]
    gtf = tmp_path / "slice.gtf"
    gtf.write_text("".join(f"{record.contig}\tfixture\t{kind}\t{a}\t{b}\t.\t{record.strand}\t.\t{attrs}{extra}\n"
                           for kind, a, b, extra in rows))
    fasta = tmp_path / "slice.fa"
    fasta.write_text(f'>{rna["id"]}\n{rna["sequence"]}\n')
    genome = Genome(c.genome_build, c.annotation_name, c.annotation_version,
                    gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name(species)
    genome.assembly_accession = c.assembly_accession
    genome.strain = c.strain
    genome.strain_taxon_id = c.strain_taxon_id
    try:
        genome.index()
        parent = genome.transcript_by_id(rna["id"])
        bound = record.bind(parent)
        assert bound.protein_sequence == fields["protein_sequence"]
        assert bound.context.key == c.key
        assert bound.transcript is parent and parent.biotype == rna["biotype"]
        assert bound.to_dict()["sequence_fingerprints"]["protein_sha256"]
        # These fixtures are annotation/sequence checks, not new MS evidence.
        assert not bound.evidence.has_evidence("ms")
        with pytest.raises(ValueError, match="translation table"):
            replace(record, translation_table=4).bind(parent)
    finally:
        genome.close()
