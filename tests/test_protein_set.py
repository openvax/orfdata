"""ProteinSet is a generic sequence input; its use is recorded by Proteome."""

from dataclasses import replace

import pytest

from orfdata import ProteinSet, Proteome, ReferenceContext


@pytest.fixture
def human():
    return ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")


@pytest.fixture
def mouse():
    return ReferenceContext(10090, "mouse", "GRCm39", "Ensembl", "115")


def test_protein_set_reads_fasta_and_adds_an_exact_reference(tmp_path, human):
    path = tmp_path / "selected.fa"
    path.write_text(">A first protein\nMPEP*\n>B\nMPEP\n")
    proteins = ProteinSet.from_fasta(
        path,
        name="selected-human-proteins",
        reference=human,
        source="curated fixture",
    )
    assert proteins.proteins == (("A", "MPEP"), ("B", "MPEP"))
    assert proteins.taxon_id == 9606
    assert proteins.source_sha256 is not None

    proteome = Proteome(human, ()).with_protein_set(proteins, role="transgene")
    assert len(proteome.proteins) == 2
    assert {protein.kind for protein in proteome.proteins} == {"protein_set"}
    assert {protein.sequence for protein in proteome.proteins} == {"MPEP"}
    assert proteome.layers[-1]["component_role"] == "transgene"
    assert proteome.layers[-1]["source_reference"] == human.key
    assert proteome.with_protein_set(proteins, role="transgene").identity == proteome.identity

    snapshot = tmp_path / "selected.json.gz"
    proteome.save(snapshot)
    loaded = Proteome.load(snapshot)
    assert loaded.identity == proteome.identity
    assert loaded.layers[-1]["provenance"] == "curated fixture"


def test_donor_set_keeps_its_source_separate_from_its_target(human, mouse, tmp_path):
    donor = ProteinSet(
        name="selected-human-transgene",
        proteins=(("DONOR-A", "MPEP"), ("DONOR-B", "MPEP")),
        source="construct manifest 7",
        taxon_id=9606,
        reference=human,
        source_sha256="b" * 64,
    )
    host = Proteome(mouse, ())
    with pytest.raises(ValueError, match="explicit target assertion"):
        host.with_protein_set(donor, role="donor")
    with pytest.raises(ValueError, match="explicit target assertion"):
        host.with_protein_set(
            donor,
            role="donor",
            target_taxon_ids=(9606,),
            target_source="construct manifest 7",
        )

    combined = host.with_protein_set(
        donor,
        role="donor",
        target_taxon_ids=(10090,),
        target_source="mouse line manifest 3",
    )
    assert len(combined.by_sequence()["MPEP"]) == 2
    assert all(protein.taxon_id == 9606 for protein in combined.proteins)
    layer = combined.layers[-1]
    assert layer["source_reference"] == human.key
    assert layer["target_taxon_ids"] == [10090]
    assert layer["target_source"] == "mouse line manifest 3"
    assert layer["source_sha256"] == "b" * 64

    path = tmp_path / "mouse-with-donor.json.gz"
    combined.save(path)
    loaded = Proteome.load(path)
    assert loaded.reference == mouse
    assert loaded.proteins[0].taxon_id == 9606
    assert loaded.layers[-1] == layer
    with pytest.raises(ValueError, match="Conflicting protein set"):
        combined.with_protein_set(
            replace(donor, proteins=(("DONOR-A", "MCHANGED"),)),
            role="donor",
            target_taxon_ids=(10090,),
            target_source="mouse line manifest 3",
        )


def test_same_protein_set_can_be_used_as_viral_or_nonviral(human):
    proteins = ProteinSet(
        name="synthetic-sequences",
        proteins=(("P1", "MAA"),),
        source="fixture",
        taxon_id=999999,
    )
    with pytest.raises(ValueError, match="sourced host taxa"):
        Proteome(human, ()).with_viral(proteins)
    viral = Proteome(human, ()).with_viral(
        proteins,
        host_taxon_ids=(9606,),
        host_source="host-range fixture",
    )
    assert viral.proteins[0].kind == "viral"
    assert viral.layers[-1]["role"] == "viral"
    assert viral.with_viral(
        proteins,
        host_taxon_ids=(9606,),
        host_source="host-range fixture",
    ).identity == viral.identity

    component = Proteome(human, ()).with_protein_set(
        proteins,
        role="custom",
        target_taxon_ids=(9606,),
        target_source="selection fixture",
    )
    assert component.proteins[0].kind == "protein_set"
    assert component.layers[-1]["role"] == "component"
