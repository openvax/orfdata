"""Native protein/source contexts map explicitly without acquiring false identity."""

import json
import sqlite3
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from orfdata import (
    Assay,
    DatasetRegistry,
    MappedMeasurement,
    Measurement,
    QuantitativeDataset,
    ReferenceContext,
    ReferenceMapping,
    Sample,
    ScopedTarget,
)
from orfdata.cli import run


@pytest.fixture
def mapped(experiment):
    manifest = experiment[0]
    native, samples, assays, observations, mappings = [], [], [], [], []
    for index, destination in enumerate(manifest.references):
        source = ReferenceContext(destination.taxon_id, destination.species, None, "uniprot",
                                  sequence_sha256=str(index + 1) * 64)
        native.append(source)
        samples.append(Sample("sample", source, attributes={"sample_breed": "source reported"}))
        assays.append(Assay("ms", "sample", "digest_ms", source,
                            search_reference=source.sequence_sha256))
        target = ScopedTarget(source, "protein", "P1", "uniprot:source")
        observations.append(Measurement("value", target, "sample", "ms", "protein_abundance",
                                        index + 10, "iBAQ", "fixture", "1",
                                        provenance={"source_asset": "input.tsv", "row": index + 1}))
        annotation = ScopedTarget(destination, "protein", "P1", "ensembl:110")
        mappings.append(ReferenceMapping("map", "fixture", "1", target, destination, (annotation,),
                                         "unique", "exact_sequence_search",
                                         {"sequence": {"source_sha256": "a" * 64,
                                                       "target_sha256": "a" * 64,
                                                       "all_candidates_checked": True}},
                                         provenance={"source_asset": "input.tsv"}))
    manifest = replace(manifest, references=(*native, *manifest.references))
    return manifest, dict(samples=samples, assays=assays, measurements=observations, reference_mappings=mappings)


def write(tmp_path, mapped, **changes):
    manifest, records = mapped
    return QuantitativeDataset.create(tmp_path / "mapped.sqlite", manifest, **dict(records, **changes))


def test_published_reference_keys_stay_stable():
    fixture = Path(__file__).parent / "data/references/context-keys-v0.4.json"
    for row in json.loads(fixture.read_text()):
        context = ReferenceContext(**row["context"])
        assert context.key == row["key"]
        assert context.is_genomic


def test_incomplete_source_requires_pins_and_never_binds_a_genome():
    with pytest.raises(ValueError, match="pinned"):
        ReferenceContext(9615, "dog", None, "uniprot")
    source = ReferenceContext(9615, "dog", None, "uniprot", sequence_sha256="a" * 64)
    assert not source.is_genomic
    assert source.key != replace(source, sequence_sha256="b" * 64).key
    assert source.to_dict()["genome_build"] is None
    assert ReferenceContext(**source.to_dict()) == source
    genome = SimpleNamespace(species="dog", reference_name="ROS_Cfam_1.0", annotation_name="ensembl")
    with pytest.raises(ValueError, match="cannot bind"):
        source.validate_genome(genome)
    with pytest.raises(ValueError, match="cannot bind"):
        ReferenceContext.from_genome(genome, context=source)
    study = ReferenceContext(4932, "yeast", None, None, source_sha256="c" * 64)
    assert study.annotation_name is None and not study.is_genomic
    assert replace(study, source_sha256="d" * 64).key != study.key
    with pytest.raises(ValueError, match="assembly unavailable"):
        ReferenceContext.from_genome(SimpleNamespace(species="dog", taxon_id=9615,
                                                      reference_name=None, annotation_name="uniprot",
                                                      sequence_sha256="a" * 64))


def test_four_species_mapping_retains_native_observations(tmp_path, mapped):
    dataset = write(tmp_path, mapped)
    assert len(list(dataset.reference_mappings())) == 4
    for index, mapping in enumerate(mapped[1]["reference_mappings"]):
        target = mapping.targets[0]
        assert list(dataset.measurements(target=target)) == []
        result = list(dataset.mapped_measurements(target))
        assert len(result) == 1
        view = result[0]
        assert view.measurement == mapped[1]["measurements"][index]
        assert view.measurement.target.reference.genome_build is None
        assert view.requested_target.reference.genome_build is not None
        assert view.attribution == "unique"
        assert view.mapping.source.reference.taxon_id == target.reference.taxon_id
        assert MappedMeasurement.from_dict(view.to_dict()) == view
        assert ReferenceMapping.from_dict(mapping.to_dict()) == mapping
        assert list(dataset.mapped_measurements(target, require_sequence=True)) == result
    with pytest.raises(ValueError, match="cannot override"):
        list(dataset.mapped_measurements(target, reference=mapping.source.reference))


def test_identifier_lookup_is_not_sequence_validation(tmp_path, mapped):
    mapping = replace(mapped[1]["reference_mappings"][0], method="identifier_lookup",
                      validation={"identifier": "unique source ID", "sequence": "not checked"})
    dataset = write(tmp_path, mapped, reference_mappings=[mapping])
    assert len(list(dataset.mapped_measurements(mapping.targets[0]))) == 1
    assert list(dataset.mapped_measurements(mapping.targets[0], require_sequence=True)) == []
    assert not mapping.sequence_validated


def test_shared_group_candidate_is_context_not_an_independent_protein(tmp_path, mapped):
    mapping = mapped[1]["reference_mappings"][0]
    group = replace(mapping.source, kind="protein_group")
    with pytest.raises(ValueError, match="group observation"):
        replace(mapping, source=group)
    unresolved = replace(mapping, source=group, status="unresolved", reason="reported representative only",
                          validation={"membership": "not independently resolved"})
    observation = replace(mapped[1]["measurements"][0], target=group)
    dataset = write(tmp_path, mapped, measurements=[observation], reference_mappings=[unresolved])
    target = unresolved.targets[0]
    assert list(dataset.mapped_measurements(target)) == []
    view = next(dataset.mapped_measurements(target, mapping_statuses=["unresolved"]))
    assert view.attribution == "group_membership_context"
    assert view.measurement.target.kind == "protein_group"
    assert view.measurement.value == observation.value
    assert not hasattr(view, "value")


def test_mapping_ambiguity_keeps_all_candidates(tmp_path, mapped):
    mapping = mapped[1]["reference_mappings"][0]
    candidate = replace(mapping.targets[0], identifier="P2")
    mapping = replace(mapping, targets=(*mapping.targets, candidate), status="ambiguous")
    dataset = write(tmp_path, mapped, reference_mappings=[mapping])
    for target in mapping.targets:
        assert list(dataset.mapped_measurements(target)) == []
        view = next(dataset.mapped_measurements(target, mapping_statuses=["ambiguous"]))
        assert len(view.mapping.targets) == 2
        assert view.attribution == "ambiguous"
    with pytest.raises(ValueError, match="candidate count"):
        replace(mapping, status="unique")


def test_gene_rna_cannot_map_to_orf_and_orfs_need_direct_validation(mapped):
    mapping = mapped[1]["reference_mappings"][0]
    orf = replace(mapping.targets[0], kind="orf")
    with pytest.raises(ValueError, match="target resolution"):
        replace(mapping, source=replace(mapping.source, kind="gene"), targets=(orf,))
    with pytest.raises(ValueError, match="sequence, frame and position"):
        replace(mapping, targets=(orf,))
    validated = replace(mapping, targets=(orf,), validation=dict(mapping.validation,
                                                                  frame={"method": "verified reading frame"},
                                                                  position={"method": "verified genomic origin"}))
    assert validated.sequence_validated
    with pytest.raises(ValueError, match="sequence, frame and position"):
        replace(validated, validation=dict(validated.validation, sequence={"source_sha256": "a" * 64,
                                                                          "target_sha256": "b" * 64,
                                                                          "all_candidates_checked": True}))


def test_cross_species_strains_and_wrong_release_are_rejected(tmp_path, mapped):
    mapping, other = mapped[1]["reference_mappings"][:2]
    with pytest.raises(ValueError, match="cross species"):
        replace(mapping, target_reference=other.target_reference, targets=other.targets)
    with pytest.raises(ValueError, match="destination reference"):
        replace(mapping, targets=(replace(mapping.targets[0], reference=replace(mapping.target_reference,
                                                                                annotation_version="111")),))
    source = replace(mapping.source, reference=replace(mapping.source.reference, strain="strain-A"))
    destination = replace(mapping.target_reference, strain="strain-B")
    with pytest.raises(ValueError, match="reference strains"):
        replace(mapping, source=source, target_reference=destination,
                targets=(replace(mapping.targets[0], reference=destination),))
    with pytest.raises(ValueError, match="different dataset release"):
        write(tmp_path, mapped, reference_mappings=[replace(mapping, dataset_version="2")])
    with pytest.raises(sqlite3.IntegrityError):
        write(tmp_path, mapped, reference_mappings=[mapping, replace(mapping, mapping_id="contradictory")])


def test_unmapped_records_and_old_storage_remain_inspectable(tmp_path, mapped):
    original = mapped[1]["reference_mappings"][0]
    missing = replace(original, status="unmapped", targets=(), reason="no source sequence match")
    dataset = write(tmp_path, mapped, reference_mappings=[missing])
    assert list(dataset.reference_mappings()) == [missing]
    assert list(dataset.mapped_measurements(original.targets[0])) == []
    with sqlite3.connect(dataset.path) as connection:
        connection.execute("DROP TABLE reference_assignment")
        connection.execute("DROP TABLE reference_mapping")
        connection.execute("PRAGMA user_version=2")
    old = QuantitativeDataset(dataset.path)
    assert list(old.reference_mappings()) == []
    assert len(list(old.measurements())) == 4


def test_mapped_cli_preserves_source_identity(tmp_path, mapped, capsys):
    from test_measurements import definition

    dataset = write(tmp_path, mapped)
    metadata = definition(dataset)
    cache_root = tmp_path / "cache"
    DatasetRegistry(metadata, cache_root=cache_root).download("fixture", show_progress=False)
    path = tmp_path / "registry.json"
    path.write_text(json.dumps(metadata))
    target = mapped[1]["reference_mappings"][0].targets[0]
    run(["datasets", "--registry", str(path), "--cache-root", str(cache_root),
         "query", "fixture", "--mapped", "--require-sequence", "--reference", target.reference.key,
         "--namespace", target.namespace, "--kind", target.kind, "--target", target.identifier])
    result = json.loads(capsys.readouterr().out)
    assert result["measurement"]["target"]["reference"]["genome_build"] is None
    assert result["requested_target"]["reference"]["genome_build"] == "GRCh38"
    assert result["attribution"] == "unique"


def test_arbitrary_taxa_use_the_same_source_reference_contract():
    source = ReferenceContext(999999, "syntheticus_example", None, None, source_sha256="a" * 64)
    destination = ReferenceContext(999999, "syntheticus_example", "assembly-1", "custom", "1")
    target = ScopedTarget(destination, "gene", "gene", "custom:1")
    mapping = ReferenceMapping("id", "data", "1", ScopedTarget(source, "gene", "gene", "source:1"),
                               destination, (target,), "unique", "source_identifier_lookup",
                               {"identifier": "exact source ID", "sequence": None})
    assert mapping.source.reference.taxon_id == 999999
    assert not mapping.sequence_validated
