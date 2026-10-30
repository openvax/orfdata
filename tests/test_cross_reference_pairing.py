"""RNA and MS keep independently declared references, including colliding IDs."""

import json
import sqlite3
from dataclasses import replace

import pytest

from orfdata import (
    AssayPair,
    DatasetRegistry,
    MatchedMeasurements,
    Measurement,
    MeasurementReference,
    QuantitativeDataset,
    TargetRelation,
)
from orfdata.cli import run


@pytest.fixture(params=range(4), ids=("human", "mouse", "dog", "yeast"))
def native_pair(request, experiment):
    manifest, samples, assays, observations = experiment
    index = request.param
    rna = replace(observations[index], value=25)
    reference = replace(rna.target.reference, genome_build=None,
                        annotation_name="source-fasta", annotation_version=None,
                        sequence_sha256="a" * 64)
    protein = replace(rna, target=replace(rna.target, reference=reference,
                                         kind="protein", namespace="source-fasta"),
                      quantity="protein_abundance", value=100, unit="iBAQ")
    pair = AssayPair("pair", rna.target.reference, "fixture", "1", "rna", "rna",
                     "same_material", {"source": "synthetic aliquot crosswalk"},
                     "concurrent", right_reference=reference)
    relation = TargetRelation("encoding", "fixture", "1", rna.target, protein.target,
                              "canonical_translation", "one_to_one",
                              {"source": "synthetic encoding relationship"})
    manifest = replace(manifest, references=(rna.target.reference, reference))
    values = dict(samples=(samples[index], replace(samples[index], reference=reference)),
                  assays=(assays[index], replace(assays[index], reference=reference, kind="digest_ms")),
                  measurements=(rna, protein), pairings=(pair,), relations=(relation,))
    return manifest, values


def create(tmp_path, native_pair, **changes):
    manifest, values = native_pair
    return QuantitativeDataset.create(tmp_path / "native.sqlite", manifest, **dict(values, **changes))


def test_native_references_pair_and_derive_without_identifier_collision(tmp_path, native_pair):
    dataset = create(tmp_path, native_pair)
    match, = dataset.matched_measurements(allowed_bases=["same_material"])
    rna, protein = native_pair[1]["measurements"]
    assert match.left == (rna,) and match.right == (protein,)
    assert rna.measurement_id == protein.measurement_id  # Source IDs really collide.
    assert rna.target.reference.key != protein.target.reference.key
    selected, = dataset.matched_measurements(
        allowed_bases=["same_material"], pair_ids=["pair"],
        left_target=rna.target, right_target=protein.target,
        left_filters={"unit": "TPM"}, right_filters={"unit": "iBAQ"})
    assert selected == match
    assert list(dataset.relations(right_target=protein.target)) == [match.relation]
    assert list(dataset.relations(reference=protein.target.reference, right_target=protein.target)) == []
    assert MatchedMeasurements.from_dict(json.loads(json.dumps(match.to_dict()))) == match
    ratio = match.ratio(protein, rna, measurement_id="ratio", method="native ratio",
                        normalization={"method": "none; synthetic units"})
    assert ratio.value == 4 and ratio.unit == "(iBAQ)/(TPM)"
    assert ratio.input_measurements == ()
    assert ratio.input_references == tuple(MeasurementReference.from_measurement(value)
                                           for value in (protein, rna))
    assert Measurement.from_dict(json.loads(json.dumps(ratio.to_dict()))) == ratio
    output = QuantitativeDataset.create(tmp_path / "derived.sqlite", native_pair[0],
                                         **dict(native_pair[1], measurements=(rna, protein, ratio)))
    assert ratio in tuple(output.measurements(target=protein.target))
    with sqlite3.connect(output.path) as connection:
        assert connection.execute("SELECT count(DISTINCT input_reference_key) FROM input_measurement").fetchone()[0] == 2


def test_cross_reference_pair_cli_retains_native_scopes(tmp_path, native_pair, capsys):
    from test_measurements import definition

    dataset = create(tmp_path, native_pair)
    definitions = definition(dataset)
    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps(definitions))
    cache = tmp_path / "cache"
    DatasetRegistry(definitions, cache_root=cache).download("fixture", show_progress=False)
    capsys.readouterr()
    run(["datasets", "--registry", str(registry), "--cache-root", str(cache),
         "pairs", "fixture", "--basis", "same_material", "--format", "jsonl"])
    row = json.loads(capsys.readouterr().out)
    assert row["rna_abundance"][0]["target"]["reference"]["genome_build"] is not None
    assert row["protein_abundance"][0]["target"]["reference"]["genome_build"] is None
    assert row["pairing"]["right_reference"]["sequence_sha256"] == "a" * 64
    selections = []
    for side, value in zip(("left", "right"), native_pair[1]["measurements"]):
        target = value.target
        selections.extend([f"--{side}-target", target.identifier, f"--{side}-kind", target.kind,
                           f"--{side}-namespace", target.namespace, f"--{side}-reference", target.reference.key,
                           f"--{side}-unit", value.unit])
    run(["datasets", "--registry", str(registry), "--cache-root", str(cache),
         "pairs", "fixture", "--basis", "same_material", "--pair-id", "pair"] + selections + ["--format", "jsonl"])
    assert json.loads(capsys.readouterr().out) == row


def test_cross_reference_foreign_keys_and_biological_identity(tmp_path, native_pair):
    samples = native_pair[1]["samples"]
    with pytest.raises(ValueError, match="biological-unit"):
        create(tmp_path, native_pair, samples=(samples[0], replace(samples[1], biological_unit_id="other")))
    with pytest.raises(sqlite3.IntegrityError):
        create(tmp_path, native_pair, assays=native_pair[1]["assays"][:1])
    # A same-named RNA assay cannot satisfy the missing protein assay's foreign key.
    assert not (tmp_path / "native.sqlite").exists()


def test_scoped_input_lineage_rejects_wrong_release_species_and_cycles(tmp_path, native_pair):
    rna, protein = native_pair[1]["measurements"]
    source = MeasurementReference.from_measurement(rna)
    derived = replace(protein, measurement_id="derived", origin="derived", provenance={},
                      input_references=(source,), derivation={"method": "synthetic transformation"})
    invalid = [replace(source, dataset_version="2"), replace(source, dataset_id="other"),
               replace(source, reference=replace(source.reference, taxon_id=999999, species="other")),
               MeasurementReference.from_measurement(derived)]
    for reference in invalid:
        with pytest.raises(ValueError):
            replace(derived, input_references=(reference,))
    with pytest.raises(ValueError, match="not both"):
        replace(derived, input_measurements=("other",))
    with pytest.raises(ValueError, match="unique"):
        replace(derived, input_references=(source, source))
    with pytest.raises(ValueError, match="observed"):
        replace(protein, input_references=(source,))
    with pytest.raises(ValueError, match="precede outputs"):
        create(tmp_path, native_pair, measurements=(protein, derived, rna))
    with pytest.raises(ValueError, match="precede outputs"):
        create(tmp_path, native_pair, measurements=(rna, protein, replace(
            derived, input_references=(replace(source, measurement_id="absent"),))))
    cyclic = replace(derived, input_references=(replace(source, measurement_id="later"),))
    later = replace(rna, measurement_id="later", origin="derived", provenance={},
                    input_references=(MeasurementReference.from_measurement(cyclic),),
                    derivation={"method": "cycle fixture"})
    with pytest.raises(ValueError, match="precede outputs"):
        create(tmp_path, native_pair, measurements=(rna, protein, cyclic, later))


def test_gene_group_context_is_visible_but_cannot_produce_specific_ratio(tmp_path, native_pair):
    values = native_pair[1]
    rna, protein = values["measurements"]
    group = replace(protein, target=replace(protein.target, kind="protein_group"))
    relation = replace(values["relations"][0], right=group.target,
                       relation="gene_protein_group_context", cardinality="unknown")
    dataset = create(tmp_path, native_pair, measurements=(rna, group), relations=(relation,))
    match, = dataset.matched_measurements(allowed_bases=["same_material"])
    assert match.right == (group,)
    with pytest.raises(ValueError, match="Parent RNA/group"):
        match.ratio(group, rna, measurement_id="ratio", method="invalid", normalization={"method": "none"})
    with pytest.raises(ValueError, match="unique protein"):
        replace(relation, cardinality="one_to_one")


def test_cross_species_pair_and_target_relation_rejected(native_pair):
    pair = native_pair[1]["pairings"][0]
    other = replace(pair.right_reference, taxon_id=999999, species="other")
    with pytest.raises(ValueError, match="cross species"):
        replace(pair, right_reference=other)
    relation = native_pair[1]["relations"][0]
    with pytest.raises(ValueError, match="cross species"):
        replace(relation, right=replace(relation.right, reference=other))


def test_old_pair_record_defaults_both_sides_to_original_reference(experiment):
    pair = AssayPair("old", experiment[0].references[0], "fixture", "1", "rna", "ms",
                     "same_material", {"source": "fixture"})
    record = pair.to_dict()
    del record["right_reference"]
    assert AssayPair.from_dict(record) == pair
