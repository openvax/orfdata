"""No manufactured sample matches, ORF specificity, replicates or molecule ratios."""

import json
import sqlite3
from dataclasses import replace

import pytest

from orfdata import (
    Assay,
    AssayPair,
    DatasetRegistry,
    MatchedMeasurements,
    Measurement,
    OrthologyMapping,
    QuantitativeDataset,
    ScopedTarget,
    TargetRelation,
    compare_orthologs,
)
from orfdata.cli import run
from orfdata.pairing import PAIRING_BASES


@pytest.fixture
def paired(experiment, tmp_path):
    manifest, samples, assays, observations = experiment
    proteins, pairs, relations = [], [], []
    assays = list(assays)
    for index, reference in enumerate(manifest.references):
        target = ScopedTarget(reference, "protein", "same_protein", "source-fasta:1")
        assays.append(Assay("protein", "sample", "digest_ms", reference, run_id="run-1",
                            technical_replicate="1", channel="126", plex="plex-1", bridge="131C",
                            normalization={"method": "source TMT bridge"}))
        proteins.append(Measurement("protein-value", target, "sample", "protein", "protein_abundance",
                                    100 + index, "iBAQ", "fixture", "1",
                                    provenance={"source_asset": "input.tsv", "row": index + 9}))
        pairs.append(AssayPair("pair", reference, "fixture", "1", "rna", "protein", "same_material",
                               {"source_asset": "input.tsv", "source_pair_id": index}, "concurrent"))
        relations.append(TargetRelation("encoding", "fixture", "1", observations[index].target,
                                         target, "canonical_translation", "one_to_one",
                                         {"source_asset": "input.tsv", "method": "test fixture"}))
    values = dict(samples=samples, assays=assays, measurements=(*observations, *proteins),
                  pairings=pairs, relations=relations)
    return manifest, values


def create(tmp_path, paired, **changes):
    manifest, values = paired
    return QuantitativeDataset.create(tmp_path / "paired.sqlite", manifest, **dict(values, **changes))


def test_four_species_native_observations_and_tmt_lineage(tmp_path, paired):
    dataset = create(tmp_path, paired)
    matches = list(dataset.matched_measurements(allowed_bases=["same_material"]))
    assert len(matches) == 4
    for match in matches:
        assert match.left[0].target.kind == "gene"
        assert match.right[0].target.kind == "protein"
        assert match.left[0].target.reference.key == match.right[0].target.reference.key
        assert len(match.left) == len(match.right) == 1
        assert match.left[0].unit == "TPM" and match.right[0].unit == "iBAQ"
        assert AssayPair.from_dict(match.pairing.to_dict()) == match.pairing
        assert TargetRelation.from_dict(match.relation.to_dict()) == match.relation
    assay = next(a for a in dataset.assays() if a.kind == "digest_ms")
    assert (assay.plex, assay.channel, assay.bridge, assay.technical_replicate) == ("plex-1", "126", "131C", "1")
    with pytest.raises(TypeError):
        list(dataset.matched_measurements())
    with pytest.raises(ValueError, match="across species"):
        list(dataset.matched_measurements(allowed_bases=["same_material"], pair_ids=["pair"]))


@pytest.mark.parametrize("basis", PAIRING_BASES)
def test_pairing_basis_must_be_allowed_explicitly(tmp_path, paired, basis):
    pair = replace(paired[1]["pairings"][0], basis=basis)
    dataset = create(tmp_path, paired, pairings=[pair])
    assert len(list(dataset.matched_measurements(allowed_bases=[basis]))) == 1
    excluded = "same_material" if basis != "same_material" else "independent_cultures"
    assert list(dataset.matched_measurements(allowed_bases=[excluded])) == []
    with pytest.raises(ValueError, match="valid allowed"):
        list(dataset.matched_measurements(allowed_bases=[]))


@pytest.mark.parametrize("temporal", ("unknown", "longitudinal"))
def test_time_mismatch_requires_opt_in(tmp_path, paired, temporal):
    pair = replace(paired[1]["pairings"][0], temporal_relation=temporal)
    dataset = create(tmp_path, paired, pairings=[pair])
    assert list(dataset.matched_measurements(allowed_bases=["same_material"])) == []
    assert len(list(dataset.matched_measurements(allowed_bases=["same_material"],
                                                allowed_temporal=[temporal]))) == 1
    with pytest.raises(ValueError, match="nonzero interval"):
        replace(pair, temporal_relation="concurrent", interval=1, interval_unit="day")


def test_missing_mates_and_multiple_values_do_not_create_replicates(tmp_path, paired):
    pair = replace(paired[1]["pairings"][0], right_assay_id=None)
    original = paired[1]["measurements"][0]
    extra = replace(original, measurement_id="rna-count", unit="counts", value=99)
    dataset = create(tmp_path, paired, pairings=[pair], measurements=(*paired[1]["measurements"], extra))
    match = next(dataset.matched_measurements(allowed_bases=["same_material"]))
    assert len(match.left) == 2 and match.right == ()
    assert {v.sample_id for v in match.left} == {"sample"}
    assert {v.assay_id for v in match.left} == {"rna"}


def test_same_tissue_does_not_make_a_pair(tmp_path, paired):
    samples = tuple(replace(sample, tissue="liver") for sample in paired[1]["samples"])
    dataset = create(tmp_path, paired, samples=samples, pairings=[])
    assert list(dataset.matched_measurements(allowed_bases=PAIRING_BASES)) == []


def test_pairing_cannot_disagree_with_known_donors_or_release(tmp_path, paired):
    sample = replace(paired[1]["samples"][0], sample_id="other", biological_unit_id="other-animal")
    assays = [replace(a, sample_id="other") if a.assay_id == "protein" and a.reference == sample.reference else a
              for a in paired[1]["assays"]]
    with pytest.raises(ValueError, match="biological-unit"):
        create(tmp_path, paired, samples=(*paired[1]["samples"], sample), assays=assays)
    pair = replace(paired[1]["pairings"][0], dataset_version="different")
    with pytest.raises(ValueError, match="different dataset release"):
        create(tmp_path, paired, pairings=[pair])
    pair = replace(paired[1]["pairings"][0], right_assay_id="missing")
    with pytest.raises(sqlite3.IntegrityError):
        create(tmp_path, paired, pairings=[pair])


def test_native_ratios_keep_units_inputs_and_undefined_denominators(tmp_path, paired):
    dataset = create(tmp_path, paired)
    matches = list(dataset.matched_measurements(allowed_bases=["same_material"]))
    for match in matches:
        numerator, denominator = match.right[0], match.left[0]
        value = match.ratio(numerator, denominator, measurement_id="ratio", method="native_ratio",
                            normalization={"protein": "source bridge", "rna": "source TPM"})
        assert value.unit == "(iBAQ)/(TPM)"
        assert value.quantity == "protein_rna_ratio"
        assert value.input_measurements == (numerator.measurement_id, denominator.measurement_id)
        assert value.derivation["uncertainty"] == "not_propagated"
        if denominator.value == 0:
            assert value.value is None and value.status == "not_measured"
            assert value.derivation["unresolved_reason"] == "nonpositive_denominator"
        else:
            assert value.value == numerator.value / denominator.value
    original = matches[0]
    for status in ("below_detection", "missing", "not_measured", "not_mapped"):
        denominator = replace(original.left[0], value=None, status=status)
        match = replace(original, left=(denominator,))
        value = match.ratio(match.right[0], denominator, measurement_id="ratio", method="test",
                            normalization={"method": "unchanged"})
        assert value.value is None and value.status == "not_measured"
    logged = replace(original.left[0], scale="log", log_base=2, value=1)
    match = replace(original, left=(logged,))
    with pytest.raises(ValueError, match="must be linear"):
        match.ratio(match.right[0], logged, measurement_id="ratio", method="test",
                    normalization={"method": "unchanged"})


def test_parent_rna_and_ambiguous_groups_cannot_become_orf_ratios(paired):
    pair, relation = paired[1]["pairings"][0], paired[1]["relations"][0]
    left, right = paired[1]["measurements"][0], paired[1]["measurements"][4]
    orf = replace(right.target, kind="orf")
    relation = replace(relation, right=orf, relation="parent_rna_context", cardinality="one_to_many")
    match = MatchedMeasurements(pair, relation, (left,), (replace(right, target=orf),))
    assert match.left[0].target.kind == "gene"
    with pytest.raises(ValueError, match="Parent RNA/group"):
        match.ratio(match.right[0], left, measurement_id="ratio", method="test", normalization={"method": "none"})
    group = replace(left.target, kind="protein_group")
    group_relation = replace(relation, left=group, relation="shared_group_member")
    match = MatchedMeasurements(pair, group_relation, (replace(left, target=group),), match.right)
    with pytest.raises(ValueError, match="Parent RNA/group"):
        match.ratio(match.right[0], match.left[0], measurement_id="ratio", method="test", normalization={"method": "none"})


def test_duplicate_isoform_attribution_is_not_one_to_one(tmp_path, paired):
    relation = paired[1]["relations"][0]
    other = replace(relation, relation_id="isoform", right=replace(relation.right, identifier="isoform2"))
    with pytest.raises(ValueError, match="Conflicting one-to-one"):
        create(tmp_path, paired, relations=[relation, other])
    with pytest.raises(ValueError, match="Conflicting one-to-one"):
        create(tmp_path, paired, relations=[relation, replace(other, cardinality="one_to_many")])
    pair = paired[1]["pairings"][0]
    left, right = paired[1]["measurements"][0], paired[1]["measurements"][4]
    match = MatchedMeasurements(pair, replace(relation, cardinality="one_to_many"), (left,), (right,))
    with pytest.raises(ValueError, match="one-to-one attribution"):
        match.ratio(right, left, measurement_id="ratio", method="test", normalization={"method": "none"})


def test_gene_and_transcript_links_are_distinct_attribution_axes(tmp_path, paired):
    gene_relation = paired[1]["relations"][0]
    transcript_relation = replace(gene_relation, relation_id="transcript-encoding",
                                   left=replace(gene_relation.left, kind="transcript", identifier="rna"))
    dataset = create(tmp_path, paired, relations=[gene_relation, transcript_relation])
    assert len(list(dataset.relations())) == 2


def test_same_target_across_samples_stays_a_separate_comparison(tmp_path, paired):
    manifest, values = paired
    original = values["measurements"][0]
    sample = replace(values["samples"][0], sample_id="second", biological_unit_id="animal-2")
    assay = Assay("rna-2", "second", "rna_seq", sample.reference)
    other = replace(original, measurement_id="second-value", sample_id="second", assay_id="rna-2", value=2)
    relation = TargetRelation("same", "fixture", "1", original.target, original.target,
                               "same_target", "one_to_one", {"identity": "exact scoped target"})
    pair = replace(values["pairings"][0], right_assay_id="rna-2", basis="biological_replicate")
    dataset = create(tmp_path, paired, samples=(*values["samples"], sample),
                      assays=(*values["assays"], assay), measurements=(*values["measurements"], other),
                      pairings=[pair], relations=[relation])
    match = next(dataset.matched_measurements(allowed_bases=["biological_replicate"]))
    assert match.left[0].target == match.right[0].target
    assert match.left[0].sample_id != match.right[0].sample_id


def test_old_storage_has_no_inferred_relationships(tmp_path, paired):
    dataset = create(tmp_path, paired)
    with sqlite3.connect(dataset.path) as connection:
        connection.execute("DROP TABLE assay_pair")
        connection.execute("DROP TABLE target_relation")
        connection.execute("PRAGMA user_version=1")
    old = QuantitativeDataset(dataset.path)
    assert list(old.pairings()) == list(old.relations()) == []
    assert list(old.matched_measurements(allowed_bases=PAIRING_BASES)) == []
    assert len(list(old.measurements())) == 8


def test_orthology_keeps_both_species_and_never_transfers_abundance(tmp_path, paired):
    dataset = create(tmp_path, paired)
    observations = paired[1]["measurements"]
    mappings = [OrthologyMapping(str(index), observations[0].target, observations[index].target,
                                paired[0].sources[0], "test:1", "one_to_one",
                                {"definition": "synthetic fixture, not biological evidence"},
                                {"method": "test fixture"}) for index in (1, 2, 3)]
    comparisons = list(compare_orthologs(dataset, dataset, mappings))
    assert len(comparisons) == 3
    for index, comparison in enumerate(comparisons, 1):
        assert comparison.left == (observations[0],)
        assert comparison.right == (observations[index],)
        assert comparison.left[0].value != comparison.right[0].value
        assert OrthologyMapping.from_dict(comparison.mapping.to_dict()) == comparison.mapping
    assert list(compare_orthologs(dataset, dataset, [replace(mappings[0], cardinality="one_to_many")])) == []
    alternate = replace(mappings[0], mapping_id="other", right=replace(mappings[0].right, identifier="other"))
    with pytest.raises(ValueError, match="Conflicting one-to-one"):
        list(compare_orthologs(dataset, dataset, [mappings[0], alternate]))
    with pytest.raises(ValueError, match="sequence, frame and position"):
        replace(mappings[0], left=replace(mappings[0].left, kind="orf"), right=replace(mappings[0].right, kind="orf"))
    orf_mapping = replace(mappings[0], left=replace(mappings[0].left, kind="orf"),
                           right=replace(mappings[0].right, kind="orf"),
                           evidence={"sequence": "pinned alignment", "frame": "same frame", "position": "synteny"})
    assert list(compare_orthologs(dataset, dataset, [orf_mapping])) == []
    comparison = next(compare_orthologs(dataset, dataset, [orf_mapping], allow_orfs=True))
    assert comparison.left == comparison.right == ()  # Gene values are not ORF values.


def test_paired_table_cli_preserves_rna_protein_and_missing_ribo(tmp_path, paired, capsys):
    from test_measurements import definition

    dataset = create(tmp_path, paired)
    metadata = definition(dataset)
    cache_root = tmp_path / "cache"
    DatasetRegistry(metadata, cache_root=cache_root).download("fixture", show_progress=False)
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(metadata))
    run(["datasets", "--registry", str(registry_path), "--cache-root", str(cache_root),
         "pairs", "fixture", "--basis", "same_material"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(rows) == 4
    assert all(row["ribosome_occupancy"] == [] for row in rows)
    assert all(row["rna_abundance"][0]["unit"] == "TPM" for row in rows)
    assert all(row["protein_abundance"][0]["unit"] == "iBAQ" for row in rows)
