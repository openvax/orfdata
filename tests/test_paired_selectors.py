"""Exact native selections keep scientific scope while skipping unrelated data."""

import gzip
import json
from dataclasses import replace
from pathlib import Path

import pytest
from test_pairing import create
from test_pairing import paired as paired
from test_wang import build
from test_wang import source as source

from orfdata import DatasetRegistry, Measurement, QuantitativeDataset, ScopedTarget, TargetRelation
from orfdata.cli import run
from orfdata.pairing import PAIRING_BASES


def test_pyensembl_gene_query_reads_only_its_native_pair(source, tmp_path, monkeypatch):
    from pyensembl import Genome

    import orfdata.dataset as storage

    dataset = build(source, tmp_path / "data")
    reference, = dataset.manifest.references
    annotation = tmp_path / "annotation.gtf"
    annotation.write_bytes(gzip.decompress(source[1]["annotation.gtf.gz"].read_bytes()))
    genome = Genome(reference_name="GRCh38", annotation_name="ensembl", annotation_version=83,
                    gtf_path_or_url=str(annotation), cache_directory_path=str(tmp_path / "genome"))
    genome.species = "human"
    genome.index()
    target = ScopedTarget.from_pyensembl(genome.gene_by_id("ENSG00000004059"),
                                        kind="gene", namespace="ensembl", context=reference)
    relation, = dataset.relations(left_target=target, right_target=target)
    pair_id = "rna_genes:protein_genes:V122"
    expected, = dataset.matched_measurements(allowed_bases=["adjacent_donor_tissue"],
                                             pair_ids=[pair_id], relation_ids=[relation.relation_id])
    decoded = {"relations": 0, "observations": 0}
    original_relation = TargetRelation.from_dict

    def relation_from_dict(cls, value):
        decoded["relations"] += 1
        return original_relation(value)

    def observation(**values):
        decoded["observations"] += 1
        return Measurement(**values)

    monkeypatch.setattr(TargetRelation, "from_dict", classmethod(relation_from_dict))
    monkeypatch.setattr(storage, "Measurement", observation)
    match, = dataset.matched_measurements(
        allowed_bases=["adjacent_donor_tissue"], pair_ids=[pair_id],
        left_target=target, right_target=target,
        left_filters={"quantity": "rna_abundance", "unit": "FPKM", "scale": "linear"},
        right_filters={"quantity": "protein_abundance", "unit": "iBAQ"})
    assert match == expected
    assert decoded == {"relations": 1, "observations": 2}
    assert match.left[0].value == 51.09 and match.right[0].value == 112559620.9
    assert match.pairing.basis == "adjacent_donor_tissue" and match.relation.relation == "same_target"


def test_native_group_selection_does_not_upgrade_attribution(source, tmp_path):
    dataset = build(source, tmp_path)
    reference, = dataset.manifest.references
    gene = ScopedTarget(reference, "gene", "ENSG00000004059", "ensembl")
    group = next(r for r in dataset.relations(left_target=gene) if r.relation == "gene_protein_group_context")
    options = dict(allowed_bases=["adjacent_donor_tissue"], pair_ids=["rna_genes:protein_groups:V122"])
    original, = dataset.matched_measurements(**options, relation_ids=[group.relation_id])
    selected, = dataset.matched_measurements(**options, left_target=gene, right_target=group.right,
                                             left_filters={"unit": "FPKM"}, right_filters={"unit": "iBAQ"})
    assert selected == replace(original, right=tuple(v for v in original.right if v.unit == "iBAQ"))
    assert selected.right[0].target.kind == "protein_group"
    with pytest.raises(ValueError, match="group membership"):
        selected.ratio(selected.right[0], selected.left[0], measurement_id="invalid",
                       method="invalid", normalization={"scope": "test"})
    assert list(dataset.matched_measurements(allowed_bases=["same_material"], left_target=gene)) == []


def test_target_selection_stays_separate_for_four_species(paired, tmp_path):
    dataset = create(tmp_path, paired)
    relations = paired[1]["relations"]
    for relation in relations:
        matches = list(dataset.matched_measurements(allowed_bases=["same_material"],
                                                    left_target=relation.left, right_target=relation.right,
                                                    pair_ids=["pair"], relation_ids=[relation.relation_id]))
        assert len(matches) == 1 and matches[0].relation == relation
        assert matches[0].left[0].target == relation.left
    with pytest.raises(ValueError, match="same species"):
        list(dataset.matched_measurements(allowed_bases=["same_material"],
                                           left_target=relations[0].left, right_target=relations[1].right))
    with pytest.raises(ValueError, match="Left target"):
        dataset.relations(reference=relations[1].left.reference, left_target=relations[0].left)
    wrong = replace(relations[0].left, reference=replace(relations[0].left.reference, annotation_version="other"))
    with pytest.raises(ValueError, match="not present"):
        dataset.relations(left_target=wrong)
    with pytest.raises(ValueError, match="across species"):
        list(dataset.matched_measurements(allowed_bases=["same_material"],
                                           right_target=relations[0].right, pair_ids=["pair"]))
    with pytest.raises(TypeError, match="ScopedTarget"):
        dataset.relations(left_target="same_gene")


def test_native_filters_keep_multiple_values_and_missing_mates(paired, tmp_path):
    first = paired[1]["measurements"][0]
    extra = replace(first, measurement_id="counts", value=99, unit="counts")
    logged = replace(first, measurement_id="log", value=2, scale="log", log_base=2)
    dataset = create(tmp_path, paired, measurements=(*paired[1]["measurements"], extra, logged))
    options = dict(allowed_bases=["same_material"], left_target=first.target)
    all_values, = dataset.matched_measurements(**options)
    assert len(all_values.left) == 3
    selected, = dataset.matched_measurements(**options, left_filters={"unit": "TPM", "scale": "log", "log_base": 2})
    assert selected.left == (logged,) and selected.right == all_values.right
    absent, = dataset.matched_measurements(**options, right_filters={"unit": "copies/cell"})
    assert absent.right == () and absent.left == all_values.left
    empty, = dataset.matched_measurements(**options, left_filters={"unit": "unavailable"}, right_filters={"unit": "unavailable"})
    assert empty.left == empty.right == () and empty.pairing == all_values.pairing
    assert list(dataset.matched_measurements(**options, relation_ids=[])) == []
    ids = [f"missing:{i}" for i in range(1500)]
    selected, = dataset.matched_measurements(**options, pair_ids=ids + ["pair"],
                                             relation_ids=ids + [all_values.relation.relation_id])
    assert selected == all_values


@pytest.mark.parametrize("filters,error", [
    ({"unit": ""}, ValueError), ({"scale": "other"}, ValueError),
    ({"log_base": 1}, ValueError), ({"scale": "linear", "log_base": 2}, ValueError),
    ({"quantity": "other"}, ValueError), ({"status": "other"}, ValueError),
    ({"reference": "other"}, ValueError), ({"assay_id": "other"}, ValueError),
    ([], TypeError),
])
def test_invalid_native_filters_fail_even_for_empty_legacy_data(filters, error):
    dataset = QuantitativeDataset(Path(__file__).parent / "data/storage/v1.sqlite")
    for side in ("left_filters", "right_filters"):
        with pytest.raises(error):
            list(dataset.matched_measurements(allowed_bases=PAIRING_BASES, **{side: filters}))


@pytest.mark.parametrize("version", (1, 2, 3))
def test_paired_selectors_read_legacy_releases(version):
    dataset = QuantitativeDataset(Path(__file__).parent / f"data/storage/v{version}.sqlite")
    matches = list(dataset.matched_measurements(allowed_bases=PAIRING_BASES))
    if version == 1:
        target = next(dataset.measurements()).target
        assert list(dataset.matched_measurements(allowed_bases=PAIRING_BASES, left_target=target)) == []
    else:
        assert matches
        for match in matches:
            selected, = dataset.matched_measurements(
                allowed_bases=PAIRING_BASES, pair_ids=[match.pairing.pair_id],
                relation_ids=[match.relation.relation_id],
                left_target=match.relation.left, right_target=match.relation.right)
            assert selected == match
            filtered, = dataset.matched_measurements(
                allowed_bases=PAIRING_BASES, pair_ids=[match.pairing.pair_id],
                relation_ids=[match.relation.relation_id],
                left_target=match.relation.left, right_target=match.relation.right,
                left_filters={"unit": match.left[0].unit})
            assert filtered == replace(match, left=tuple(v for v in match.left if v.unit == match.left[0].unit))
    with pytest.raises(TypeError, match="sequences"):
        list(dataset.matched_measurements(allowed_bases=PAIRING_BASES, pair_ids="pair"))


def test_cli_scoped_target_and_native_filters(paired, tmp_path, capsys):
    from test_measurements import definition

    dataset = create(tmp_path, paired)
    definitions = definition(dataset)
    registry_file = tmp_path / "registry.json"
    registry_file.write_text(json.dumps(definitions))
    cache = tmp_path / "cache"
    DatasetRegistry(definitions, cache_root=cache).download("fixture", show_progress=False)
    target = paired[1]["relations"][0].left
    prefix = ["datasets", "--registry", str(registry_file), "--cache-root", str(cache),
              "pairs", "fixture", "--basis", "same_material", "--pair-id", "pair"]
    args = ["--left-target", target.identifier, "--left-kind", target.kind,
            "--left-namespace", target.namespace, "--left-reference", target.reference.key]
    run(prefix + args + ["--left-unit", "TPM", "--left-scale", "linear", "--right-unit", "iBAQ"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(rows) == 1 and rows[0]["rna_abundance"][0]["unit"] == "TPM"
    assert rows[0]["protein_abundance"][0]["unit"] == "iBAQ"
    with pytest.raises(SystemExit) as exc:
        run(prefix + ["--left-target", target.identifier])
    assert exc.value.code == 2
    assert "Left target queries require" in capsys.readouterr().err
