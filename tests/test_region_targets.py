"""Source windows stay distinct from complete ORFs across storage and queries."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest
from test_measurements import definition, write_dataset

from orfdata import (
    Assay,
    DatasetManifest,
    Measurement,
    OrthologyMapping,
    QuantitativeDataset,
    ReferenceContext,
    Sample,
    ScopedTarget,
    SourceAsset,
    SourceRNA,
    TranscriptORF,
)
from orfdata.cli import run
from orfdata.pairing import AssayPair, MatchedMeasurements, TargetRelation


def test_native_regions_roundtrip_and_cli(tmp_path, capsys):
    path = Path(__file__).parent / "data/regions/ingolia.json"
    source = json.loads(path.read_text())
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    reference = ReferenceContext(4932, "yeast", None, "SGD", "2008-06-22",
                                 source_sha256=sha)
    manifest = DatasetManifest(
        "regions", "fixture", "Native Ingolia quantification slices", (reference,),
        (SourceAsset("slice.json", path.as_uri(), sha, path.stat().st_size,
                     "GEO factual source slice; Ingolia et al. 2009", redistributable=True),),
        ("https://doi.org/10.1126/science.1168978",), ("rna_seq", "ribo_seq"),
        dict(source_features=4, source_assays=2, genomic_mapping="unresolved"),
        dict(scope="native source-region round-trip; not a full source importer"))
    samples, assays, values = [], [], []
    for entry in source:
        sid = entry["sample"]
        rna = sid == "GSM346117"
        samples.append(Sample(sid, reference, attributes=dict(strain="BY4741")))
        assays.append(Assay(sid, sid, "rna_seq" if rna else "ribo_seq", reference))
        for record in entry["rows"]:
            row = record["native"]
            target = ScopedTarget(reference, "region", row["yorf"], "Ingolia2009:feature")
            for column, unit in (("count", "counts"), ("norm", "RPKM"), ("dens", "reads/kb")):
                values.append(Measurement(
                    f"{sid}:{record['line']}:{column}", target, sid, sid,
                    "rna_abundance" if rna else "ribosome_occupancy", float(row[column]),
                    unit, manifest.dataset_id, manifest.version,
                    provenance=dict(source_asset="slice.json", upstream=entry["source"],
                                    row=record["line"], column=column, native_value=row[column],
                                    feature_length_nt=int(row["len"]), total_cds_reads=int(row["total"]))))
    dataset = QuantitativeDataset.create(tmp_path / "source.sqlite", manifest,
                                         samples=samples, assays=assays, measurements=values)
    assert dataset.counts["measurement"] == 24 and dataset.counts["target"] == 4
    assert {m.measurement_id: m for m in dataset.measurements()} == {m.measurement_id: m for m in values}
    target = ScopedTarget(reference, "region", "YAL005C_5'utr", "Ingolia2009:feature")
    expected = list(dataset.measurements(target=target, unit="RPKM"))
    assert len(expected) == 2 and {m.quantity for m in expected} == {"rna_abundance", "ribosome_occupancy"}
    assert not list(dataset.measurements(target=replace(target, kind="orf")))
    assert all(Measurement.from_dict(m.to_dict()) == m for m in expected)
    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps(definition(dataset)))
    args = ["datasets", "--registry", str(registry), "--cache-root", str(tmp_path / "cache")]
    run(args + ["download", "regions", "--quiet"])
    capsys.readouterr()
    run(args + ["query", "regions", "--reference", reference.key, "--kind", "region",
                "--namespace", target.namespace, "--target", target.identifier, "--unit", "RPKM"])
    assert [Measurement.from_dict(json.loads(line)) for line in capsys.readouterr().out.splitlines()] == expected


def test_region_ids_do_not_cross_species_or_reference(tmp_path, experiment):
    observations = tuple(replace(m, target=replace(m.target, kind="region")) for m in experiment[3])
    dataset = write_dataset(tmp_path / "regions.sqlite", experiment, measurements=observations)
    for value in observations:
        assert list(dataset.measurements(target=value.target)) == [value]
        wrong = replace(value.target, reference=replace(value.target.reference, annotation_version="111"))
        with pytest.raises(ValueError, match="not present"):
            list(dataset.measurements(target=wrong))


def test_same_region_ratio_does_not_become_an_orf_estimate(experiment):
    rna = experiment[3][-1]
    rna = replace(rna, target=replace(rna.target, kind="region"), value=2, unit="RPKM")
    rpf = replace(rna, measurement_id="rpf", assay_id="rpf", quantity="ribosome_occupancy", value=6)
    reference = rna.target.reference
    pair = AssayPair("pair", reference, "fixture", "1", "rna", "rpf", "same_material",
                     dict(scope="synthetic assay pairing"))
    relation = TargetRelation("region", "fixture", "1", rna.target, rpf.target, "same_target",
                              "one_to_one", dict(scope="same synthetic region"))
    ratio = MatchedMeasurements(pair, relation, (rna,), (rpf,)).ratio(
        rpf, rna, measurement_id="ratio", method="explicit source-window density ratio",
        normalization=dict(scope="same synthetic quantification window"))
    assert ratio.value == 3 and ratio.target.kind == "region"
    assert ratio.quantity == "translation_efficiency" and ratio.origin == "derived"
    with pytest.raises(ValueError, match="gene or transcript"):
        ScopedTarget.from_pyensembl(object(), kind="region", namespace="test")
    with pytest.raises(ValueError, match="canonical genes or independently supported ORFs"):
        OrthologyMapping("regions", rna.target,
                         replace(rna.target, reference=experiment[0].references[0]),
                         experiment[0].sources[0], "test", "one_to_one", {"score": 1}, {"source": "synthetic"})


def test_region_to_partial_candidates_keeps_measurements_and_uncertainty(experiment, tmp_path):
    manifest, samples, assays, observations = experiment
    reference = manifest.references[0]
    region = ScopedTarget(reference, "region", "window", "synthetic:regions")
    rna = SourceRNA("rna", "ATGGCTAAATAG", reference, "synthetic")
    candidates = [TranscriptORF.at(rna, start, end=9, start_complete=None, stop_complete=False,
                                   orf_id=f"candidate-{start}", source="synthetic") for start in (0, 3)]
    targets = [ScopedTarget.from_orf(orf, namespace="synthetic:candidates") for orf in candidates]
    assert [target.identifier for target in targets] == [orf.orf_id for orf in candidates]
    assert all(orf.record is None for orf in candidates)
    left = replace(observations[0], target=region, value=2, unit="RPKM")
    right = [replace(left, measurement_id=f"rpf-{i}", assay_id="rpf", target=target,
                     quantity="ribosome_occupancy", value=6) for i, target in enumerate(targets)]
    relations = [TargetRelation(
        f"context-{i}", manifest.dataset_id, manifest.version, region, target,
        "region_orf_context", "one_to_many",
        dict(source_asset="input.tsv", basis="source-reported candidate within measured window",
             source_window=[0, 12], coding_interval=orf.cdna_interval,
             start_complete=orf.start_complete, stop_complete=orf.stop_complete,
             frame=orf.coding_phase, translation_status="candidate"))
        for i, (target, orf) in enumerate(zip(targets, candidates))]
    pair = AssayPair("pair", reference, manifest.dataset_id, manifest.version, "rna", "rpf",
                     "same_material", dict(scope="synthetic source pairing"), "concurrent")
    dataset = QuantitativeDataset.create(
        tmp_path / "candidates.sqlite", manifest, samples=samples,
        assays=(*assays, Assay("rpf", "sample", "ribo_seq", reference)),
        measurements=(left, *right), pairings=(pair,), relations=relations)
    assert list(dataset.measurements(target=region)) == [left]
    assert list(dataset.relations(left_target=region)) == relations
    matches = list(dataset.matched_measurements(allowed_bases=["same_material"]))
    assert len(matches) == 2
    for match in matches:
        assert match.left == (left,)
        assert len(match.right) == 1 and match.right[0].target.kind == "orf"
        restored = MatchedMeasurements.from_dict(json.loads(json.dumps(match.to_dict())))
        assert restored == match and restored.relation.evidence["start_complete"] is None
        # Even a source-declared single candidate cannot receive a window's ratio.
        for cardinality in ("one_to_many", "one_to_one", "unknown"):
            candidate = replace(match, relation=replace(match.relation, cardinality=cardinality))
            with pytest.raises(ValueError, match="region context"):
                candidate.ratio(candidate.right[0], left, measurement_id="invalid-ratio",
                                method="test", normalization=dict(scope="synthetic"))
    with pytest.raises(ValueError, match="source region"):
        replace(relations[0], left=targets[0])
    with pytest.raises(ValueError, match="source region"):
        replace(relations[0], right=region)
    with pytest.raises(ValueError, match="require mapping evidence"):
        replace(relations[0], evidence={})
