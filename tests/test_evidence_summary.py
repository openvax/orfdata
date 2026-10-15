"""Aggregate counts retain source identity and cannot become resolved evidence."""

import csv
import hashlib
import io
import json
from dataclasses import replace
from pathlib import Path

import pytest
from test_reference_orf import catalog, model, mouse_transcript  # noqa: F401

from orfdata import (
    EvidenceSummary,
    ORFCatalog,
    ORFEvidence,
    ORFRecord,
    ReferenceContext,
    ScopedTarget,
    load_catalog,
    load_uorf_catalog,
)
from orfdata.cli import run

EXAMPLES = json.loads((Path(__file__).parent / "data/aggregate_evidence.json").read_text())["examples"]


def summary(**kwargs):
    """Synthetic evidence for the synthetic mouse binding fixture only."""
    context = ReferenceContext(10090, "mouse", "GRCm39", "synthetic", "1")
    values = dict(source="synthetic", release="1", target=ScopedTarget(context, "protein", "native", "fixture"),
                  protein_sha256=hashlib.sha256(b"MKPG").hexdigest(), method="ms", metric="support count",
                  value=2, definition="Synthetic test count, not experimental data", scope="protein_sequence",
                  source_url="https://example.test/fixture", source_sha256="a" * 64)
    values.update(kwargs)
    return EvidenceSummary(**values)


@pytest.mark.parametrize("example", EXAMPLES, ids=lambda e: e["species"] + ":" + e["native"]["transcript accession"])
def test_native_multispecies_counts_round_trip_without_resolved_observations(example):
    summaries = tuple(EvidenceSummary.from_dict(row) for row in example["summaries"])
    ms, ribo = summaries
    assert ms.value == int(example["native"]["MS score"])
    assert ribo.value == int(example["native"]["TE Translation Event score"])
    assert ms.protein_sha256 == hashlib.sha256(example["protein_sequence"].encode()).hexdigest()
    assert ms.target.reference.taxon_id == {"human": 9606, "mouse": 10090, "yeast": 4932}[example["species"]]
    if example["species"] == "yeast":
        assert ms.target.reference.strain_taxon_id == 559292
    for s in summaries:
        assert EvidenceSummary.from_dict(json.loads(json.dumps(s.to_dict()))) == s
    evidence = ORFEvidence(summaries=summaries)
    assert evidence.peptides == evidence.ribosome_studies == ()
    assert evidence.has_evidence("ms") and not evidence.has_evidence("proteomics")
    assert evidence.has_evidence("ribosome") == (ribo.value > 0)
    for flag in ("require_resolved", "require_unique", "require_reviewed"):
        assert not evidence.has_evidence("ms", **{flag: True})
    assert not evidence.has_evidence("ribosome", require_resolved=True)


@pytest.mark.parametrize("value", (None, 0, 2))
def test_unknown_zero_and_positive_counts_remain_distinct(value):
    s = summary(value=value)
    assert s.to_dict()["value"] is value
    assert EvidenceSummary.from_dict(s.to_dict()).value == value
    assert ORFEvidence(summaries=(s,)).has_evidence("ms") == (value == 2)


@pytest.mark.parametrize("value", (-1, 0.5, 2.0, True, "2", float("nan"), float("inf")))
def test_arbitrary_scores_and_invalid_counts_are_rejected(value):
    with pytest.raises(ValueError, match="support count"):
        summary(value=value)


def test_shared_native_entity_deduplicates_across_transcript_assignments():
    first, second = EXAMPLES[:2]
    assert first["native"]["transcript accession"] != second["native"]["transcript accession"]
    a, b = (EvidenceSummary.from_dict(x["summaries"][0]) for x in (first, second))
    assert a == b and a.key == b.key
    assert ORFEvidence(summaries=(a, b)).summaries == (a,)
    with pytest.raises(ValueError, match="Conflicting"):
        ORFEvidence(summaries=(a, replace(b, value=b.value + 1)))


def test_reference_aliases_deduplicate_but_species_and_native_entities_do_not():
    s = summary()
    alias = replace(s, target=replace(s.target, reference=replace(s.target.reference, species="mus_musculus")))
    assert s.key == alias.key and ORFEvidence(summaries=(s, alias)).summaries == (s,)
    another = replace(s, target=replace(s.target, identifier="another-native-entity"))
    human = replace(s, target=replace(s.target, reference=replace(s.target.reference, taxon_id=9606, species="human")))
    assert len(ORFEvidence(summaries=(s, another, human)).summaries) == 3


def test_catalogue_preserves_one_summary_per_entity_and_rejects_conflicts(tmp_path):
    s = summary()
    a = ORFRecord.from_model(model(), ORFEvidence(summaries=(s,)))
    b = replace(a, orf_id="other", transcript_id="other-transcript")
    data = ORFCatalog([a, b], catalog(a).metadata)
    assert data.evidence_summaries == (s,)
    assert data.query(evidence="ms") == (a, b)
    assert data.query(evidence="ms", require_resolved=True) == ()
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(data.to_dict()))
    restored = load_catalog(path)
    assert restored.records == data.records and restored.evidence_summaries == (s,)
    assert ORFRecord.from_model(a.model, a.evidence) == a
    conflicting = replace(b, summaries=(replace(s, definition="Another counting scope"),))
    with pytest.raises(ValueError, match="Conflicting"):
        ORFCatalog([a, conflicting], data.metadata)


def test_summary_cannot_attach_to_another_species_or_protein_or_orf():
    s = summary()
    with pytest.raises(ValueError, match="species"):
        ORFRecord.from_model(model(taxon_id=9606, species="human"), ORFEvidence(summaries=(s,)))
    with pytest.raises(ValueError, match="protein sequences"):
        ORFRecord.from_model(model(protein_sequence="MKPA"), ORFEvidence(summaries=(s,)))
    native = replace(s, scope="orf", target=replace(s.target, kind="orf", identifier="test-orf"))
    ORFRecord.from_model(model(), ORFEvidence(summaries=(native,)))
    with pytest.raises(ValueError, match="ORF identities"):
        ORFRecord.from_model(model(orf_id="another"), ORFEvidence(summaries=(native,)))
    wrong_ref = replace(native, target=replace(native.target, reference=replace(native.target.reference, genome_build="GRCm38")))
    with pytest.raises(ValueError, match="assemblies"):
        ORFRecord.from_model(model(), ORFEvidence(summaries=(wrong_ref,)))
    for changes in (dict(genome_build="GRCm38"), dict(strain="C57BL/6J")):
        other = replace(s, target=replace(s.target, reference=replace(s.target.reference, **changes)))
        with pytest.raises(ValueError, match="assemblies|strain contexts"):
            ORFRecord.from_model(model(), ORFEvidence(summaries=(other,)))


def test_bound_reference_and_base_annotation_never_invent_variant_peptides(mouse_transcript):  # noqa: F811
    s = summary()
    record = ORFRecord.from_model(model(), ORFEvidence(summaries=(s,)))
    bound = record.bind(mouse_transcript)
    assert bound.transcript is mouse_transcript and bound.evidence.summaries == (s,)
    assert bound.to_dict()["evidence"]["summaries"] == (s.to_dict(),)
    assert bound.evidence.peptides == ()
    annotation = record.map_genomic_position(205, genome_build="GRCm39")
    assert annotation["reference_peptides"] == []
    assert not annotation["reference_ms_support"]
    assert annotation["mutant_translation_evidence"] == "not_assessed"


def test_resolved_filters_preserve_old_catalogue_behavior_and_exports():
    data = load_catalog()
    assert data.evidence_summaries == ()
    for method in ("ms", "ribosome", "proteomics", "immunopeptidomics", "ms_unclassified"):
        assert data.query(evidence=method) == data.query(evidence=method, require_resolved=True)
    record = data.get("c7riboseqorf56")
    assert "summaries" not in record.to_dict() and "summaries" not in record.evidence.to_dict()
    assert record.has_evidence("immunopeptidomics", require_resolved=True, require_unique=True)
    with pytest.raises(ValueError, match="requires an evidence method"):
        data.query(require_resolved=True)


def test_cli_resolution_and_native_summary_export(tmp_path, capsys):
    s = summary()
    record = ORFRecord.from_model(model(), ORFEvidence(summaries=(s,)))
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog(record).to_dict()))
    run(["--catalog", str(path), "--evidence", "ms"])
    result = json.loads(capsys.readouterr().out)
    assert result["records"][0]["summaries"] == [s.to_dict()]
    run(["--catalog", str(path), "--evidence", "ms", "--require-resolved"])
    assert json.loads(capsys.readouterr().out)["records"] == []


def test_tsv_exports_native_summaries_including_unknown_zero_and_empty(tmp_path, capsys):
    # Regression for #58: DictWriter previously rejected the summaries field.
    summaries = (summary(), summary(metric="unknown count", value=None),
                 summary(method="ribosome", metric="study count", value=0))
    supported = ORFRecord.from_model(model(), ORFEvidence(summaries=summaries))
    plain = ORFRecord.from_model(model(orf_id="plain"))
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(ORFCatalog([supported, plain], catalog(plain).metadata).to_dict()))
    run(["--catalog", str(path), "--format", "tsv"])
    rows = csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t")
    exported = {row["orf_id"]: json.loads(row["summaries_json"]) for row in rows}
    assert exported == {model().orf_id: [s.to_dict() for s in summaries], "plain": []}
    restored = [EvidenceSummary.from_dict(s) for s in exported[model().orf_id]]
    assert restored == list(summaries) and [s.value for s in restored] == [2, None, 0]
    assert {s.key for s in restored} == {s.key for s in summaries}


def test_legacy_upstream_tsv_keeps_columns_and_rejects_aggregate_evidence(tmp_path, capsys):
    run(["--gene", "TPST1", "--format", "tsv"], catalog_loader=lambda path: load_uorf_catalog())
    header = capsys.readouterr().out.splitlines()[0].split("\t")
    assert header[-1] == "initiation_context_json" and "summaries_json" not in header
    record = ORFRecord.from_model(model(orf_type="uORF"), ORFEvidence(summaries=(summary(),)))
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog(record).to_dict()))
    with pytest.raises(SystemExit):
        run(["--catalog", str(path), "--format", "tsv"], catalog_loader=load_uorf_catalog)
    assert "use --format json" in capsys.readouterr().err


def test_method_scope_source_pins_and_nonhuman_immunopeptidomics():
    s = summary(method="immunopeptidomics")
    assert ORFEvidence(summaries=(s,)).has_evidence("immunopeptidomics")
    assert ORFEvidence(summaries=(s,)).has_evidence("ms")
    assert not ORFEvidence(summaries=(s,)).has_evidence("proteomics")
    for changes in (dict(source=""), dict(release=""), dict(definition=""),
                    dict(source_sha256="unknown"), dict(protein_sha256="A" * 64),
                    dict(source_url="local.txt"), dict(scope="gene"), dict(method="confidence")):
        with pytest.raises(ValueError):
            replace(s, **changes)
    ribo = ORFEvidence(ribosome_studies=("named-study",), summaries=(summary(method="ribosome"),))
    assert ribo.has_evidence("ribosome", require_resolved=True)
    assert not ribo.has_evidence("ms")
    with pytest.raises(ValueError, match="MS filters"):
        ribo.has_evidence("ribosome", require_unique=True)
