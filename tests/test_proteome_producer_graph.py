"""Portable producer records preserve evidence and no-sequence outcomes offline."""

from dataclasses import replace

import pytest

from orfdata import (
    ProducerInfo,
    ProducerLink,
    ProducerRecord,
    Proteome,
    ProteomeEntry,
    ReferenceContext,
)


@pytest.fixture
def reference():
    return ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")


def test_offline_graph_has_alternative_translations_shared_evidence_and_outcomes(reference, tmp_path):
    baseline = Proteome(reference, [ProteomeEntry("ENSP1", "MAA", "annotated", 9606,
                                                 "Ensembl:115")])
    dna = ProteomeEntry.from_producer("MTA", taxon_id=9606, producer=ProducerInfo(
        "reference_cds", "dna_prediction", "varcode", "V1", sample_id="sample-A",
        reference=reference))
    rna = ProteomeEntry.from_producer("MTA", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "ASM1", sample_id="sample-A",
        reference=reference, start_complete=False, stop_complete=None))
    version = "fixture.v1"
    translation_a = ProducerRecord("isovar", version, "translation", "T1",
                                   sample_id="sample-A", payload={"cdna": "AUGACA"})
    translation_b = ProducerRecord("isovar", version, "translation", "T2",
                                   sample_id="sample-A", payload={"cdna": "AUGACG"})
    evidence = ProducerRecord("isovar", version, "evidence", "reads-1",
                              sample_id="sample-A", support="rna_sequence",
                              scope={"sample": "sample-A", "read_set": "R1"},
                              payload={"read_count": 7})
    partial = ProducerRecord("isovar", version, "event", "E-partial",
                             sample_id="sample-A", status="incomplete", complete=False,
                             payload={"start_known": False, "stop_known": None})
    zero = ProducerRecord("varcode", "fixture.v1", "event", "E-zero",
                          sample_id="sample-A", status="zero_product", complete=True,
                          payload={"reason": "no_coding_product"})
    capped = ProducerRecord("isovar", version, "event", "E-capped",
                            sample_id="sample-A", status="capped", complete=False,
                            limits={"max_hypotheses": 4}, payload={"remaining": "unknown"})
    links = (
        ProducerLink(translation_a.key, rna.accession, "entry", "translates_to"),
        ProducerLink(translation_b.key, rna.accession, "entry", "translates_to"),
        ProducerLink(evidence.key, translation_a.key, "record", "supports_translation"),
        ProducerLink(evidence.key, translation_b.key, "record", "supports_translation"),
        ProducerLink(partial.key, rna.accession, "entry", "produces"),
    )
    records = (translation_a, translation_b, evidence, partial, zero, capped)
    external = Proteome.from_entries((dna, rna), reference=reference,
                                     records=records, links=links)
    joint = baseline.merge(external)
    assert len(joint.by_sequence()["MTA"]) == 2
    assert sum(record.kind == "evidence" for record in joint.records) == 1
    assert len([link for link in joint.links if link.source == evidence.key]) == 2
    assert not any(link.source in (zero.key, capped.key) and link.target_kind == "entry"
                   for link in joint.links)
    assert joint.merge(external).identity == joint.identity
    path = tmp_path / "joint.json.gz"
    joint.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == joint.identity
    assert loaded.records == joint.records and loaded.links == joint.links
    assert loaded.records[-1].limits["max_hypotheses"] == 4
    assert loaded.records[2].support == "rna_sequence"
    assert loaded.records[2].scope["read_set"] == "R1"
    changed_event = replace(capped, limits={"max_hypotheses": 5})
    changed_evidence = replace(evidence, scope={"sample": "sample-A", "read_set": "R2"})
    assert joint.identity != baseline.merge(Proteome.from_entries(
        (dna, rna), reference=reference,
        records=records[:-1] + (changed_event,), links=links)).identity
    assert joint.identity != baseline.merge(Proteome.from_entries(
        (dna, rna), reference=reference,
        records=(translation_a, translation_b, changed_evidence, partial, zero, capped),
        links=links)).identity


def test_graph_rejects_conflicts_and_missing_links(reference):
    entry = ProteomeEntry.from_producer("MAA", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P1", reference=reference))
    record = ProducerRecord("future.tool", "schema.9", "translation", "T1",
                            payload={"nested": {"unrecognized": [1, 2]}})
    collection = Proteome.from_entries([entry], reference=reference, records=[record],
                                       links=[ProducerLink(record.key, entry.accession,
                                                           "entry", "translates_to")])
    assert collection.records[0].payload["nested"]["unrecognized"] == (1, 2)
    assert collection.with_records([record], links=collection.links).identity == collection.identity
    with pytest.raises(ValueError, match="Conflicting producer record identity"):
        collection.with_records([replace(record, payload={"nested": {"unrecognized": [3]}})])
    with pytest.raises(ValueError, match="unknown target"):
        collection.with_records(links=[ProducerLink(record.key, "missing", "entry", "produces")])
    with pytest.raises(ValueError, match="explicit support"):
        ProducerRecord("isovar", "1", "evidence", "R1")
    with pytest.raises(ValueError, match="explicit outcome"):
        ProducerRecord("isovar", "1", "event", "E1")


def test_filter_and_follow_on_entry_addition_keep_valid_graph(reference):
    entry = ProteomeEntry.from_producer("MAA", taxon_id=9606, producer=ProducerInfo(
        "nc_orf", "rna_assembly", "isovar", "P1", reference=reference))
    record = ProducerRecord("isovar", "1", "translation", "T1")
    collection = Proteome.from_entries([entry], reference=reference, records=[record],
                                       links=[ProducerLink(record.key, entry.accession,
                                                           "entry", "translates_to")])
    assert collection.with_entries([]).links == collection.links
    filtered = collection.filter(lambda _: False)
    assert filtered.records == (record,) and filtered.links == ()
