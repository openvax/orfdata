"""End-to-end release auditing checks values and indexes against source rows."""

import json
import sqlite3
from dataclasses import replace

import pytest
from test_pfizer_evidence import source as source
from test_pfizer_native import inputs as inputs

from orfdata.sources.pfizer_evidence import EvidenceStore
from orfdata.sources.pfizer_native import import_native
from tools import verify_dog_dataset


@pytest.fixture
def assembled(inputs, source, tmp_path, monkeypatch):
    paths, manifest, _ = inputs
    dataset = import_native(paths, tmp_path / "dataset.sqlite", manifest)
    plan = tmp_path / "audit-fixture-plan.json"
    plan.write_text(json.dumps(EvidenceStore(paths["evidence.sqlite"]).metadata["source"]["selection"]["plan"]))
    monkeypatch.setattr(verify_dog_dataset, "PLAN_PATH", plan)
    return dataset, source[0], paths


def check(assembled, **options):
    dataset, source_path, paths = assembled
    return verify_dog_dataset.verify(dataset.path, source_path, paths["search.fasta"], paths["protein.tsv"],
                                     lineage_path=paths["lineage.json"], **options)


def test_all_native_fields_indexes_and_public_measurements_are_checked(assembled):
    receipt = check(assembled)
    assert receipt["state"] == "verified"
    assert receipt["native_counts"]["TargetPsms"] == 3
    assert receipt["native_counts"]["WorkflowInputFiles"] == 97
    assert receipt["links"] == 3 and receipt["unresolved_links"] == 0
    assert receipt["projection_vectors"] == 12 and receipt["logical_cells"] == 429
    assert receipt["sampled_quantitative_records"] == 6 and receipt["sampled_common_measurements"] == 429
    assert receipt["independently_decoded_native_cells"] == 429
    assert receipt["search_proteins_compared"] == 10 and receipt["summary_native_matches"] == 1
    assert receipt["native_group_size_histogram"] == {1: 3}
    assert receipt["native_group_count_mismatches"] == []
    assert receipt["independently_searched_peptides"] == 6
    assert check(assembled, sample_records=1)["sampled_quantitative_records"] == 2


def test_independent_numeric_decoder_detects_a_shared_converter_error(assembled, monkeypatch):
    converter = verify_dog_dataset.NativeQuantification
    original = converter.peptide_abundances

    def shifted(self, *args, **kwargs):
        for measurement in original(self, *args, **kwargs):
            yield replace(measurement, value=measurement.value + 1) if measurement.value is not None else measurement

    monkeypatch.setattr(converter, "peptide_abundances", shifted)
    with pytest.raises(ValueError, match="Independent native numeric bytes"):
        check(assembled)


@pytest.mark.parametrize("failure,reason", [
    ("scalar", "typed fields"), ("vector", "typed fields"), ("source_key", "identity/sequence"),
    ("link", "Indexed membership"), ("projection", "quantitative target"),
    ("search_protein", "search protein"), ("candidate", "Candidate status"),
    ("summary", "Summary representative/native"), ("protein_index", "protein/search identity"),
    ("sample", "samples/assays differ"), ("lineage", "Retained lineage differs"),
])
def test_data_and_index_corruption_is_detected(assembled, failure, reason):
    dataset, _, _ = assembled
    with sqlite3.connect(dataset.path) as db:
        peptide = db.execute("SELECT r.i FROM evidence_record r JOIN evidence_table t ON t.i=r.table_i "
                             "WHERE t.name='TargetPeptideGroups' ORDER BY r.source_rowid LIMIT 1").fetchone()[0]
        if failure == "scalar":
            db.execute("UPDATE evidence_record SET attributes=json_set(attributes,'$.Sequence','CHANGED') WHERE i=?", (peptide,))
        elif failure == "vector":
            peptide, original = db.execute("SELECT v.record_i,v.values_data FROM evidence_vector v "
                                           "JOIN evidence_record r ON r.i=v.record_i JOIN evidence_table t ON t.i=r.table_i "
                                           "WHERE t.name='TargetPsms' AND v.column_name='Abundances' AND length(v.values_data)>0 LIMIT 1").fetchone()
            db.execute("UPDATE evidence_vector SET values_data=? WHERE record_i=? AND column_name='Abundances'",
                       (bytes([original[0] ^ 1]) + original[1:], peptide))
        elif failure == "source_key":
            db.execute("UPDATE evidence_record SET source_key='[\"9999\"]' WHERE i=?", (peptide,))
        elif failure == "link":
            db.execute("UPDATE evidence_link SET left_i=NULL WHERE record_i=(SELECT min(record_i) FROM evidence_link)")
        elif failure == "projection":
            db.execute("UPDATE peptide_projection SET unit='changed unit' WHERE record_i=?", (peptide,))
        elif failure == "search_protein":
            db.execute("UPDATE search_protein SET record=json_set(record,'$.header','changed') WHERE accession=(SELECT min(accession) FROM search_protein)")
        elif failure == "candidate":
            sequence, = db.execute("SELECT sequence FROM search_peptide_hit LIMIT 1").fetchone()
            db.execute("DELETE FROM search_peptide_hit WHERE sequence=?", (sequence,))
        elif failure == "summary":
            db.execute("DELETE FROM summary_native")
        elif failure == "sample":
            db.execute("UPDATE sample SET record=json_set(record,'$.tissue','changed') WHERE id=(SELECT min(id) FROM sample)")
        elif failure == "lineage":
            db.execute("UPDATE metadata SET value=json_set(value,'$.lineage.channels[0].tissue','changed') WHERE name='measurement_storage'")
        else:
            db.execute("UPDATE native_protein_search SET status='sequence_mismatch'")
    with pytest.raises(ValueError, match=reason):
        check(assembled)


def test_original_source_and_complete_selection_pins_are_required(assembled, tmp_path, monkeypatch):
    empty = tmp_path / "wrong-plan.json"
    empty.write_text('{"tables": []}')
    monkeypatch.setattr(verify_dog_dataset, "PLAN_PATH", empty)
    with pytest.raises(ValueError, match="complete pinned"):
        check(assembled)
    _, source_path, _ = assembled
    with sqlite3.connect(source_path) as db:
        db.execute("UPDATE TargetPsms SET Sequence='CHANGED'")
    with pytest.raises(ValueError, match="native file differs"):
        check(assembled)
