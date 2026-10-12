"""Real old-format fixtures and compact-storage query/lineage invariants."""

import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, Measurement, QuantitativeDataset

FIXTURE = Path(__file__).parent / "data/storage"


@pytest.mark.parametrize("version", (1, 2, 3))
def test_released_storage_readers_and_roundtrip(version, tmp_path):
    old = QuantitativeDataset(FIXTURE / f"v{version}.sqlite")
    expected = json.loads((FIXTURE / f"v{version}.json").read_text())
    measurements = [Measurement.from_dict(m) for m in expected["measurements"]]
    assert old.manifest == DatasetManifest.from_dict(expected["manifest"])
    def key(measurement):
        return measurement.target.reference.key, measurement.measurement_id
    assert sorted(old.measurements(), key=key) == sorted(measurements, key=key)
    new = QuantitativeDataset.create(tmp_path / "v4.sqlite", old.manifest, samples=old.samples(),
                                    assays=old.assays(), measurements=measurements, mappings=old.mappings(),
                                    pairings=old.pairings(), relations=old.relations(),
                                    reference_mappings=old.reference_mappings())
    assert list(new.measurements()) == list(old.measurements())
    assert list(new.samples()) == list(old.samples()) and list(new.assays()) == list(old.assays())
    assert list(new.pairings()) == list(old.pairings()) and list(new.relations()) == list(old.relations())
    assert list(new.reference_mappings()) == list(old.reference_mappings())
    for reference in old.manifest.references:
        for quantity in (None, "rna_abundance", "protein_abundance"):
            for status in (None, "measured", "below_detection", "missing"):
                filters = dict(reference=reference, quantity=quantity, status=status)
                assert list(new.measurements(**filters)) == list(old.measurements(**filters))
        for assay in (None, "rna", "protein", "missing"):
            filters = dict(reference=reference, sample_id="sample", assay_id=assay)
            assert list(new.measurements(**filters)) == list(old.measurements(**filters))
        for value in old.measurements(reference=reference):
            assert list(new.measurements(target=value.target)) == list(old.measurements(target=value.target))
    with sqlite3.connect(new.path) as connection:
        assert connection.execute("PRAGMA user_version").fetchone()[0] == 4
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert connection.execute("SELECT count(*) FROM observation_input").fetchone()[0] == 4


def test_compact_storage_preserves_nondefault_measurement_fields(experiment, tmp_path):
    manifest, samples, assays, values = experiment
    original = values[0]
    values = [replace(original, value=0.0, scale="log", log_base=2, uncertainty=0.0,
                      uncertainty_type="standard_error", provenance={"source_asset": "input.tsv", "nested": [False, 0, None, {}, []]}),
              replace(original, measurement_id="censored", value=None, status="below_detection", detection_limit=0.0),
              replace(original, measurement_id="predicted", value=2.5, origin="predicted",
                      input_measurements=(original.measurement_id,), derivation={"method": "synthetic fixture"})]
    dataset = QuantitativeDataset.create(tmp_path / "fields.sqlite", manifest, samples=samples,
                                         assays=assays, measurements=values)
    assert {m.measurement_id: m.to_dict() for m in dataset.measurements()} == {m.measurement_id: m.to_dict() for m in values}
    with sqlite3.connect(dataset.path) as connection:
        records = [json.loads(row[0]) for row in connection.execute("SELECT record FROM observation")]
    assert all("dataset_id" not in r and "target" not in r and "measurement_id" not in r for r in records)
    assert any(r.get("uncertainty") == 0.0 for r in records)


def test_internal_foreign_keys_reject_cross_reference_links(experiment, tmp_path):
    manifest, samples, assays, values = experiment
    dataset = QuantitativeDataset.create(tmp_path / "scope.sqlite", manifest, samples=samples,
                                         assays=assays, measurements=values)
    with sqlite3.connect(dataset.path) as connection:
        connection.execute("PRAGMA foreign_keys=ON")
        target_i = connection.execute("SELECT i FROM target WHERE reference_i=2").fetchone()[0]
        with pytest.raises(sqlite3.IntegrityError):
            connection.execute("UPDATE observation SET target_i=? WHERE reference_i=1", (target_i,))
        assay_i = connection.execute("SELECT i FROM assay WHERE reference_i=2").fetchone()[0]
        with pytest.raises(sqlite3.IntegrityError):
            connection.execute("UPDATE observation SET assay_i=? WHERE reference_i=1", (assay_i,))
