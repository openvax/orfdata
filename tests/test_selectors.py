"""Native series selection preserves values, lineage and target resolution."""

import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry, Measurement, QuantitativeDataset
from orfdata.cli import run
from orfdata.sources.pfizer import import_protein


@pytest.fixture(params=("human", "mouse", "dog"))
def native(request, tmp_path):
    root = Path(__file__).parent / "data/pfizer" / request.param
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    return import_protein({a.asset_id: root / a.asset_id for a in manifest.sources}, tmp_path / "native.sqlite", manifest)


def test_native_unit_values_and_source_identity(native):
    original = list(native.measurements(quantity="protein_abundance"))
    assert len(original) == 10
    selected = {}
    for unit in ("iBAQ", "rIBAQ"):
        values = list(native.measurements(quantity="protein_abundance", unit=unit, scale="linear"))
        assert len(values) == 5 and values == [m for m in original if m.unit == unit]
        assert any(m.value == 0 and m.status == "measured" for m in values)
        assert any(m.value is None and m.status == "missing" for m in values)
        selected[unit] = {(m.sample_id, m.assay_id, m.provenance["line"]) for m in values}
    assert selected["iBAQ"] == selected["rIBAQ"]  # Same five source rows, two series.
    for unit in ("copies/cell", "ibaq", "iBAQ' OR 1=1 --"):
        assert list(native.measurements(unit=unit)) == []
    assert list(native.measurements(unit="iBAQ", scale="log")) == []
    mapping = next(m for m in native.reference_mappings() if m.targets)
    target = mapping.targets[0]
    assert list(native.mapped_measurements(target, unit="iBAQ")) == []
    views = list(native.mapped_measurements(target, mapping_statuses=["unresolved"], unit="iBAQ", scale="linear"))
    assert views and all(v.measurement.unit == "iBAQ" and v.attribution == "group_membership_context" for v in views)


def test_cli_native_and_mapped(native, tmp_path, capsys):
    from test_measurements import definition

    definitions = definition(native)
    path = tmp_path / "registry.json"
    path.write_text(json.dumps(definitions))
    cache = tmp_path / "cache"
    DatasetRegistry(definitions, cache_root=cache).download(native.manifest.dataset_id, show_progress=False)
    prefix = ["datasets", "--registry", str(path), "--cache-root", str(cache), "query", native.manifest.dataset_id]
    run(prefix + ["--unit", "iBAQ", "--scale", "linear"] + ["--format", "jsonl"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(rows) == 5 and all(r["unit"] == "iBAQ" for r in rows)
    target = next(m.targets[0] for m in native.reference_mappings() if m.targets)
    run(prefix + ["--mapped", "--mapping-status", "unresolved", "--reference", target.reference.key,
                  "--target", target.identifier, "--kind", target.kind, "--namespace", target.namespace,
                  "--unit", "rIBAQ", "--scale", "linear"] + ["--format", "jsonl"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert rows and all(r["measurement"]["unit"] == "rIBAQ" for r in rows)
    run(prefix + ["--unit", "copies/cell"] + ["--format", "jsonl"])
    assert capsys.readouterr().out == ""


def test_parquet_uses_the_same_selection(native, tmp_path):
    parquet = pytest.importorskip("pyarrow.parquet")
    output = tmp_path / "selected.parquet"
    native.export_parquet(output, unit="iBAQ", scale="linear")
    rows = parquet.read_table(output).to_pylist()
    assert len(rows) == 5 and all(r["unit"] == "iBAQ" for r in rows)
    assert [Measurement.from_dict(json.loads(r["measurement_json"])) for r in rows] == list(native.measurements(unit="iBAQ"))
    native.export_parquet(output, unit="unavailable")
    assert parquet.read_table(output).num_rows == 0


@pytest.mark.parametrize("version", (1, 2, 3))
def test_legacy_exact_series_filters(version):
    dataset = QuantitativeDataset(Path(__file__).parent / f"data/storage/v{version}.sqlite")
    original = list(dataset.measurements())
    for unit in ("TPM", "iBAQ", "fixture-units"):
        assert list(dataset.measurements(unit=unit, scale="linear")) == [m for m in original if m.unit == unit]
    assert list(dataset.measurements(scale="log")) == []
    assert list(dataset.measurements(log_base=2)) == []


@pytest.fixture
def transformed(experiment, tmp_path):
    manifest, samples, assays, values = experiment
    value = values[0]
    observations = [replace(value, measurement_id="linear", value=0),
                    replace(value, measurement_id="log2", value=-2, scale="log", log_base=2),
                    replace(value, measurement_id="log10", value=-1, scale="log", log_base=10)]
    return QuantitativeDataset.create(tmp_path / "log.sqlite", manifest, samples=samples, assays=assays, measurements=observations)


def test_log_base_selection_never_transforms(transformed):
    assert [m.value for m in transformed.measurements(unit="TPM", scale="linear")] == [0]
    assert {m.value for m in transformed.measurements(scale="log")} == {-1, -2}
    assert [m.value for m in transformed.measurements(unit="TPM", log_base=2)] == [-2]
    assert [m.value for m in transformed.measurements(scale="log", log_base=10)] == [-1]
    assert list(transformed.measurements(log_base=3)) == []


def test_cli_log_base_and_contradictory_selection(transformed, tmp_path, capsys):
    from test_measurements import definition

    definitions = definition(transformed)
    path = tmp_path / "registry.json"
    path.write_text(json.dumps(definitions))
    cache = tmp_path / "cache"
    DatasetRegistry(definitions, cache_root=cache).download("fixture", show_progress=False)
    prefix = ["datasets", "--registry", str(path), "--cache-root", str(cache), "query", "fixture"]
    run(prefix + ["--unit", "TPM", "--scale", "log", "--log-base", "2"] + ["--format", "jsonl"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(rows) == 1 and rows[0]["value"] == -2 and rows[0]["log_base"] == 2
    with pytest.raises(SystemExit):
        run(prefix + ["--scale", "linear", "--log-base", "2"] + ["--format", "jsonl"])
    assert "cannot have a logarithm base" in capsys.readouterr().err


@pytest.mark.parametrize("filters", [
    {"unit": ""}, {"unit": "  "}, {"unit": 1}, {"scale": "log2"}, {"log_base": 0},
    {"log_base": 1}, {"log_base": -2}, {"log_base": True}, {"log_base": float("inf")},
    {"log_base": float("nan")}, {"scale": "linear", "log_base": 2},
])
def test_invalid_selectors(transformed, filters):
    with pytest.raises((ValueError, TypeError)):
        list(transformed.measurements(**filters))
    target = next(transformed.measurements()).target
    with pytest.raises((ValueError, TypeError)):
        list(transformed.mapped_measurements(target, **filters))
