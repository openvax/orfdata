"""Source windows and native normalization do not imply complete translated ORFs."""

import csv
import gzip
import io
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry, ScopedTarget, build_dataset
from orfdata.cli import run
from orfdata.sources.ingolia import import_ribo

ROOT = Path(__file__).parent / "data/ingolia"


def test_packaged_yeast_release_discovery_is_offline(monkeypatch):
    import requests

    from orfdata.sources.ingolia import ribo_manifest

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Discovery contacted network"))
    registry = DatasetRegistry()
    manifest = ribo_manifest()
    assert manifest.dataset_id == "ingolia-yeast-rna-ribo"
    assert manifest.coverage["observations"] == 329952
    assert manifest in registry.list(species="yeast", assay="ribo_seq")
    assert manifest in registry.list(species="yeast", assay="rna_seq")
    assert manifest not in registry.list(species="human")
    assert manifest not in registry.list(species="mouse")
    assert manifest not in registry.list(species="dog")
    assert len(manifest.sources) == 17 and len(manifest.references) == 1
    assert manifest.references[0].taxon_id == 4932
    assert not manifest.references[0].is_genomic


@pytest.fixture
def source():
    manifest = DatasetManifest.from_dict(json.loads((ROOT / "manifest.json").read_text()))
    return manifest, {a.asset_id: ROOT / a.asset_id for a in manifest.sources}


def build(source, tmp_path):
    manifest, inputs = source
    return build_dataset(manifest, lambda i, o: import_ribo(i, o, manifest), source_paths=inputs,
                         cache_root=tmp_path, download_options={"show_progress": False})


def test_all_reported_values_and_normalization_inputs_are_preserved(source, tmp_path):
    dataset = build(source, tmp_path)
    expected = json.loads((ROOT / "expected.json").read_text())
    values = {m.measurement_id: m for m in dataset.measurements()}
    assert len(values) == 192 and dataset.counts["target"] == 8
    for sid, sample in expected.items():
        for record in sample["rows"]:
            row, line = record["native"], record["fixture_line"]
            count_id = f"{sid}:{line}:count"
            for column, unit in (("count", "counts"), ("dens", "reads/kb"), ("norm", "RPKM")):
                value = values[f"{sid}:{line}:{column}"]
                assert value.value == float(row[column]) and value.status == "measured"
                assert value.unit == unit and value.target.kind == "region"
                assert value.target.identifier == row["yorf"] and value.target.reference == source[0].references[0]
                assert value.target.namespace == "Ingolia2009:feature"
                assert value.sample_id == value.assay_id == sid
                assert value.provenance["native_value"] == row[column]
                assert value.provenance["feature_length_nt"] == int(row["len"])
                assert value.provenance["total_cds_reads"] == int(row["total"])
                assert value.provenance["line"] == line and value.provenance["column"] == column
                if column == "count":
                    assert value.origin == "observed" and value.input_measurements == ()
                else:
                    assert value.origin == "derived" and value.input_measurements == (count_id,)
                    assert values[count_id].target == value.target
                    assert value.derivation["source_value_preserved"] is True
                    assert value.derivation["reported_decimal_places"] == 3
    assert sum(m.origin == "derived" for m in values.values()) == 128
    assert any(m.value == 0 for m in values.values())
    assert all(m.value > 0 for m in values.values() if m.sample_id == "GSM346118" and m.unit == "counts")
    assert {m.quantity for m in values.values()} == {"rna_abundance", "ribosome_occupancy"}
    assert {m.target.identifier for m in values.values()} >= {"HRA1", "Q0080", "15S_rRNA", "YAL001C_intron", "YAL005C_5'utr"}
    assert list(dataset.pairings()) == list(dataset.reference_mappings()) == []
    assert list(dataset.matched_measurements(allowed_bases=["same_material"])) == []
    assert not source[0].references[0].is_genomic
    with pytest.raises(ValueError, match="Incomplete source reference"):
        source[0].references[0].validate_genome(object())
    samples = {s.sample_id: s for s in dataset.samples()}
    assert samples["GSM346114"].attributes["source_title"] == "Footprints-rich-2"
    assert all(s.biological_unit_id is None and s.biological_replicate is None for s in samples.values())
    assert all(s.attributes["biological_strain"] == "BY4741" for s in samples.values())
    assert all(not any("contact" in key for key in s.attributes["source_geo"]) for s in samples.values())
    assert all(a.run_id is None and not a.calibration["absolute_molecule_counts"] for a in dataset.assays())
    with pytest.raises(ValueError, match="Reference context"):
        list(dataset.measurements(target=ScopedTarget(replace(source[0].references[0], annotation_version="current"),
                                                      "region", "HRA1", "Ingolia2009:feature")))


@pytest.mark.parametrize("change,error", [
    ("negative", "nonnegative integers"), ("fractional", "nonnegative integers"),
    ("nan", "normalization"), ("density", "normalization"), ("length", "lengths disagree"),
    ("denominator", "normalization denominator"), ("duplicate", "Duplicate"),
    ("missing", "inventory reconciliation"), ("unknown", "inventory reconciliation"),
    ("header", "columns"),
])
def test_changed_source_rows_fail_atomically(source, tmp_path, change, error):
    manifest, inputs = source
    asset = "GSM346114.quant.txt.gz"
    with gzip.open(inputs[asset], "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns, rows = list(reader.fieldnames), list(reader)
    if change in ("negative", "fractional"):
        rows[0]["count"] = "-1" if change == "negative" else "1.5"
    elif change in ("nan", "density"):
        rows[0]["dens"] = "NaN" if change == "nan" else "0.0"
    elif change == "length":
        rows[0]["len"] = str(int(rows[0]["len"]) + 1)
    elif change == "denominator":
        rows[0]["total"] = str(int(rows[0]["total"]) + 1)
    elif change == "duplicate":
        rows.insert(1, dict(rows[0]))
    elif change == "missing":
        rows.pop()
    elif change == "unknown":
        rows[0]["yorf"] = "unreported-feature"
    else:
        columns[1] = "wrong"
        for row in rows:
            row["wrong"] = row.pop("norm")
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    altered = tmp_path / "altered.gz"
    altered.write_bytes(gzip.compress(buffer.getvalue().encode(), mtime=0))
    output = tmp_path / "bad.sqlite"
    with pytest.raises(ValueError, match=error):
        import_ribo(dict(inputs, **{asset: altered}), output, manifest)
    assert not output.exists()


@pytest.mark.parametrize("change", ("title", "organism", "biosample", "reference"))
def test_source_identity_changes_fail(source, tmp_path, change):
    manifest, inputs = source
    if change == "reference":
        manifest = replace(manifest, references=(replace(manifest.references[0], annotation_version="current"),))
    else:
        asset = "GSM346111.metadata.txt"
        old, new = {"title": ("Footprints-rich-1", "Footprints-rich-2"),
                    "organism": ("Saccharomyces cerevisiae", "Homo sapiens"),
                    "biosample": ("SAMN02195614", "SAMN00000000")}[change]
        altered = tmp_path / "metadata.txt"
        altered.write_text(inputs[asset].read_text().replace(old, new))
        inputs = dict(inputs, **{asset: altered})
    with pytest.raises(ValueError, match="source"):
        import_ribo(inputs, tmp_path / "bad.sqlite", manifest)
    assert not (tmp_path / "bad.sqlite").exists()


def test_offline_registry_source_rebuild_and_region_cli(source, tmp_path, monkeypatch, capsys):
    import requests
    from test_measurements import definition

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Offline source path contacted network"))
    dataset = build(source, tmp_path / "original")
    assert build(source, tmp_path / "original").path == dataset.path
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definition(dataset)))
    registry = DatasetRegistry(definition(dataset), cache_root=tmp_path / "cache")
    name = dataset.manifest.dataset_id
    registry.build(name, source_paths=source[1], download_options={"show_progress": False})
    assert registry.inspect(name).status == "missing" and registry.inspect_build(name).status == "available"
    assert dict(registry.open_build(name).counts) == dict(dataset.counts)
    ref, = dataset.manifest.references
    args = ["datasets", "--registry", str(config), "--cache-root", str(tmp_path / "cache"),
            "query", name, "--built", "--reference", ref.key, "--kind", "region",
            "--namespace", "Ingolia2009:feature", "--target", "HRA1", "--unit", "RPKM"]
    run(args + ["--format", "jsonl"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(rows) == 8 and {r["quantity"] for r in rows} == {"rna_abundance", "ribosome_occupancy"}
    assert all(r["target"]["kind"] == "region" and r["origin"] == "derived" for r in rows)
