"""Source builds and published artifacts remain distinct, explicit offline origins."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry
from orfdata.cli import run
from orfdata.sources.pfizer import import_protein


@pytest.fixture
def experiment(tmp_path):
    from test_measurements import definition

    root = Path(__file__).parent / "data/pfizer/dog"
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    inputs = {}
    for asset in manifest.sources:
        path = tmp_path / asset.asset_id
        path.write_bytes((root / asset.asset_id).read_bytes())
        inputs[asset.asset_id] = path
    artifact = import_protein(inputs, tmp_path / "published.sqlite", manifest)
    definitions = definition(artifact)
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    return registry, manifest, inputs, artifact, definitions


def test_explicit_offline_built_origin_and_force_reader_preservation(experiment, monkeypatch):
    import requests

    registry, manifest, inputs, artifact, _ = experiment
    name = manifest.dataset_id
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Read/build contacted network"))
    assert registry.inspect_build(name).status == "missing"
    with pytest.raises(FileNotFoundError):
        registry.open_build(name)
    built = registry.build(name, source_paths=inputs)
    checked = registry.inspect_build(name)
    assert checked.status == "available" and not checked.verified
    assert all(x["verified"] and x["origin"].startswith("file:") for x in checked.sources.values())
    assert registry.inspect(name).status == "missing"
    with pytest.raises(FileNotFoundError):
        registry.open(name)
    registry.download(name, source_paths={"dataset.sqlite": artifact.path}, show_progress=False)
    published = registry.open(name)
    assert published.path != built.path
    assert registry.open_build(name).path == built.path
    assert [x.to_dict() for x in published.measurements()] == [x.to_dict() for x in built.measurements()]
    replaced = registry.build(name, source_paths=inputs, force=True)
    assert replaced.path != built.path and built.path.exists()
    assert registry.open_build(name).path == replaced.path
    for path in inputs.values():
        path.unlink()  # The build is usable after its temporary overrides disappear.
    assert list(registry.open_build(name).measurements()) == list(built.measurements())
    assert registry.open(name).path == published.path


@pytest.mark.parametrize("change", ("source", "transform", "output", "embedded_manifest"))
def test_corruption_and_identity_changes_are_rejected(experiment, change):
    registry, manifest, inputs, _, _ = experiment
    name = manifest.dataset_id
    built = registry.build(name, source_paths=inputs)
    receipt_path = built.path.parent / ".datacache-manifest.json"
    receipt = json.loads(receipt_path.read_text())
    if change == "source":
        asset = next(iter(receipt["sources"]))
        for spec in (receipt["sources"][asset], receipt["definition"]["sources"][asset]):
            spec["sha256"] = "a" * 64
    elif change == "transform":
        receipt["definition"]["transform"]["options"]["title"] = "Changed transform"
    elif change == "output":
        data = bytearray(built.path.read_bytes())
        data[-1] ^= 1
        built.path.write_bytes(data)
    else:
        # A self-consistent build receipt still cannot change embedded identity.
        other = replace(manifest, title="Wrong embedded manifest")
        built.path.unlink()
        import_protein(inputs, built.path, other)
        receipt["outputs"]["dataset.sqlite"].update(
            sha256=hashlib.sha256(built.path.read_bytes()).hexdigest(), size=built.path.stat().st_size)
    receipt_path.write_text(json.dumps(receipt))
    if change != "embedded_manifest":
        assert registry.inspect_build(name).status == "invalid"
    with pytest.raises(Exception, match="(?i)(invalid|manifest)"):
        registry.open_build(name)


def test_changed_source_manifest_has_separate_build_cache(experiment, tmp_path):
    registry, manifest, inputs, _, definitions = experiment
    registry.build(manifest.dataset_id, source_paths=inputs)
    changed = json.loads(json.dumps(definitions))
    changed[manifest.dataset_id]["versions"][manifest.version]["manifest"]["sources"][0]["sha256"] = "a" * 64
    other = DatasetRegistry(changed, cache_root=tmp_path / "cache")
    assert other.inspect_build(manifest.dataset_id).status == "missing"
    assert registry.inspect_build(manifest.dataset_id, verify_files=False).status == "available"


def test_cli_build_then_explicit_offline_query_and_inspection(experiment, tmp_path, capsys, monkeypatch):
    import requests

    _, manifest, inputs, _, definitions = experiment
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definitions))
    prefix = ["datasets", "--registry", str(config), "--cache-root", str(tmp_path / "cache")]
    name = manifest.dataset_id
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("CLI contacted network"))
    run(prefix + ["inspect", name, "--built"])
    assert json.loads(capsys.readouterr().out)["status"] == "missing"
    overrides = [arg for name, path in inputs.items() for arg in ("--source-path", f"{name}={path}")]
    run(prefix + ["build", name, "--quiet"] + overrides)
    built = json.loads(capsys.readouterr().out)
    run(prefix + ["inspect", name, "--built"])
    inspected = json.loads(capsys.readouterr().out)
    assert inspected["status"] == "available" and inspected["origin"] == "source_build"
    assert not inspected["verified"]
    run(prefix + ["query", name, "--built", "--unit", "iBAQ"])
    records = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert len(records) == built["counts"]["measurement"] // 2
    run(prefix + ["pairs", name, "--built", "--basis", "same_material"])
    assert capsys.readouterr().out == ""
    run(prefix + ["inspect", name])
    assert json.loads(capsys.readouterr().out)["status"] == "missing"
