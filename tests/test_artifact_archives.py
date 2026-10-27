"""Split downloads retain one independently verified, offline dataset artifact."""

import hashlib
import io
import json
import shutil
import tarfile
from dataclasses import replace
from pathlib import Path

import pytest
from test_measurements import definition, write_dataset

from orfdata import DatasetManifest, DatasetRegistry


def split_archive(source, directory, *, member=None, extra=False):
    directory.mkdir()
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as archive:
        archive.add(source, arcname=member or source.name)
        if extra:
            archive.add(source, arcname="unexpected.sqlite")
    payload = stream.getvalue()
    parts, paths = [], []
    for index, start in enumerate(range(0, len(payload), 256)):
        block = payload[start:start + 256]
        path = directory / f"archive.part{index:04d}"
        path.write_bytes(block)
        paths.append(path)
        parts.append(dict(url="https://example.invalid/" + path.name,
                          sha256=hashlib.sha256(block).hexdigest(), size=len(block)))
    return dict(format="tar.gz", parts=parts, sha256=hashlib.sha256(payload).hexdigest(), size=len(payload)), paths


def archive_release(dataset, tmp_path, **options):
    archive, paths = split_archive(dataset.path, tmp_path / "parts", member="dataset.sqlite", **options)
    definitions = definition(dataset)
    release = definitions[dataset.manifest.dataset_id]["versions"][dataset.manifest.version]
    release["archive"] = archive
    release["assets"]["dataset.sqlite"] = {k: release["assets"]["dataset.sqlite"][k] for k in ("sha256", "size")}
    return definitions, paths


def test_split_and_historical_single_file_releases_reopen_offline(tmp_path, experiment, monkeypatch):
    import requests

    old = write_dataset(tmp_path / "old.sqlite", experiment)
    old_expected = list(old.measurements())
    updated = (replace(experiment[0], version="2"), experiment[1], experiment[2],
               [replace(m, dataset_version="2") for m in experiment[3]])
    new = write_dataset(tmp_path / "new.sqlite", updated)
    expected = list(new.measurements())
    definitions, parts = archive_release(new, tmp_path)
    definitions["fixture"]["versions"].update(definition(old)["fixture"]["versions"])
    cache = tmp_path / "cache"
    registry = DatasetRegistry(definitions, cache_root=cache)
    assert registry.inspect("fixture").status == "missing" and not cache.exists()
    assert registry.download_metadata("fixture")["archive"] == definitions["fixture"]["versions"]["2"]["archive"]
    with pytest.raises(FileNotFoundError):
        registry.open("fixture")
    registry.download("fixture", "1", source_paths={"dataset.sqlite": old.path}, show_progress=False)
    first = registry.download("fixture", source_paths={"dataset.sqlite": parts}, show_progress=False)["dataset.sqlite"]
    second = registry.download("fixture", force=True, source_paths={"dataset.sqlite": parts}, show_progress=False)["dataset.sqlite"]
    assert first != second and first.exists()
    assert registry.inspect("fixture").verified
    assert registry.inspect("fixture", verify_files=False).status == "available"
    assert not registry.inspect("fixture", verify_files=False).verified
    for path in [*parts, old.path, new.path]:
        path.unlink()
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Offline query contacted network"))
    assert list(registry.open("fixture").measurements()) == expected
    assert list(registry.open("fixture", "1").measurements()) == old_expected


@pytest.mark.parametrize("change", ["member_hash", "manifest", "reversed", "missing_part", "extra_member"])
def test_rejected_archives_never_publish_a_dataset(tmp_path, experiment, change):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions, parts = archive_release(dataset, tmp_path, extra=change == "extra_member")
    release = definitions["fixture"]["versions"]["1"]
    if change == "member_hash":
        release["assets"]["dataset.sqlite"]["sha256"] = "0" * 64
    elif change == "manifest":
        release["manifest"]["title"] = "Wrong source scope"
    elif change == "reversed":
        release["archive"]["parts"].reverse()
        parts.reverse()  # Each part still matches; their assembly must fail.
    elif change == "missing_part":
        parts[-1].unlink()
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    with pytest.raises((ValueError, OSError)):
        registry.download("fixture", source_paths={"dataset.sqlite": parts}, show_progress=False)
    assert registry.inspect("fixture").status == "missing"


def test_forced_failed_download_keeps_the_previous_reader(tmp_path, experiment):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions, parts = archive_release(dataset, tmp_path)
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    registry.download("fixture", source_paths={"dataset.sqlite": parts}, show_progress=False)
    old = registry.open("fixture")
    parts[0].write_bytes(b"changed archive part")
    with pytest.raises(ValueError):
        registry.download("fixture", force=True, source_paths={"dataset.sqlite": parts}, show_progress=False)
    assert registry.open("fixture").path == old.path
    assert list(old.measurements()) == list(dataset.measurements())


def test_cache_receipt_cannot_repin_a_modified_payload(tmp_path, experiment):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions, parts = archive_release(dataset, tmp_path)
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    path = registry.download("fixture", source_paths={"dataset.sqlite": parts}, show_progress=False)["dataset.sqlite"]
    data = bytearray(path.read_bytes())
    data[-1] ^= 1
    path.write_bytes(data)
    receipt_path = path.parent / ".datacache-manifest.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["files"]["dataset.sqlite"]["sha256"] = hashlib.sha256(data).hexdigest()
    receipt_path.write_text(json.dumps(receipt))
    assert registry.inspect("fixture").status == "invalid"
    with pytest.raises(ValueError, match="invalid"):
        registry.open("fixture")


@pytest.mark.parametrize("override", [[], "one-path", {"not-dataset.sqlite": []}])
def test_invalid_part_overrides_do_not_create_a_cache(tmp_path, experiment, override):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions, _ = archive_release(dataset, tmp_path)
    cache = tmp_path / "cache"
    registry = DatasetRegistry(definitions, cache_root=cache)
    supplied = override if isinstance(override, dict) else {"dataset.sqlite": override}
    with pytest.raises(ValueError):
        registry.download("fixture", source_paths=supplied)
    assert not cache.exists()


@pytest.mark.parametrize("change", ["format", "unverified_part", "duplicate_url", "total_size", "local_path", "boolean_size"])
def test_unpinned_or_inconsistent_archive_definitions_are_rejected(tmp_path, experiment, change):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions, _ = archive_release(dataset, tmp_path)
    archive = definitions["fixture"]["versions"]["1"]["archive"]
    if change == "format":
        archive["format"] = "unknown"
    elif change == "unverified_part":
        del archive["parts"][0]["sha256"]
    elif change == "duplicate_url":
        archive["parts"][1]["url"] = archive["parts"][0]["url"]
    elif change == "total_size":
        archive["size"] += 1
    elif change == "local_path":
        archive["parts"][0]["path"] = "/unexpected-local-source"
    else:
        archive["parts"][0]["size"] = True
    with pytest.raises(ValueError):
        DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("local_override", [False, True])
def test_archived_source_build_reuses_offline_without_component_cache(tmp_path, monkeypatch, local_override):
    import requests

    root = Path(__file__).parent / "data/pfizer/dog"
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    archived = next(s for s in manifest.sources if s.asset_id == "protein.tsv")
    archive, parts = split_archive(root / archived.asset_id, tmp_path / "parts")
    if not local_override:
        for part, path in zip(archive["parts"], parts):
            part["url"] = path.as_uri()
    manifest = replace(manifest, sources=tuple(replace(s, metadata=dict(archive=archive)) if s == archived else s
                                               for s in manifest.sources))
    definitions = {manifest.dataset_id: dict(default_version=manifest.version, versions={manifest.version: dict(
        manifest=manifest.to_dict(), assets={}, distribution=dict(access_terms="Fixture", redistributable=False))})}
    cache = tmp_path / "cache"
    registry = DatasetRegistry(definitions, cache_root=cache)
    inputs = {s.asset_id: root / s.asset_id for s in manifest.sources if s.asset_id != archived.asset_id}
    if local_override:
        inputs[archived.asset_id] = parts
    built = registry.build(manifest.dataset_id, source_paths=inputs, download_options=dict(show_progress=False))
    expected = list(built.measurements())
    for path in parts:
        path.unlink()
    shutil.rmtree(cache / "orfdata" / "components")
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Cached build contacted network"))
    assert registry.build(manifest.dataset_id).path == built.path
    assert list(registry.open_build(manifest.dataset_id).measurements()) == expected
