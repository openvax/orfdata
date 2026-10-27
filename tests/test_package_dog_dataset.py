"""Verified canine artifacts install and replay offline without build inputs."""

import json
from pathlib import Path

import pytest
from test_pfizer_evidence import source as source
from test_pfizer_native import inputs as inputs
from test_verify_dog_dataset import assembled as assembled
from test_verify_dog_dataset import check

from orfdata import DatasetRegistry
from orfdata.search_reference import SearchReference
from orfdata.sources.pfizer_native import file_sha256
from tools import package_dog_dataset, verify_dog_dataset


@pytest.fixture
def accepted(assembled, tmp_path, monkeypatch):
    receipt = check(assembled)
    receipt["tool_sha256"] = file_sha256(verify_dog_dataset.__file__)
    path = tmp_path / "acceptance.json"
    path.write_text(json.dumps(receipt))
    # Explicitly substitute the regression selection, never advertise it as
    # having passed the complete production plan.
    monkeypatch.setattr(package_dog_dataset, "PLAN_PATH", verify_dog_dataset.PLAN_PATH)
    return assembled, path


def pack(accepted, output):
    assembled, receipt = accepted
    return package_dog_dataset.package(assembled[0].path, receipt, output,
                                       url="https://example.invalid/regression.sqlite.gz", data_release="regression-only",
                                       max_bytes=16 * 1024 * 1024)


def test_packaged_artifact_and_public_evidence_replay_without_build_inputs(accepted, tmp_path, monkeypatch):
    (dataset, _, inputs), _ = accepted
    expected = list(dataset.measurements())
    expected_group = dataset.source_evidence().summary_group("A0A5F4CET6")
    archive = tmp_path / "dataset.sqlite.gz"
    result = pack(accepted, archive)
    registry = DatasetRegistry(result["registry"], cache_root=tmp_path / "installed")
    registry.download(dataset.manifest.dataset_id, source_paths={"dataset.sqlite": archive}, show_progress=False)
    # Inputs are private fixture files except the shared FASTA regression slice.
    # Remove only the private files, and make FASTA opening fail if attempted.
    for name in ("evidence.sqlite", "lineage.json", "protein.tsv"):
        Path(inputs[name]).unlink()

    def forbidden(*args, **kwargs):
        raise AssertionError("Ordinary artifact queries must not open source FASTA")

    monkeypatch.setattr(SearchReference, "__init__", forbidden)
    installed = registry.open(dataset.manifest.dataset_id)
    assert file_sha256(installed.path) == result["compression"]["installed_sha256"]
    assert list(installed.measurements()) == expected
    assert installed.source_evidence().summary_group("A0A5F4CET6") == expected_group
    assert installed.source_evidence().search_protein("A0A5F4CET6") is not None
    assert registry.inspect_reference(dataset.manifest.dataset_id).status == "missing"
    assert result["compression"]["download_sha256"] == file_sha256(archive)


def test_stale_acceptance_and_wrong_selection_cannot_be_packaged(accepted, tmp_path, monkeypatch):
    _, receipt = accepted
    original = receipt.read_text()
    for key in ("tool_sha256", "dataset_sha256"):
        altered = json.loads(original)
        altered[key] = "0" * 64
        receipt.write_text(json.dumps(altered))
        with pytest.raises(ValueError, match="matching successful"):
            pack(accepted, tmp_path / "failed.gz")
        assert not (tmp_path / "failed.gz").exists()
    receipt.write_text(original)
    monkeypatch.setattr(package_dog_dataset, "PLAN_PATH", Path(package_dog_dataset.__file__).with_name("dog_result_plan.json"))
    with pytest.raises(ValueError, match="complete pinned"):
        pack(accepted, tmp_path / "partial.gz")


def test_gzip_bounds_reproducibility_and_existing_output_preservation(tmp_path):
    source = tmp_path / "source.bin"
    source.write_bytes(bytes(range(256)) * 100)
    with pytest.raises(ValueError, match="byte limit"):
        package_dog_dataset.compress_verified(source, tmp_path / "failed.gz", max_bytes=100)
    assert not (tmp_path / "failed.gz").exists()
    results = [package_dog_dataset.compress_verified(source, tmp_path / name, max_bytes=1024 * 1024)
               for name in ("one.gz", "two.gz")]
    assert results[0] == results[1]
    with pytest.raises(FileExistsError):
        package_dog_dataset.compress_verified(source, tmp_path / "one.gz", max_bytes=1024 * 1024)
    assert file_sha256(tmp_path / "one.gz") == results[0]["download_sha256"]


@pytest.mark.parametrize("changed", ["dataset", "receipt"])
def test_verified_inputs_cannot_change_during_packaging(accepted, tmp_path, monkeypatch, changed):
    original = package_dog_dataset.compress_verified

    def change_then_compress(source_path, *args, **kwargs):
        path = Path(source_path) if changed == "dataset" else accepted[1]
        with path.open("ab") as handle:
            handle.write(b"\n")
        return original(source_path, *args, **kwargs)

    monkeypatch.setattr(package_dog_dataset, "compress_verified", change_then_compress)
    output = tmp_path / "changed.gz"
    with pytest.raises(ValueError, match="inputs changed"):
        pack(accepted, output)
    assert not output.exists()
