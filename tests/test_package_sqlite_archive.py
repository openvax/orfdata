"""Archive byte preservation, deterministic parts and bounded failure cleanup."""

import os

import pytest
from test_measurements import definition, write_dataset

from orfdata import DatasetRegistry
from tools import package_sqlite_archive as packer


def pack(source, destination, **limits):
    return packer.archive_verified(source, destination, member="dataset.sqlite",
                                   url_prefix="https://example.invalid/release/", part_bytes=512,
                                   max_bytes=limits.get("max_bytes", 256000))


def test_deterministic_parts_install_with_the_existing_managed_archive_reader(tmp_path, experiment):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    first = pack(dataset.path, tmp_path / "one")
    second = pack(dataset.path, tmp_path / "two")
    assert first["archive"] == second["archive"]
    assert len(first["archive"]["parts"]) > 1
    assert all(0 < p["size"] <= 512 for p in first["archive"]["parts"])
    definitions = definition(dataset)
    release = definitions["fixture"]["versions"]["1"]
    release["archive"] = first["archive"]
    release["assets"]["dataset.sqlite"] = dict(size=first["installed_size"], sha256=first["installed_sha256"])
    paths = [tmp_path / "one" / p["url"].rsplit("/", 1)[-1] for p in first["archive"]["parts"]]
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    registry.download("fixture", source_paths={"dataset.sqlite": paths}, show_progress=False)
    assert list(registry.open("fixture").measurements()) == list(dataset.measurements())


def test_output_limit_and_existing_directory_preservation(tmp_path):
    source = tmp_path / "source.bin"
    source.write_bytes(os.urandom(4096))
    with pytest.raises(ValueError, match="byte limit"):
        pack(source, tmp_path / "failed", max_bytes=64)
    assert list((tmp_path / "failed").iterdir()) == []
    existing = tmp_path / "existing"
    existing.mkdir()
    marker = existing / "owned-by-user"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError):
        pack(source, existing)
    assert marker.read_text() == "preserve"


@pytest.mark.parametrize("failure", ["interrupt", "changed_source"])
def test_interrupted_or_changed_source_removes_only_owned_parts(tmp_path, monkeypatch, failure):
    source = tmp_path / "source.bin"
    source.write_bytes(os.urandom(16384))
    if failure == "interrupt":
        original = packer._PartsWriter.write

        def interrupted(self, block):
            if self.count > 1000:
                raise KeyboardInterrupt("interrupted archive")
            return original(self, block)

        monkeypatch.setattr(packer._PartsWriter, "write", interrupted)
        expected = KeyboardInterrupt
    else:
        original = packer._PartsWriter.close

        def changed(self):
            original(self)
            with source.open("ab") as handle:
                handle.write(b"local source edit")

        monkeypatch.setattr(packer._PartsWriter, "close", changed)
        expected = ValueError
    with pytest.raises(expected):
        pack(source, tmp_path / "failed")
    assert source.exists()
    assert list((tmp_path / "failed").iterdir()) == []


@pytest.mark.parametrize("options", [dict(part_bytes=2**31), dict(part_bytes=True), dict(max_bytes=0),
                                     dict(member="../dataset.sqlite"), dict(url_prefix="https://example.invalid/no-slash"),
                                     dict(expected_sha256="0" * 64)])
def test_invalid_archive_configuration_writes_nothing(tmp_path, options):
    source = tmp_path / "source.bin"
    source.write_bytes(b"small source")
    arguments = dict(member="dataset.sqlite", url_prefix="https://example.invalid/", part_bytes=1024, max_bytes=1024)
    arguments.update(options)
    with pytest.raises(ValueError):
        packer.archive_verified(source, tmp_path / "output", **arguments)
    assert not (tmp_path / "output").exists()
