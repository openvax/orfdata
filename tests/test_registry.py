"""Shared registry base and build-receipt helpers (#185)."""

from pathlib import Path
from types import SimpleNamespace

import datacache
import pytest

from orfdata import CatalogRegistry, DatasetRegistry
from orfdata.cache import content_digest, get_data_dir
from orfdata.registry import (
    ReleaseRegistry,
    build_published_registry,
    build_store_path,
    inspect_build_receipt,
    resolve_species_taxon,
)


def test_catalog_and_dataset_registries_share_release_registry_base():
    assert issubclass(CatalogRegistry, ReleaseRegistry)
    assert issubclass(DatasetRegistry, ReleaseRegistry)


def test_resolve_species_taxon_handles_aliases_and_conflicts():
    assert resolve_species_taxon() is None
    assert resolve_species_taxon(taxon_id=9606) == 9606
    assert resolve_species_taxon(species="human") == 9606
    assert resolve_species_taxon(species="human", taxon_id=9606) == 9606
    with pytest.raises(ValueError, match="Unknown species alias"):
        resolve_species_taxon(species="unicorn")
    with pytest.raises(ValueError, match="Species and taxon_id disagree"):
        resolve_species_taxon(species="human", taxon_id=10090)


def test_build_published_registry_falls_back_default_version_and_filters(tmp_path):
    definitions = {
        "example": {
            "default_version": "2",
            "versions": {
                "1": {
                    "assets": {
                        "item.bin": {
                            "url": "https://example.invalid/v1.bin",
                            "sha256": "a" * 64,
                            "size": 1,
                        }
                    }
                },
                "2": {"assets": {}},
                "3": {
                    "archive": {"parts": []},
                    "assets": {
                        "item.bin": {
                            "url": "https://example.invalid/v3.bin",
                            "sha256": "b" * 64,
                            "size": 1,
                        }
                    },
                },
            },
        }
    }
    registry = build_published_registry(
        definitions,
        "items",
        cache_root=tmp_path,
        include_release=lambda release: "archive" not in release,
    )
    inspected = registry.inspect("example", "1")
    expected_path = (
        get_data_dir(tmp_path)
        / "items"
        / "example"
        / "1"
        / content_digest(definitions["example"]["versions"]["1"])
    )
    assert inspected.status == "missing"
    assert Path(inspected.path) == expected_path


def test_inspect_build_receipt_validates_inventory_integrity_and_transform(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("hello")
    source_bytes = source.read_bytes()
    import hashlib

    sha = hashlib.sha256(source_bytes).hexdigest()
    size = len(source_bytes)

    class StubManifest:
        schema = "test.manifest.v1"
        version = "1"

        def __init__(self, *, sha256=sha, size_bytes=size, tag="ok"):
            self.sources = (SimpleNamespace(asset_id="src", sha256=sha256, size=size_bytes),)
            self._tag = tag

        def to_dict(self):
            return {"id": "stub", "version": self.version, "tag": self._tag}

    manifest = StubManifest()
    store = build_store_path("test-imports", "stub", manifest, cache_root=tmp_path)
    assert inspect_build_receipt(store, manifest, "out.bin").status == "missing"

    datacache.materialize(
        store,
        {"src": {"path": str(source), "sha256": sha, "size": size}},
        transform=dict(version=manifest.schema, options=manifest.to_dict()),
        outputs={"out.bin": {}},
        builder=lambda inputs, outputs: Path(outputs["out.bin"]).write_bytes(b"done"),
    )
    valid = inspect_build_receipt(store, manifest, "out.bin")
    assert valid.status == "available" and valid.verified is False

    wrong_output = inspect_build_receipt(
        store, manifest, "other.bin", inventory_error="custom inventory error"
    )
    assert wrong_output.status == "invalid"
    assert "custom inventory error" in str(wrong_output.error)

    wrong_pin = inspect_build_receipt(
        store, StubManifest(sha256="0" * 64), "out.bin", integrity_error="custom integrity error"
    )
    assert wrong_pin.status == "invalid"
    assert "custom integrity error" in str(wrong_pin.error)

    wrong_transform = inspect_build_receipt(
        store, StubManifest(tag="changed"), "out.bin", transform_error="custom transform error"
    )
    assert wrong_transform.status == "invalid"
    assert "custom transform error" in str(wrong_transform.error)
