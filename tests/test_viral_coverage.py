"""The published GDC inventory and bundled coding references have distinct scopes."""

import gzip
import hashlib
from importlib.resources import files
from pathlib import Path

import pytest

from orfdata import (
    ViralPanelEntry,
    download_viral_inventory,
    viral_panel_coverage,
    viral_references,
)
from tools.build_gdc_viral_coverage import SOURCE_SHA256, build

FIXTURE = Path(__file__).parent / "data/viral/GRCh83.d1.vd1_virus_decoy.txt.gz"


def test_coverage_rebuilds_from_the_pinned_gdc_inventory():
    raw = gzip.decompress(FIXTURE.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA256
    packaged = files("orfdata").joinpath(
        "data/viral/gdc-grch38-d1-vd1-coverage.tsv").read_bytes()
    assert build(raw) == packaged


def test_coverage_lists_every_gdc_row_and_preserves_missingness():
    entries = viral_panel_coverage()
    assert len(entries) == 200
    assert all(isinstance(entry, ViralPanelEntry) for entry in entries)
    assert len({entry.abbreviation for entry in entries}) == 200
    assert sum(entry.abbreviation.startswith("HPV") for entry in entries) == 189
    assert {entry.status for entry in entries} == {"source_cds_bundled", "not_bundled"}
    assert sum(entry.status == "source_cds_bundled" for entry in entries) == 1
    assert all(entry.bundled_reference_id is None and entry.coding_models is None
               for entry in entries if entry.status == "not_bundled")
    hpv16, = (entry for entry in entries if entry.abbreviation == "HPV16")
    assert hpv16.refseq == "NC_001526"  # The inventory omits the version.
    assert hpv16.pave_id == "HPV16REF"
    assert hpv16.bundled_reference_id == "GRCh38.d1.vd1:HPV16"
    assert hpv16.coding_models == 9


def test_coverage_agrees_with_released_gdc_coding_references():
    bundled = {ref.id: ref.coding_models for ref in viral_references()
               if ref.id.startswith("GRCh38.d1.vd1:")}
    covered = {entry.bundled_reference_id: entry.coding_models
               for entry in viral_panel_coverage() if entry.bundled_reference_id}
    assert covered == bundled
    with pytest.raises(ValueError, match="Unknown viral panel"):
        viral_panel_coverage("GRCh38.d1.v1")


def test_coverage_filters_by_availability_and_source_name():
    hpv16, = viral_panel_coverage(available=True)
    assert hpv16.available and hpv16.abbreviation == "HPV16"
    assert len(viral_panel_coverage(available=False)) == 199
    cmv, = viral_panel_coverage(search="nc_006273.2")
    assert cmv.abbreviation == "CMV ( HHV-5)" and not cmv.available
    hpv18, = viral_panel_coverage(available=False, virus="hpv18")
    assert hpv18.name == "Human papillomavirus 18"
    for bad in ("", "   ", True):
        with pytest.raises(ValueError, match="search"):
            viral_panel_coverage(search=bad)
    with pytest.raises(ValueError, match="available"):
        viral_panel_coverage(available="yes")
    with pytest.raises(ValueError, match="virus"):
        viral_panel_coverage(virus="")


def test_inventory_download_uses_datacache_with_the_published_pin(monkeypatch, tmp_path):
    import orfdata.datasets as datasets

    called = {}

    class FakeDataCache:
        def fetch_file(self, url, **kwargs):
            called.update(url=url, **kwargs)
            return tmp_path / "inventory.txt"

    monkeypatch.setattr(datasets, "_datacache", lambda: FakeDataCache())
    monkeypatch.setattr(datasets, "get_data_dir", lambda root: tmp_path / "orfdata")
    assert download_viral_inventory(cache_root=tmp_path, force=True, timeout=12) == (
        tmp_path / "inventory.txt")
    assert called["url"] == "https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt"
    assert called["expected_sha256"] == SOURCE_SHA256
    assert called["expected_size"] == called["max_bytes"] == len(gzip.decompress(FIXTURE.read_bytes()))
    assert called["raw"] and called["record_provenance"] and called["force"]
    assert called["timeout"] == 12
