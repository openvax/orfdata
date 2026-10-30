"""The published GDC inventory and bundled coding references have distinct scopes."""

import gzip
import hashlib
from importlib.resources import files
from pathlib import Path

import pandas as pd
import pytest

from orfdata import (
    download_viral_inventory,
    viral_contexts,
    viral_panel_coverage,
    viral_references,
)
from orfdata.viral_coverage import CONTEXT_COLUMNS, COVERAGE_COLUMNS
from tools.build_gdc_viral_coverage import SOURCE_SHA256, build

FIXTURE = Path(__file__).parent / "data/viral/GRCh83.d1.vd1_virus_decoy.txt.gz"


def test_coverage_rebuilds_from_the_pinned_gdc_inventory():
    raw = gzip.decompress(FIXTURE.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA256
    packaged = files("orfdata").joinpath(
        "data/viral/gdc-grch38-d1-vd1-coverage.tsv").read_bytes()
    assert build(raw) == packaged


def test_coverage_lists_every_gdc_row_and_preserves_missingness():
    inventory = viral_panel_coverage()
    assert len(inventory) == 200
    assert inventory["abbreviation"].nunique() == 200
    assert int(inventory["abbreviation"].str.startswith("HPV").sum()) == 189
    assert set(inventory["status"]) == {"source_cds_bundled", "not_bundled"}
    assert int((inventory["status"] == "source_cds_bundled").sum()) == 1
    # A row with no bundled reference says so by leaving both fields empty.
    unbundled = inventory[inventory["status"] == "not_bundled"]
    assert (unbundled["bundled_reference_id"] == "").all()
    assert unbundled["coding_models"].isna().all()
    hpv16 = inventory.set_index("abbreviation").loc["HPV16"]
    assert hpv16["refseq"] == "NC_001526"  # The inventory omits the version.
    assert hpv16["pave_id"] == "HPV16REF"
    assert hpv16["bundled_reference_id"] == "GRCh38.d1.vd1:HPV16"
    assert hpv16["coding_models"] == 9


def test_coverage_agrees_with_released_gdc_coding_references():
    bundled = {ref.id: ref.coding_models for ref in viral_references()
               if ref.id.startswith("GRCh38.d1.vd1:")}
    inventory = viral_panel_coverage()
    released = inventory[inventory["bundled_reference_id"] != ""]
    assert dict(zip(released["bundled_reference_id"], released["coding_models"])) == bundled
    with pytest.raises(ValueError, match="Unknown viral panel"):
        viral_panel_coverage("GRCh38.d1.v1")


def test_coverage_filters_by_availability_and_source_name():
    hpv16 = viral_panel_coverage(available=True).squeeze()
    assert hpv16["available"] and hpv16["abbreviation"] == "HPV16"
    assert len(viral_panel_coverage(available=False)) == 199
    cmv = viral_panel_coverage(search="nc_006273.2").squeeze()
    assert cmv["abbreviation"] == "CMV ( HHV-5)" and not cmv["available"]
    hpv18 = viral_panel_coverage(available=False, virus="hpv18").squeeze()
    assert hpv18["name"] == "Human papillomavirus 18"
    for bad in ("", "   ", True):
        with pytest.raises(ValueError, match="search"):
            viral_panel_coverage(search=bad)
    with pytest.raises(ValueError, match="available"):
        viral_panel_coverage(available="yes")
    with pytest.raises(ValueError, match="virus"):
        viral_panel_coverage(virus="")


def test_a_coverage_query_hands_back_a_frame_callers_can_work_on():
    inventory = viral_panel_coverage()
    assert isinstance(inventory, pd.DataFrame)
    assert inventory.attrs["panel"] == "GRCh38.d1.vd1"
    # Editing a result never reaches the packaged table behind it.
    inventory.loc[inventory.index[0], "name"] = "scribbled on"
    assert viral_panel_coverage().at[0, "name"] == "Human Cytomegalovirus"


def test_every_frame_carries_exactly_the_columns_it_documents():
    """The documented column set is the API, so a filtered frame keeps it too."""
    inventory = viral_panel_coverage()
    assert tuple(inventory.columns) == COVERAGE_COLUMNS
    assert tuple(viral_panel_coverage(available=True).columns) == COVERAGE_COLUMNS
    assert tuple(viral_panel_coverage(virus="nothing-matches-this").columns) == COVERAGE_COLUMNS
    assert tuple(viral_contexts().columns) == CONTEXT_COLUMNS
    # The two derived columns say what they promise about the packaged ones.
    assert (inventory["available"] == (inventory["bundled_reference_id"] != "")).all()
    assert inventory["contexts"].map(type).eq(tuple).all()
    assert viral_contexts()["aliases"].map(type).eq(tuple).all()


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
    assert called["resume"] == "auto"
