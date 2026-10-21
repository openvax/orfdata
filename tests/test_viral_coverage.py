"""The published GDC inventory and bundled coding references have distinct scopes."""

import gzip
import hashlib
from importlib.resources import files
from pathlib import Path

import pytest

from orfdata import ViralPanelEntry, viral_panel_coverage, viral_references
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
