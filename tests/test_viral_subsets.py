"""Named viral subsets resolve against the GDC inventory, or fail loudly."""

import pytest

from orfdata import normalize_viral_subsets, viral_panel_coverage, viral_subset_abbreviations
from orfdata.viral_subsets import _check_abbreviations, viral_subsets

PANEL = set(viral_panel_coverage()["abbreviation"])


def test_every_published_subset_resolves_to_real_inventory_rows():
    rows = viral_subsets()
    assert [row["name"] for row in rows] == [
        "cancer-associated", "cancer-causing", "cancer-integrated-driver"]
    for row in rows:
        resolved = viral_subset_abbreviations(row["name"])
        assert set(resolved) <= PANEL
        # The advertised count is what the subset actually selects.
        assert row["references"] == len(resolved) == len(row["abbreviations"])
        assert row["source_urls"] and row["limitations"]
    counts = {row["name"]: row["references"] for row in rows}
    assert counts == {"cancer-associated": 37, "cancer-causing": 20,
                      "cancer-integrated-driver": 15}


def test_the_iarc_tiers_nest():
    causing = set(viral_subset_abbreviations("cancer-causing"))
    associated = set(viral_subset_abbreviations("cancer-associated"))
    assert causing < associated
    assert "HPV16" in causing  # The single most consequential member of the tier.
    assert viral_subset_abbreviations(["cancer-causing", "cancer-integrated-driver"]) == tuple(
        sorted(causing | set(viral_subset_abbreviations("cancer-integrated-driver"))))


def test_a_misspelled_abbreviation_is_a_packaging_error_not_a_silent_filter():
    with pytest.raises(ValueError, match="outside the GDC inventory"):
        _check_abbreviations("test", ["HPV16", "HPV1 6"], PANEL)
    with pytest.raises(ValueError, match="repeats"):
        _check_abbreviations("test", ["HPV16", "HPV16"], PANEL)
    assert _check_abbreviations("test", ["HPV16", "HBV"], PANEL) is None


def test_subset_names_and_aliases_normalize():
    assert normalize_viral_subsets(None) == ()
    assert normalize_viral_subsets("Cancer-Causing") == ("cancer-causing",)
    assert normalize_viral_subsets(["integrated-direct-driver"]) == ("cancer-integrated-driver",)
    assert normalize_viral_subsets(["cancer-causing", "cancer-causing"]) == ("cancer-causing",)
    with pytest.raises(ValueError, match="Unknown viral subset"):
        normalize_viral_subsets(["cancer-causin"])
    with pytest.raises(ValueError, match="nonempty names"):
        normalize_viral_subsets([" "])
