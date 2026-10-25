"""A second collection — another species, another inventory — is data, not a code path."""

import csv
import hashlib
from dataclasses import replace

import pandas as pd
import pytest
from conftest import _genbank, serve_records

from orfdata import (
    download_viral_references,
    viral_panel_coverage,
    viral_reference_summary,
    viral_registry,
)
from orfdata.record_sources import GENBANK_SOURCE, reference_source
from orfdata.reference_ledger import ledger_table, write_ledger
from orfdata.reference_table import PackagedTable, _packaged, write_table
from orfdata.viral_coverage import _coverage
from orfdata.viral_registry import (
    INVENTORY_COLUMNS,
    INVENTORY_NUMERIC,
    INVENTORY_SCHEMA,
    ViralCollection,
)

PANEL = "MTEST.v1"
RECORDS = {
    "MTEST001.1": _genbank("MTEST001.1", "ATG" + "GCT" * 8 + "TAA" + "ACGT" * 6, [
        ("1..30", {"gene": "gagX", "product": "protein gagX", "codon_start": 1,
                   "translation": "MAAAAAAAA"}),
    ], name="MTEST001", description="Mouse test virus 1", taxon=11101),
    "MTEST002.1": _genbank("MTEST002.1", "ATG" + "GCT" * 5 + "TAA" + "ACGT" * 9, [
        ("1..21", {"gene": "envX", "product": "protein envX", "codon_start": 1,
                   "translation": "MAAAAA"}),
    ], name="MTEST002", description="Mouse test virus 2", taxon=11102),
}


def _inventory_row(index, abbreviation, accession):
    return dict(index=index, name=f"Mouse test virus {index}", abbreviation=abbreviation,
                refseq="", genbank=accession, pave_id="", pave_notes="",
                bundled_reference_id="", coding_models="", status="not_bundled")


@pytest.fixture
def mouse_collection(tmp_path, monkeypatch):
    """Register a whole second collection from packaged tables and nothing else."""
    packaged = tmp_path / "packaged"
    packaged.mkdir()

    class _Anchor:
        """Serve this collection's tables from the fixture, the rest from the package."""

        def joinpath(self, resource):
            local = packaged / resource.rsplit("/", 1)[-1]
            if local.exists():
                return local
            from importlib.resources import files as packaged_files

            return packaged_files("orfdata").joinpath(resource)

    monkeypatch.setattr("orfdata.reference_table.files", lambda _package: _Anchor())
    _packaged.cache_clear()

    inventory = write_table(
        pd.DataFrame([_inventory_row(1, "MTV1", "MTEST001.1"),
                      _inventory_row(2, "MTV2", "MTEST002.1")]),
        schema=INVENTORY_SCHEMA, columns=INVENTORY_COLUMNS,
        metadata={"panel": PANEL, "source": "https://example.invalid/mouse-inventory"})
    (packaged / "mtest-inventory.tsv").write_bytes(inventory)
    collection = ViralCollection(
        name="mtest-viruses", title="Mouse test viral screen",
        target_species="mus_musculus", target_taxon_id=10090, panel=PANEL,
        panel_url="https://example.invalid/mouse-inventory",
        panel_sha256="0" * 64, panel_size=1, panel_filename="mouse-inventory.txt",
        contexts=(), primary_source=GENBANK_SOURCE, references=2,
        inventory=PackagedTable(resource="data/viral/mtest-inventory.tsv",
                                sha256=hashlib.sha256(inventory).hexdigest(),
                                schema=INVENTORY_SCHEMA, columns=INVENTORY_COLUMNS,
                                numeric=INVENTORY_NUMERIC))

    serve_records(monkeypatch, RECORDS)
    ledger = write_ledger(
        reference_source(GENBANK_SOURCE).build_ledger(_coverage(collection), timeout=1),
        panel=PANEL, source=GENBANK_SOURCE, source_url="https://example.invalid/efetch",
        checked="2026-01-01T00:00:00Z")
    (packaged / "mtest-genbank.tsv").write_bytes(ledger)
    collection = replace(collection, ledgers={GENBANK_SOURCE: ledger_table(
        "data/viral/mtest-genbank.tsv", hashlib.sha256(ledger).hexdigest())})
    monkeypatch.setitem(viral_registry._COLLECTIONS, collection.name, collection)
    yield collection
    _packaged.cache_clear()


def test_a_second_collection_downloads_without_a_single_new_code_path(mouse_collection, tmp_path):
    inventory = viral_panel_coverage(collection="mtest-viruses")
    assert inventory["abbreviation"].to_list() == ["MTV1", "MTV2"]
    # The panel name resolves to the same collection.
    assert viral_panel_coverage(PANEL).equals(inventory)

    bundle = download_viral_references(species="mouse", collection="mtest-viruses",
                                       cache_root=str(tmp_path / "cache"))
    assert bundle.collection == "mtest-viruses" and bundle.target_species == "mus_musculus"
    rows = {row["abbreviation"]: row for row in
            csv.DictReader(bundle.genome_table.open(), delimiter="\t")}
    assert set(rows) == {"MTV1", "MTV2"}
    # The taxon and the accession version come from this collection's own ledger.
    assert rows["MTV1"]["viral_taxon_id"] == "11101"
    assert rows["MTV2"]["accession_version"] == "MTEST002.1"

    summary = viral_reference_summary(species="mouse", collection="mtest-viruses",
                                      cache_root=str(tmp_path / "cache"))
    assert summary["genome_count"] == 2 and summary["transcript_count"] == 2
    assert summary["target_taxon_id"] == 10090


def test_a_collection_is_refused_for_a_species_it_is_not_for(mouse_collection, tmp_path):
    with pytest.raises((ValueError, SystemExit), match="different target"):
        download_viral_references(species="human", collection="mtest-viruses",
                                  cache_root=str(tmp_path / "cache"))
