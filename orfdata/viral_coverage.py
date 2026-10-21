"""Offline metadata for the published GDC viral reference inventory."""

import csv
import hashlib
import io
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

PANEL = "GRCh38.d1.vd1"
INVENTORY_URL = "https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt"
INVENTORY_SHA256 = "f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56"
_COVERAGE_SHA256 = "0bfa5ffc184eafcdb615b00cbe79acdd6e5e4f22ff2f919a2c31d62eb4cc1492"
_RESOURCE = "data/viral/gdc-grch38-d1-vd1-coverage.tsv"
_HEADER = ("index", "name", "abbreviation", "refseq", "genbank", "pave_id",
           "pave_notes", "bundled_reference_id", "coding_models", "status")


@dataclass(frozen=True)
class ViralPanelEntry:
    """One GDC inventory row and its current ORFdata coding-reference status.

    Inventory accessions are reproduced as published and may omit versions.
    ``not_bundled`` means that ORFdata has not released a qualified coding
    reference for this row; it says nothing about whether source annotations
    exist elsewhere. ``source_cds_bundled`` does not claim an exhaustive ORF set.
    """

    index: int
    name: str
    abbreviation: str
    refseq: str | None
    genbank: str | None
    pave_id: str | None
    pave_notes: str | None
    status: str
    bundled_reference_id: str | None
    coding_models: int | None


@lru_cache(maxsize=1)
def _coverage():
    raw = files("orfdata").joinpath(_RESOURCE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != _COVERAGE_SHA256:
        raise ValueError("Bundled GDC viral coverage checksum mismatch")
    lines = raw.decode("utf-8").splitlines()
    expected_prefix = (
        "# schema: orfdata.viral_panel_coverage.v1",
        f"# panel: {PANEL}",
        f"# source: {INVENTORY_URL}",
        f"# source_sha256: {INVENTORY_SHA256}",
    )
    if tuple(lines[:4]) != expected_prefix:
        raise ValueError("Unsupported GDC viral coverage metadata")
    reader = csv.DictReader(io.StringIO("\n".join(lines[4:])), delimiter="\t")
    if tuple(reader.fieldnames or ()) != _HEADER:
        raise ValueError("Unsupported GDC viral coverage columns")
    entries = tuple(ViralPanelEntry(
        index=int(row["index"]), name=row["name"], abbreviation=row["abbreviation"],
        refseq=row["refseq"] or None, genbank=row["genbank"] or None,
        pave_id=row["pave_id"] or None, pave_notes=row["pave_notes"] or None,
        status=row["status"], bundled_reference_id=row["bundled_reference_id"] or None,
        coding_models=int(row["coding_models"]) if row["coding_models"] else None,
    ) for row in reader)
    if len(entries) != 200 or tuple(entry.index for entry in entries) != tuple(range(1, 201)):
        raise ValueError("Incomplete GDC viral coverage inventory")
    return entries


def viral_panel_coverage(panel=PANEL):
    """Return all GDC viral inventory rows, including references not yet bundled.

    This is an inventory and coverage query, not a selection of protein
    sequences. The released coding references are available separately through
    ``viral_references()``.
    """
    if panel != PANEL:
        raise ValueError(f"Unknown viral panel {panel!r}; available panel: {PANEL}")
    return _coverage()
