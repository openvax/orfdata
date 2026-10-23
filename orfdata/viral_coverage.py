"""Offline metadata for the published GDC viral reference inventory."""

import csv
import hashlib
import io
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

PANEL = "GRCh38.d1.vd1"
INVENTORY_URL = "https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt"
INVENTORY_SHA256 = "f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56"
INVENTORY_SIZE = 13398
_COVERAGE_SHA256 = "0bfa5ffc184eafcdb615b00cbe79acdd6e5e4f22ff2f919a2c31d62eb4cc1492"
_RESOURCE = "data/viral/gdc-grch38-d1-vd1-coverage.tsv"
_HEADER = ("index", "name", "abbreviation", "refseq", "genbank", "pave_id",
           "pave_notes", "bundled_reference_id", "coding_models", "status")
_CONTEXTS = {
    "cancer-transcriptome-screening": {
        "title": "GDC cancer transcriptome screening references",
        "description": "Membership in the GDC GRCh38.d1.vd1 viral screening inventory",
        "source": INVENTORY_URL,
        "aliases": ("cancer", "cancer-screening", "gdc-cancer"),
    },
}
_CONTEXT_ALIASES = {
    alias: name
    for name, metadata in _CONTEXTS.items()
    for alias in metadata.get("aliases", ())
}
@dataclass(frozen=True)
class ViralContext:
    """One sourced, nonexclusive viral-reference selection context."""

    name: str
    title: str
    description: str
    source: str
    references: int
    aliases: tuple[str, ...] = ()


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
    contexts: tuple[str, ...] = ()

    @property
    def available(self):
        """Whether ORFdata bundles source CDS models for this panel row."""
        return self.bundled_reference_id is not None


def _entry_contexts(abbreviation):
    return ("cancer-transcriptome-screening",)


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
        contexts=_entry_contexts(row["abbreviation"]),
    ) for row in reader)
    if len(entries) != 200 or tuple(entry.index for entry in entries) != tuple(range(1, 201)):
        raise ValueError("Incomplete GDC viral coverage inventory")
    return entries


def viral_contexts():
    """List selectable contexts without downloading data.

    Context membership is many-to-many and describes a sourced selection use;
    it does not change a reference's accession/version identity or assert
    disease causality.
    """
    entries = _coverage()
    return tuple(ViralContext(
        name=name,
        title=metadata["title"],
        description=metadata["description"],
        source=metadata["source"],
        references=sum(name in entry.contexts for entry in entries),
        aliases=metadata.get("aliases", ()),
    ) for name, metadata in sorted(_CONTEXTS.items()))


def normalize_viral_contexts(contexts, *, allow_unknown=False):
    if contexts is None:
        return ()
    contexts = (contexts,) if isinstance(contexts, str) else tuple(contexts)
    if not contexts or any(not isinstance(value, str) or not value.strip() for value in contexts):
        raise ValueError("contexts must contain one or more nonempty names")
    normalized = tuple(dict.fromkeys(
        _CONTEXT_ALIASES.get(value.strip().lower(), value.strip().lower())
        for value in contexts))
    unknown = sorted(set(normalized) - set(_CONTEXTS))
    if unknown and not allow_unknown:
        raise ValueError(
            f"Unknown viral context(s) {', '.join(unknown)}; available: "
            f"{', '.join(sorted(_CONTEXTS))}"
        )
    return normalized


def viral_panel_coverage(panel=PANEL, *, available=None, virus=None, search=None,
                         contexts=None, exclude_contexts=None, context_match="any"):
    """Return all GDC viral inventory rows, including references not yet bundled.

    This is an inventory and coverage query, not a selection of protein
    sequences. The released coding references are available separately through
    ``viral_references()``.
    """
    if panel != PANEL:
        raise ValueError(f"Unknown viral panel {panel!r}; available panel: {PANEL}")
    if available is not None and type(available) is not bool:
        raise ValueError("available must be True, False, or None")
    if virus is not None and (not isinstance(virus, str) or not virus.strip()):
        raise ValueError("virus must be a nonempty GDC inventory abbreviation")
    if search is not None and (not isinstance(search, str) or not search.strip()):
        raise ValueError("search must be a nonempty string")
    selected_contexts = normalize_viral_contexts(contexts)
    excluded_contexts = normalize_viral_contexts(exclude_contexts)
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    query = search.casefold() if search is not None else None
    return tuple(entry for entry in _coverage()
                 if (available is None or entry.available == available)
                 and (virus is None or entry.abbreviation.casefold() == virus.casefold())
                 and (query is None or any(query in (value or "").casefold() for value in (
                     entry.name, entry.abbreviation, entry.refseq, entry.genbank, entry.pave_id)))
                 and (not selected_contexts or (
                     any(context in entry.contexts for context in selected_contexts)
                     if context_match == "any"
                     else all(context in entry.contexts for context in selected_contexts)))
                 and not any(context in entry.contexts for context in excluded_contexts))


def download_viral_inventory(*, cache_root=None, force=False, timeout=60):
    """Fetch the pinned GDC inventory metadata through DataCache.

    This downloads GDC's inventory text, not viral genomes or proteins. The
    published bytes are checked by SHA-256 and size on both download and cache
    reuse. ``viral_panel_coverage()`` works offline without this download.
    """
    from .datasets import _datacache, get_data_dir

    cache = get_data_dir(cache_root) / "viral" / "inventories" / PANEL
    return Path(_datacache().fetch_file(
        INVENTORY_URL, filename="GRCh83.d1.vd1_virus_decoy.txt", cache_root=cache,
        expected_sha256=INVENTORY_SHA256, expected_size=INVENTORY_SIZE,
        max_bytes=INVENTORY_SIZE, raw=True, record_provenance=True,
        force=force, timeout=timeout))
