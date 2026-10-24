"""Offline metadata for a collection's published viral inventory.

Every inventory query names a collection, by its own name or by the panel it
follows, and reads that collection's packaged table. The GDC human screen is
the default because it is the one collection published today.
"""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .viral_registry import (
    GDC_VIRUSES,
    ViralCollection,
    viral_collection,
    viral_collections,
)

PANEL = GDC_VIRUSES.panel
INVENTORY_URL = GDC_VIRUSES.panel_url
INVENTORY_SHA256 = GDC_VIRUSES.panel_sha256
INVENTORY_SIZE = GDC_VIRUSES.panel_size
_CONTEXTS = {
    "cancer-transcriptome-screening": {
        "title": "GDC human cancer transcriptome viral screen",
        "description": "Membership in the GDC GRCh38.d1.vd1 viral screening inventory for human analysis",
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
    """One sourced, nonexclusive target-species viral-reference selection context."""

    name: str
    title: str
    description: str
    source: str
    references: int
    aliases: tuple[str, ...] = ()


@dataclass(frozen=True)
class ViralPanelEntry:
    """One human-targeted GDC inventory row and its ORFdata coding-reference status.

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


def _collection_for(panel=None, collection=None):
    """The collection a panel name or a collection name refers to."""
    if collection is not None:
        return collection if isinstance(collection, ViralCollection) else viral_collection(
            collection)
    if panel is None:
        return GDC_VIRUSES
    for candidate in viral_collections():
        if panel in (candidate.panel, candidate.name):
            return candidate
    raise ValueError(f"Unknown viral panel {panel!r}; available: "
                     f"{', '.join(sorted(item.panel for item in viral_collections()))}")


@lru_cache(maxsize=4)
def _coverage(collection=GDC_VIRUSES):
    """A collection's inventory rows, from its packaged, checksummed table."""
    inventory = collection.inventory
    if inventory.metadata().get("panel") != collection.panel:
        raise ValueError(f"{inventory.name} is not the {collection.panel} inventory")
    # Membership in the inventory is what puts a row in the collection's contexts.
    entries = tuple(ViralPanelEntry(
        index=int(row["index"]), name=row["name"], abbreviation=row["abbreviation"],
        refseq=row["refseq"] or None, genbank=row["genbank"] or None,
        pave_id=row["pave_id"] or None, pave_notes=row["pave_notes"] or None,
        status=row["status"], bundled_reference_id=row["bundled_reference_id"] or None,
        coding_models=int(row["coding_models"]) if row["coding_models"] else None,
        contexts=collection.contexts,
    ) for row in inventory.read())
    if (len(entries) != collection.references
            or tuple(entry.index for entry in entries) != tuple(
                range(1, collection.references + 1))):
        raise ValueError(f"Incomplete {collection.panel} inventory")
    return entries


def viral_contexts(panel=None, *, collection=None):
    """List selectable contexts without downloading data.

    Context membership is many-to-many and describes a sourced selection use;
    it does not change a reference's accession/version identity or assert
    disease causality.
    """
    entries = _coverage(_collection_for(panel, collection))
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


def viral_panel_coverage(panel=None, *, collection=None, available=None, virus=None,
                         search=None, contexts=None, exclude_contexts=None,
                         context_match="any"):
    """Return every inventory row of one collection, named by panel or by collection.

    Includes references whose source CDS records are not yet bundled in ORFdata.
    This is an inventory and coverage query, not a selection of protein
    sequences. The released coding references are available separately through
    ``viral_references()``.
    """
    chosen = _collection_for(panel, collection)
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
    return tuple(entry for entry in _coverage(chosen)
                 if (available is None or entry.available == available)
                 and (virus is None or entry.abbreviation.casefold() == virus.casefold())
                 and (query is None or any(query in (value or "").casefold() for value in (
                     entry.name, entry.abbreviation, entry.refseq, entry.genbank, entry.pave_id)))
                 and (not selected_contexts or (
                     any(context in entry.contexts for context in selected_contexts)
                     if context_match == "any"
                     else all(context in entry.contexts for context in selected_contexts)))
                 and not any(context in entry.contexts for context in excluded_contexts))


def download_viral_inventory(*, collection=None, cache_root=None, force=False, timeout=60):
    """Fetch a collection's own published inventory metadata through DataCache.

    This downloads the inventory text the collection follows, not viral genomes
    or proteins. The published bytes are checked by SHA-256 and size on both
    download and cache reuse. ``viral_panel_coverage()`` works offline without
    this download.
    """
    from .datasets import _datacache, get_data_dir

    chosen = _collection_for(collection=collection)
    cache = get_data_dir(cache_root) / "viral" / "inventories" / chosen.panel
    return Path(_datacache().fetch_file(
        chosen.panel_url, filename=chosen.panel_filename, cache_root=cache,
        expected_sha256=chosen.panel_sha256, expected_size=chosen.panel_size,
        max_bytes=chosen.panel_size, raw=True, record_provenance=True,
        force=force, timeout=timeout))
