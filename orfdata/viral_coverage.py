"""Offline metadata for a collection's published viral inventory.

Every inventory query names a collection, by its own name or by the panel it
follows, and reads that collection's packaged table. The GDC human screen is
the default because it is the one collection published today.

An inventory reads as a data frame with a default integer index, one row per
reference the panel screens for, whose columns are :data:`COVERAGE_COLUMNS`:

========================  ======  ======================================================
column                    dtype   meaning
========================  ======  ======================================================
``index``                 Int64   the row's 1-based position in the published inventory
``name``                  str     the virus name as the panel publishes it
``abbreviation``          str     the panel's short label, unique within a collection
``refseq``                str     the RefSeq accession the panel lists, or ``""``
``genbank``               str     the GenBank accession the panel lists, or ``""``
``pave_id``               str     the PaVE reference the panel lists, or ``""``
``pave_notes``            str     the panel's own note about that PaVE reference
``bundled_reference_id``  str     the ORFdata coding reference for this row, or ``""``
``coding_models``         Int64   how many source CDS models that reference carries
``status``                str     ``source_cds_bundled`` or ``not_bundled``
``available``             bool    whether ``bundled_reference_id`` is set
``contexts``              object  tuple of context names this row belongs to
========================  ======  ======================================================

:func:`viral_contexts` returns a second frame, columns :data:`CONTEXT_COLUMNS`:
``name``, ``title``, ``description`` and ``source`` as str, ``references`` as
the number of inventory rows in that context, and ``aliases`` as a tuple of
the shorter names it also answers to.

Filter either with the arguments here or with pandas; ``frame.attrs`` carries
the packaged table's own header, including which panel it follows.

Inventory accessions are reproduced as published and may omit versions. A
``status`` of ``not_bundled`` means that ORFdata has not released a qualified
coding reference for that row; it says nothing about whether source
annotations exist elsewhere. ``source_cds_bundled`` does not claim an
exhaustive ORF set.
"""

from functools import lru_cache
from pathlib import Path

import pandas as pd

from .viral_registry import (
    GDC_VIRUSES,
    INVENTORY_COLUMNS,
    viral_collection,
)

PANEL = GDC_VIRUSES.panel
INVENTORY_URL = GDC_VIRUSES.panel_url
INVENTORY_SHA256 = GDC_VIRUSES.panel_sha256
INVENTORY_SIZE = GDC_VIRUSES.panel_size
# What a read inventory carries: its packaged columns, then the two this
# module derives. The module docstring says what each one means.
COVERAGE_COLUMNS = INVENTORY_COLUMNS + ("available", "contexts")
CONTEXT_COLUMNS = ("name", "title", "description", "source", "references", "aliases")
# The columns a search looks through, being the ways a row can be named.
SEARCHED = ("name", "abbreviation", "refseq", "genbank", "pave_id")
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


def _collection_for(panel=None, collection=None):
    """The collection a panel name or a collection name refers to."""
    if panel is not None and collection is not None:
        by_panel = viral_collection(panel)
        by_name = viral_collection(collection)
        if by_panel is not by_name:
            raise ValueError(
                f"Viral panel {panel!r} and collection {by_name.name!r} disagree"
            )
        return by_name
    try:
        return viral_collection(collection if collection is not None else panel)
    except ValueError as error:
        if panel is not None and collection is None:
            raise ValueError(str(error).replace("viral collection", "viral panel")) from None
        raise


@lru_cache(maxsize=4)
def _coverage(collection=GDC_VIRUSES):
    """A collection's inventory rows, from its packaged, checksummed table."""
    frame = collection.inventory()
    if frame.attrs.get("panel") != collection.panel:
        raise ValueError(f"{collection.inventory_file.name} is not the {collection.panel} inventory")
    if (len(frame) != collection.references
            or frame["index"].to_list() != list(range(1, collection.references + 1))):
        raise ValueError(f"Incomplete {collection.panel} inventory")
    # Whether ORFdata bundles source CDS models for a row is what the inventory
    # is usually asked, so it is a column rather than a rule each caller writes.
    frame["available"] = frame["bundled_reference_id"] != ""
    # Membership in the inventory is what puts a row in the collection's contexts.
    frame["contexts"] = [collection.contexts] * len(frame)
    return frame


def viral_contexts(panel=None, *, collection=None):
    """List selectable contexts without downloading data, as :data:`CONTEXT_COLUMNS`.

    Context membership is many-to-many and describes a sourced selection use;
    it does not change a reference's accession/version identity or assert
    disease causality.
    """
    held = _coverage(_collection_for(panel, collection))["contexts"]
    return pd.DataFrame([{
        "name": name,
        "title": metadata["title"],
        "description": metadata["description"],
        "source": metadata["source"],
        "references": int(held.map(lambda row, name=name: name in row).sum()),
        "aliases": metadata.get("aliases", ()),
    } for name, metadata in sorted(_CONTEXTS.items())])


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
    """Every inventory row of one collection, named by panel or by collection.

    Returns a frame whose columns are :data:`COVERAGE_COLUMNS`, described in
    this module's docstring. Includes references whose source CDS records are not yet bundled in
    ORFdata. This is an inventory and coverage query, not a selection of
    protein sequences. The released coding references are available separately
    through ``viral_references()``.
    """
    chosen = _collection_for(panel, collection)
    if available is not None and type(available) is not bool:
        raise ValueError("available must be True, False, or None")
    if virus is not None and (not isinstance(virus, str) or not virus.strip()):
        raise ValueError("virus must be a nonempty GDC inventory abbreviation")
    if search is not None and (not isinstance(search, str) or not search.strip()):
        raise ValueError("search must be a nonempty string")
    selected = normalize_viral_contexts(contexts)
    excluded = normalize_viral_contexts(exclude_contexts)
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")

    frame = _coverage(chosen)
    keep = pd.Series(True, index=frame.index)
    if selected:
        match = any if context_match == "any" else all
        keep &= frame["contexts"].map(
            lambda held: match(context in held for context in selected))
    if excluded:
        keep &= ~frame["contexts"].map(
            lambda held: any(context in held for context in excluded))
    if available is not None:
        keep &= frame["available"] == available
    if virus is not None:
        keep &= frame["abbreviation"].str.casefold() == virus.casefold()
    if search is not None:
        query = search.casefold()
        found = pd.Series(False, index=frame.index)
        for column in SEARCHED:
            found |= frame[column].str.casefold().str.contains(query, regex=False)
        keep &= found
    return frame[keep]


def download_viral_inventory(*, collection=None, cache_root=None, force=False, timeout=60):
    """Fetch a collection's own published inventory metadata through DataCache.

    This downloads the inventory text the collection follows, not viral genomes
    or proteins. The published bytes are checked by SHA-256 and size on both
    download and cache reuse. ``viral_panel_coverage()`` works offline without
    this download.
    """
    from .datasets import _datacache, _resumable, get_data_dir

    chosen = _collection_for(collection=collection)
    cache = get_data_dir(cache_root) / "viral" / "inventories" / chosen.panel
    return Path(_datacache().fetch_file(
        chosen.panel_url, filename=chosen.panel_filename, cache_root=cache,
        expected_sha256=chosen.panel_sha256, expected_size=chosen.panel_size,
        max_bytes=chosen.panel_size, raw=True, record_provenance=True,
        force=force, **_resumable(dict(timeout=timeout))))
