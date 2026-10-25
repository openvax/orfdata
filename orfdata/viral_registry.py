"""The viral reference collections ORFdata publishes, as data rather than code paths.

A collection says what it is for and where it comes from: the species whose
analysis it supports, the published inventory it follows, the sources its
molecules can be fetched from, and which source to prefer where the inventory
names two references that are not the same molecule.

Everything a collection needs offline is packaged and checksummed: the
inventory itself, and one source summary recording what each source
publishes for every entry in it. They all have the same shape, so a second
collection — another species, another inventory, another pair of sources — is
a new entry here, not a new code path.

One collection is published today, ``gdc-viruses``, for human analysis.
"""

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping, Optional

import pandas as pd

from .record_sources import GENBANK_SOURCE, PAVE_SOURCE, reference_source
from .reference_ledger import ExpectedSequence, ledger_table, read_ledger
from .reference_table import CheckedTSV

INVENTORY_SCHEMA = "orfdata.viral_panel_coverage.v1"
INVENTORY_COLUMNS = ("index", "name", "abbreviation", "refseq", "genbank", "pave_id",
                     "pave_notes", "bundled_reference_id", "coding_models", "status")
INVENTORY_NUMERIC = ("index", "coding_models")


@dataclass(frozen=True)
class ViralCollection:
    """One published collection of viral references for a target species."""

    name: str
    title: str
    target_species: str
    target_taxon_id: int
    panel: str
    panel_url: str
    panel_sha256: str
    panel_size: int
    panel_filename: str
    contexts: tuple[str, ...]
    primary_source: str
    inventory_file: CheckedTSV = field(repr=False, hash=False)
    references: int
    source_files: Mapping[str, CheckedTSV] = field(
        default_factory=dict, repr=False, hash=False
    )
    alternate_source: Optional[str] = None
    alternate_note: str = ""
    expected_sequences: tuple[ExpectedSequence, ...] = ()
    aliases: tuple[str, ...] = field(default=(), repr=False)

    def __post_init__(self):
        """Reject incomplete declarations before any download can use them."""
        if type(self.target_taxon_id) is not int or self.target_taxon_id < 1:
            raise ValueError("target_taxon_id must be a positive NCBI taxonomy ID")
        if type(self.references) is not int or self.references < 1:
            raise ValueError("references must be a positive row count")
        if type(self.panel_size) is not int or self.panel_size < 1:
            raise ValueError("panel_size must be a positive byte count")
        try:
            digest_is_hex = len(self.panel_sha256) == 64 and int(self.panel_sha256, 16) >= 0
        except (TypeError, ValueError):
            digest_is_hex = False
        if not digest_is_hex:
            raise ValueError("panel_sha256 must be a SHA-256 hex digest")
        sources = dict(self.source_files)
        for source in {self.primary_source, *sources}:
            reference_source(source)
        if self.primary_source not in sources:
            raise ValueError("A collection must include a summary for its primary source")
        if self.alternate_source is not None and self.alternate_source not in sources:
            raise ValueError("alternate_source must have a summary in this collection")
        unknown_checks = sorted({expected.source for expected in self.expected_sequences}
                                - set(sources))
        if unknown_checks:
            raise ValueError(
                f"Expected sequence sources have no summary: {', '.join(unknown_checks)}"
            )
        object.__setattr__(self, "source_files", MappingProxyType(sources))
        object.__setattr__(self, "contexts", tuple(self.contexts))
        object.__setattr__(self, "expected_sequences", tuple(self.expected_sequences))
        object.__setattr__(self, "aliases", tuple(self.aliases))

    def source(self, name=None):
        """One of this collection's sources: its primary one by default."""
        name = name or self.primary_source
        if name not in self.sources():
            raise ValueError(
                f"Source {name!r} is not part of {self.name!r}; available: "
                f"{', '.join(self.sources())}"
            )
        return reference_source(name)

    def expected_for(self, source):
        """The sequences a rebuilt source summary must still reproduce."""
        return tuple(expected for expected in self.expected_sequences
                     if expected.source == source)

    def sources(self) -> tuple[str, ...]:
        """Every source this collection draws on, primary first."""
        return (self.primary_source, *(
            name for name in self.source_files if name != self.primary_source
        ))

    def inventory(self) -> pd.DataFrame:
        """The published inventory as a data frame."""
        return self.inventory_file.read()

    def source_summary(self, source: str | None = None) -> pd.DataFrame:
        """What each entry looks like at one source, keyed by abbreviation.

        The primary source is used when ``source`` is omitted, just as it is
        for :meth:`source`. A source name not declared by this collection is
        an error rather than an empty frame that could hide a typo.
        """
        source = source or self.primary_source
        if source not in self.sources():
            raise ValueError(
                f"Source {source!r} is not part of {self.name!r}; available: "
                f"{', '.join(self.sources())}"
            )
        return read_ledger(self.source_files[source])

    def to_dict(self) -> dict[str, object]:
        """A JSON-ready description, as the other ORFdata status objects give."""
        return {
            "name": self.name,
            "title": self.title,
            "target_species": self.target_species,
            "target_taxon_id": self.target_taxon_id,
            "panel": self.panel,
            "panel_url": self.panel_url,
            "panel_sha256": self.panel_sha256,
            "panel_size": self.panel_size,
            "panel_filename": self.panel_filename,
            "contexts": list(self.contexts),
            "references": self.references,
            "primary_source": self.primary_source,
            "alternate_source": self.alternate_source,
            "alternate_note": self.alternate_note,
            "inventory_file": self.inventory_file.to_dict(),
            "source_files": {name: table.to_dict()
                             for name, table in sorted(self.source_files.items())},
            "expected_sequences": [expected.to_dict() for expected in self.expected_sequences],
            "aliases": list(self.aliases),
        }


GDC_VIRUSES = ViralCollection(
    name="gdc-viruses",
    title="GDC GRCh38.d1.vd1 viral screen",
    target_species="homo_sapiens",
    target_taxon_id=9606,
    panel="GRCh38.d1.vd1",
    panel_url="https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt",
    panel_sha256="f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56",
    panel_size=13398,
    panel_filename="GRCh83.d1.vd1_virus_decoy.txt",
    contexts=("cancer-transcriptome-screening",),
    primary_source=GENBANK_SOURCE,
    inventory_file=CheckedTSV(
        resource="data/viral/gdc-grch38-d1-vd1-coverage.tsv",
        sha256="0bfa5ffc184eafcdb615b00cbe79acdd6e5e4f22ff2f919a2c31d62eb4cc1492",
        schema=INVENTORY_SCHEMA,
        columns=INVENTORY_COLUMNS,
        numeric=INVENTORY_NUMERIC,
    ),
    references=200,
    source_files={
        GENBANK_SOURCE: ledger_table("data/viral/gdc-grch38-d1-vd1-genbank.tsv",
                                    "14a8162627e3696e734520e7a0b15541d03185c3ee5652e86f9a5ad7545b9ebe"),
        PAVE_SOURCE: ledger_table("data/viral/gdc-grch38-d1-vd1-pave.tsv",
                                 "331db58d91663aea8f0d49e8e101bb7e150fb4ced5ce16b35f889b85ee81e515"),
    },
    alternate_source=PAVE_SOURCE,
    expected_sequences=(ExpectedSequence(
        source=PAVE_SOURCE,
        abbreviation="HPV16",
        # Pinned in tools/build_hpv16_cds.py from GDC's published reference files.
        sequence_md5="f2463721818c29998ef5ab917e6df74a",
        note="the GDC decoy contig is the alternate source's molecule, which is why "
             "this collection follows it",
    ),),
    alternate_note=(
        "GDC's decoy FASTA carries the PaVE reference where it differs from the GenBank "
        "accession beside it: PaVE's HPV16REF reproduces the GDC HPV16 DNA exactly, while "
        "GenBank's K02718.1 is a different, shorter molecule."
    ),
)

_COLLECTIONS = {GDC_VIRUSES.name: GDC_VIRUSES}


def _normalized(value):
    """Normalize a human-facing collection name without changing its meaning."""
    return str(value).strip().lower().replace("_", "-")


def _taxon_id(species):
    """Resolve a species selector to a positive NCBI taxon."""
    from .reference import species_taxon

    if species is None:
        return None
    if type(species) is int:
        if species > 0:
            return species
    else:
        found = species_taxon(species)
        if found is not None:
            return found
    raise ValueError(f"Unknown species alias {species!r}")


def viral_collections(
    *, species: str | int | None = None, source: str | None = None
) -> tuple[ViralCollection, ...]:
    """Published collections, optionally limited to one species or source.

    This is discovery only: it reads no table and downloads nothing. ``species``
    accepts the same aliases as the rest of ORFdata, and ``source`` is one of
    the names returned by :meth:`ViralCollection.sources`.
    """
    taxon_id = _taxon_id(species)
    found = tuple(_COLLECTIONS[name] for name in sorted(_COLLECTIONS))
    if taxon_id is not None:
        found = tuple(item for item in found if item.target_taxon_id == taxon_id)
    if source is not None:
        source = getattr(source, "name", source)
        reference_source(source)  # A typo is an error, not an empty listing.
        found = tuple(item for item in found if source in item.sources())
    return found


def viral_collection(
    name: ViralCollection | str | None = None, *, species: str | int | None = None
) -> ViralCollection:
    """Resolve one collection by object, name, alias, panel, or target species.

    With no name, the only collection for ``species`` is selected. With
    neither argument, the package's sole default collection is returned.
    Ambiguity is always reported instead of resolved by registry order.
    """
    from .reference import canonical_species

    taxon_id = _taxon_id(species)

    if isinstance(name, ViralCollection):
        collection = name
    elif name is None:
        if taxon_id is None:
            collection = GDC_VIRUSES
        else:
            candidates = viral_collections(species=taxon_id)
            if not candidates:
                raise ValueError(
                    f"No viral collection is registered for {species!r}; "
                    f"the default {GDC_VIRUSES.name!r} is a different target"
                )
            if len(candidates) != 1:
                names = ", ".join(item.name for item in candidates) or "none"
                raise ValueError(
                    f"Expected one viral collection for {species!r}; available: {names}"
                )
            collection, = candidates
    else:
        normalized = _normalized(name)
        matches = [item for item in viral_collections()
                   if normalized in {
                       _normalized(item.name),
                       _normalized(item.panel),
                       *(_normalized(alias) for alias in item.aliases),
                   }]
        if len(matches) != 1:
            available = ", ".join(item.name for item in viral_collections()) or "none"
            raise ValueError(f"Unknown viral collection {name!r}; available: {available}")
        collection, = matches

    if taxon_id is not None and taxon_id != collection.target_taxon_id:
        raise ValueError(
            f"The {collection.name} collection is for target species "
            f"{collection.target_species} (NCBI taxon {collection.target_taxon_id}); "
            f"{species!r} is a different target."
        )
    canonical_species(collection.target_species)  # Fails loudly on an unknown species name.
    return collection
