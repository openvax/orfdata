"""The viral reference collections ORFdata publishes, as data rather than constants.

A collection says what it is for and where it comes from: the species whose
analysis it supports, the published inventory it follows, the sources its
molecules can be fetched from, and which source to prefer where the inventory
names two references that are not the same molecule.

Everything a collection needs offline is a packaged, checksummed table: the
inventory itself, and one ledger per source recording what that source
publishes for every entry in it. They all have the same shape, so a second
collection — another species, another inventory, another pair of sources — is
a new entry here, not a new code path.

One collection is published today, ``gdc-viruses``, for human analysis.
"""

from dataclasses import dataclass, field
from typing import Optional

from .record_sources import GENBANK_SOURCE, PAVE_SOURCE, reference_source
from .reference_ledger import LedgerAnchor, ledger_table, read_ledger
from .reference_table import PackagedTable

INVENTORY_SCHEMA = "orfdata.viral_panel_coverage.v1"
INVENTORY_COLUMNS = ("index", "name", "abbreviation", "refseq", "genbank", "pave_id",
                     "pave_notes", "bundled_reference_id", "coding_models", "status")


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
    inventory: PackagedTable
    references: int
    # Out of the comparison so a collection stays hashable as a cache key; two
    # collections with the same inventory are the same inventory either way.
    ledgers: dict = field(default_factory=dict, repr=False, compare=False)
    alternate_source: Optional[str] = None
    alternate_note: str = ""
    anchors: tuple[LedgerAnchor, ...] = ()
    aliases: tuple[str, ...] = field(default=(), repr=False)

    def source(self, name=None):
        """One of this collection's sources: its primary one by default."""
        return reference_source(name or self.primary_source)

    def source_anchors(self, source):
        """The molecules a rebuilt ledger for one source must still reproduce."""
        return tuple(anchor for anchor in self.anchors if anchor.source == source)

    def sources(self):
        """Every source this collection draws on, primary first."""
        return tuple(name for name in (self.primary_source, self.alternate_source) if name)

    def ledger(self, source=None):
        """What each entry looks like at one source; empty when none is packaged.

        The alternate source is the default, because that is the comparison a
        collection is usually asking for.
        """
        table = self.ledgers.get(source or self.alternate_source or self.primary_source)
        return read_ledger(table) if table is not None else {}

    def to_dict(self):
        """A JSON-ready description, as the other ORFdata status objects give."""
        return {
            "name": self.name,
            "title": self.title,
            "target_species": self.target_species,
            "target_taxon_id": self.target_taxon_id,
            "panel": self.panel,
            "contexts": list(self.contexts),
            "references": self.references,
            "primary_source": self.primary_source,
            "alternate_source": self.alternate_source,
            "alternate_note": self.alternate_note,
            "inventory": self.inventory.to_dict(),
            "ledgers": {name: table.to_dict() for name, table in sorted(self.ledgers.items())},
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
    inventory=PackagedTable(
        resource="data/viral/gdc-grch38-d1-vd1-coverage.tsv",
        sha256="0bfa5ffc184eafcdb615b00cbe79acdd6e5e4f22ff2f919a2c31d62eb4cc1492",
        schema=INVENTORY_SCHEMA,
        columns=INVENTORY_COLUMNS,
    ),
    references=200,
    ledgers={
        GENBANK_SOURCE: ledger_table("data/viral/gdc-grch38-d1-vd1-genbank.tsv",
                                    "14a8162627e3696e734520e7a0b15541d03185c3ee5652e86f9a5ad7545b9ebe"),
        PAVE_SOURCE: ledger_table("data/viral/gdc-grch38-d1-vd1-pave.tsv",
                                 "ada474a021afdaf48b89ec417c3ceeb9d8812cf2a68e5fd4141c52b4d324cf9a"),
    },
    alternate_source=PAVE_SOURCE,
    anchors=(LedgerAnchor(
        source=PAVE_SOURCE, abbreviation="HPV16",
        # Pinned in tools/build_hpv16_cds.py from GDC's published reference files.
        sequence_md5="f2463721818c29998ef5ab917e6df74a",
        note="the GDC decoy contig is the alternate source's molecule, which is why "
             "this collection follows it"),),
    alternate_note=(
        "GDC's decoy FASTA carries the PaVE reference where it differs from the GenBank "
        "accession beside it: PaVE's HPV16REF reproduces the GDC HPV16 DNA exactly, while "
        "GenBank's K02718.1 is a different, shorter molecule."
    ),
)

_COLLECTIONS = {GDC_VIRUSES.name: GDC_VIRUSES}


def viral_collections():
    """Every published viral reference collection."""
    return tuple(_COLLECTIONS[name] for name in sorted(_COLLECTIONS))


def viral_collection(name, *, species=None):
    """One published collection, checked against the species it is for."""
    from .reference import canonical_species, species_taxon

    normalized = str(name).strip().lower().replace("_", "-")
    collection = _COLLECTIONS.get(normalized)
    if collection is None:
        raise ValueError(
            f"Unknown viral collection {name!r}; available: "
            f"{', '.join(sorted(_COLLECTIONS))}"
        )
    if species is not None:
        taxon_id = species if type(species) is int else species_taxon(species)
        if taxon_id != collection.target_taxon_id:
            raise ValueError(
                f"The {collection.name} collection is for target species "
                f"{collection.target_species} (NCBI taxon {collection.target_taxon_id}); "
                f"{species!r} is a different target."
            )
    canonical_species(collection.target_species)  # Fails loudly on an unknown species name.
    return collection
