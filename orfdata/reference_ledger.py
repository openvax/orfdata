"""What each reference entry looks like at one source.

A reference inventory often names an entry more than once — an NCBI accession
and a specialist database's identifier — and those are not always the same
molecule. Which one a collection should carry is decided by what the people
who built the reference actually used, not by which identifier is listed
first.

A ledger records, per entry and per source, the identifier that source
publishes it under, the length and checksums of the molecule, how many coding
features it annotates, and the description it gives. Every source is written
down the same way, so the same ledger answers two questions offline: whether a
downloaded record is the molecule the collection expects, and what a second
source would add. A row the source could not serve records why, and pins
nothing.

Ledgers are packaged, checksummed :class:`~orfdata.reference_table.PackagedTable`
data, built by ``tools/build_reference_ledger.py``.
"""

import hashlib
from dataclasses import dataclass

from .reference_table import PackagedTable, parse_table, write_table

SCHEMA = "orfdata.reference_ledger.v2"
COLUMNS = ("index", "abbreviation", "source", "identifier", "state", "length_nt",
           "sequence_md5", "sequence_sha256", "coding_features", "taxon_id", "title")
# What the source turned out to be for an entry.
PUBLISHED = "published"        # It publishes the entry and the record checks out.
UNAVAILABLE = "unavailable"    # It does not publish the entry at all.
INCONSISTENT = "inconsistent"  # It publishes one whose own annotation does not check out.
STATES = (PUBLISHED, UNAVAILABLE, INCONSISTENT)


def ledger_table(resource, sha256):
    """The packaged table one ledger is read from."""
    return PackagedTable(resource=resource, sha256=sha256, schema=SCHEMA, columns=COLUMNS)


@dataclass(frozen=True)
class LedgerAnchor:
    """One molecule a ledger must reproduce for the claim resting on it to hold.

    A ledger is only worth pinning if the source it describes is the source the
    reference was really built from. An anchor names the one entry that settles
    that, so rebuilding the ledger re-checks the claim instead of quietly
    recording a different answer.
    """

    source: str
    abbreviation: str
    sequence_md5: str
    note: str = ""

    def check(self, entries):
        """Raise unless a freshly built ledger still reproduces this molecule."""
        entry = entries.get(self.abbreviation)
        if entry is None or entry.sequence_md5 != self.sequence_md5:
            raise ValueError(
                f"{self.source} {self.abbreviation} no longer reproduces "
                f"{self.sequence_md5}{': ' + self.note if self.note else ''}")
        return entry


@dataclass(frozen=True)
class LedgerEntry:
    """One entry's reference molecule at one source."""

    index: int
    abbreviation: str
    source: str
    identifier: str
    state: str
    length: int = 0
    sequence_md5: str = ""
    sequence_sha256: str = ""
    coding_features: int = 0
    taxon_id: int = 0
    title: str = ""

    @property
    def published(self):
        """Whether the source publishes an entry that checked out."""
        return self.state == PUBLISHED

    def matches(self, sequence):
        """Whether a molecule is the one this entry records."""
        if not self.published:
            return False
        return hashlib.md5(sequence.upper().encode("ascii")).hexdigest() == self.sequence_md5

    def verify(self, record):
        """Check a freshly fetched record against what the ledger pins."""
        digest, sha256 = record.checksums()
        if (record.length != self.length or digest != self.sequence_md5
                or sha256 != self.sequence_sha256):
            raise ValueError(
                f"{self.source} {self.identifier} no longer serves the molecule the packaged "
                f"ledger pins ({record.length} bases against {self.length})"
            )
        if len(record.coding_features) != self.coding_features:
            raise ValueError(
                f"{self.source} {self.identifier} serves {len(record.coding_features)} coding "
                f"features, not the {self.coding_features} the packaged ledger pins"
            )
        return record

    def row(self):
        """This entry as a ledger row; a row the source could not serve pins nothing."""
        measured = {"length_nt": self.length, "sequence_md5": self.sequence_md5,
                    "sequence_sha256": self.sequence_sha256,
                    "coding_features": self.coding_features, "taxon_id": self.taxon_id or ""}
        return {"index": self.index, "abbreviation": self.abbreviation, "source": self.source,
                "identifier": self.identifier, "state": self.state, "title": self.title,
                **({key: "" for key in measured} if not self.published else measured)}


def entries_from_rows(rows, name):
    """Build ledger entries from parsed rows, keyed by entry abbreviation."""
    entries = {}
    for row in rows:
        entry = LedgerEntry(
            index=int(row["index"]), abbreviation=row["abbreviation"], source=row["source"],
            identifier=row["identifier"], state=row["state"],
            length=int(row["length_nt"] or 0), sequence_md5=row["sequence_md5"],
            sequence_sha256=row["sequence_sha256"],
            coding_features=int(row["coding_features"] or 0),
            taxon_id=int(row["taxon_id"] or 0), title=row["title"])
        if entry.state not in STATES:
            raise ValueError(f"Unknown ledger state {entry.state!r} in {name}")
        if entry.abbreviation in entries:
            raise ValueError(f"Repeated ledger entry {entry.abbreviation} in {name}")
        entries[entry.abbreviation] = entry
    return entries


def read_ledger(table):
    """Read a packaged ledger, keyed by entry abbreviation."""
    return entries_from_rows(table.read(), table.name)


def read_ledger_bytes(raw, *, name="reference ledger"):
    """Read a ledger that is not packaged yet, as its builder does to check its own work."""
    return entries_from_rows(parse_table(raw, schema=SCHEMA, columns=COLUMNS, name=name)[1], name)


def write_ledger(entries, *, panel, source, source_url, checked, notes=()):
    """Serialize a ledger, as its builder tool writes it."""
    metadata = {"panel": panel, "source": f"{source} ({source_url})", "checked_utc": checked}
    metadata.update(dict(notes))
    return write_table((entry.row() for entry in entries), schema=SCHEMA, columns=COLUMNS,
                       metadata=metadata)
