"""What each reference entry looks like at one source.

A reference inventory often names an entry more than once — an NCBI accession
and a specialist database's identifier — and those are not always the same
molecule. Which one a collection should carry is decided by what the people
who built the reference actually used, not by which identifier is listed
first.

A ledger is a table, one row per entry, recording the identifier its source
publishes the entry under, the length and checksums of the molecule, how many
coding features it annotates, and the description it gives. Every source is
written down the same way, so the same ledger answers two questions offline:
whether a downloaded record is the molecule the collection expects, and what a
second source would add. A row the source could not serve records why, and
pins nothing.

These source summaries are checked :class:`~orfdata.reference_table.CheckedTSV` data shipped
inside ORFdata and built by ``tools/build_reference_ledger.py``.

A ledger reads as a data frame indexed by ``abbreviation`` (:data:`INDEX`)
whose columns are :data:`FRAME_COLUMNS`:

===================  =======  =======================================================
column               dtype    meaning
===================  =======  =======================================================
``abbreviation``     str      the inventory row this describes; the frame's index
``index``            Int64    that row's 1-based position in the inventory
``source``           str      which source this row speaks for, e.g. ``genbank``
``identifier``       str      what the source publishes the entry under, versioned
``state``            str      one of :data:`STATES`; only ``published`` pins anything
``length_nt``        Int64    bases in the molecule the source serves
``sequence_md5``     str      MD5 of the uppercased sequence, as GDC publishes it
``sequence_sha256``  str      SHA-256 of the same bases
``coding_features``  Int64    how many CDS features the source annotates on it
``taxon_id``         Int64    the NCBI taxon the record cross-references, if any
``title``            str      the source's own one-line description of the record
===================  =======  =======================================================

Every column but ``abbreviation``, ``index``, ``source``, ``identifier``,
``state`` and ``title`` is :data:`MEASURED` off the molecule, so it is blank
(``pd.NA`` or ``""``) wherever ``state`` is not ``published``. ``taxon_id`` is
also blank where the source names no taxon. :data:`COLUMNS` is the on-disk
column order, which carries ``abbreviation`` as an ordinary column.
"""

import hashlib
from dataclasses import dataclass

import pandas as pd

from .reference_table import CheckedTSV, read_table, write_table

SCHEMA = "orfdata.reference_ledger.v2"
COLUMNS = ("index", "abbreviation", "source", "identifier", "state", "length_nt",
           "sequence_md5", "sequence_sha256", "coding_features", "taxon_id", "title")
NUMERIC = ("index", "length_nt", "coding_features", "taxon_id")
# A read ledger is keyed by entry abbreviation, so that column is its index.
INDEX = "abbreviation"
FRAME_COLUMNS = tuple(column for column in COLUMNS if column != INDEX)
# What was measured off the molecule, and so is blank unless the source served one.
MEASURED = ("length_nt", "sequence_md5", "sequence_sha256", "coding_features", "taxon_id")
# What the source turned out to be for an entry.
PUBLISHED = "published"        # It publishes the entry and the record checks out.
UNAVAILABLE = "unavailable"    # It does not publish the entry at all.
INCONSISTENT = "inconsistent"  # It publishes one whose own annotation does not check out.
STATES = (PUBLISHED, UNAVAILABLE, INCONSISTENT)


def ledger_table(resource, sha256):
    """The checked TSV file one ledger is read from."""
    return CheckedTSV(resource=resource, sha256=sha256, schema=SCHEMA, columns=COLUMNS,
                      numeric=NUMERIC)


@dataclass(frozen=True)
class ExpectedSequence:
    """One sequence a rebuilt ledger must still contain.

    A ledger is only worth pinning if the source it describes is the source the
    reference was really built from. This names the one entry that settles
    that, so rebuilding the ledger re-checks the claim instead of quietly
    recording a different answer.
    """

    source: str
    abbreviation: str
    sequence_md5: str
    note: str = ""

    def to_dict(self):
        """Return the four parts of this claim as plain data."""
        return {
            "source": self.source,
            "abbreviation": self.abbreviation,
            "sequence_md5": self.sequence_md5,
            "note": self.note,
        }

    def check(self, ledger):
        """Raise unless a freshly built ledger still reproduces this molecule."""
        if (self.abbreviation not in ledger.index
                or ledger.at[self.abbreviation, "source"] != self.source
                or ledger.at[self.abbreviation, "sequence_md5"] != self.sequence_md5):
            raise ValueError(
                f"{self.source} {self.abbreviation} no longer reproduces "
                f"{self.sequence_md5}{': ' + self.note if self.note else ''}")
        return ledger.loc[self.abbreviation]


def ledger_row(ledger, abbreviation):
    """One entry's row, or ``None`` where the source has nothing for it."""
    return ledger.loc[abbreviation] if abbreviation in ledger.index else None


def published(row):
    """Whether a ledger row's source publishes an entry that checked out."""
    return row is not None and row["state"] == PUBLISHED


def same_molecule(row, sequence):
    """Whether a molecule is the one a ledger row records."""
    if row["state"] != PUBLISHED:
        return False
    return hashlib.md5(sequence.upper().encode("ascii")).hexdigest() == row["sequence_md5"]


def verify_record(row, record):
    """Check a freshly fetched record against what a ledger row pins."""
    if row["state"] != PUBLISHED:
        raise ValueError(f"{row['source']} does not publish {record.identifier}")
    if record.source != row["source"] or record.identifier != row["identifier"]:
        raise ValueError(
            f"Expected {row['source']} {row['identifier']}, got "
            f"{record.source} {record.identifier}"
        )
    if (record.length != row["length_nt"] or record.checksums() != (row["sequence_md5"],
                                                                    row["sequence_sha256"])):
        raise ValueError(
            f"{row['source']} {row['identifier']} no longer serves the molecule the packaged "
            f"ledger pins ({record.length} bases against {row['length_nt']})"
        )
    if len(record.coding_features) != row["coding_features"]:
        raise ValueError(
            f"{row['source']} {row['identifier']} serves {len(record.coding_features)} coding "
            f"features, not the {row['coding_features']} the packaged ledger pins"
        )
    return record


def _keyed(frame, name):
    """Check a ledger says something sensible, and key it by entry abbreviation."""
    unknown = sorted(set(frame["state"]) - set(STATES))
    if unknown:
        raise ValueError(f"Unknown ledger state {unknown[0]!r} in {name}")
    repeated = frame["abbreviation"][frame["abbreviation"].duplicated()]
    if len(repeated):
        raise ValueError(f"Repeated ledger entry {repeated.iloc[0]} in {name}")
    return frame.set_index(INDEX)


def read_ledger(table):
    """Read a packaged ledger: a frame indexed by :data:`INDEX`, columns :data:`FRAME_COLUMNS`."""
    return _keyed(table.read(), table.name)


def read_ledger_bytes(raw, *, name="reference ledger"):
    """Read a ledger that is not packaged yet, as its builder does to check its own work."""
    return _keyed(read_table(raw, schema=SCHEMA, columns=COLUMNS, numeric=NUMERIC, name=name),
                  name)


def write_ledger(ledger, *, panel, source, source_url, checked, notes=()):
    """Serialize a ledger, as its builder tool writes it."""
    frame = ledger.reset_index() if ledger.index.name == "abbreviation" else ledger.copy()
    # A taxon nobody told us is a blank rather than organism zero, and a row
    # the source could not serve pins nothing, whatever was measured. The
    # taxon goes first so this never has to mask over blanks it just wrote.
    frame.loc[frame["taxon_id"] == 0, "taxon_id"] = pd.NA
    unserved = frame["state"] != PUBLISHED
    frame.loc[unserved, ["length_nt", "coding_features", "taxon_id"]] = pd.NA
    frame.loc[unserved, ["sequence_md5", "sequence_sha256"]] = ""
    metadata = {"panel": panel, "source": f"{source} ({source_url})", "checked_utc": checked}
    metadata.update(dict(notes))
    return write_table(frame, schema=SCHEMA, columns=COLUMNS, metadata=metadata)


def empty_ledger():
    """A ledger with no rows, for a source that has none packaged.

    It carries the same index and columns as any other ledger, so a caller
    never has to tell "nothing packaged" apart from "nothing published".
    """
    frame = pd.DataFrame({column: pd.Series(dtype="Int64" if column in NUMERIC else "str")
                          for column in COLUMNS})
    return frame.set_index(INDEX)


def ledger_from_rows(rows, name="reference ledger"):
    """A ledger built from freshly fetched rows, in the shape :data:`FRAME_COLUMNS` describes."""
    frame = pd.DataFrame(list(rows), columns=list(COLUMNS))
    for column in COLUMNS:
        frame[column] = (frame[column].astype("Int64") if column in NUMERIC
                         else frame[column].fillna("").astype("str"))
    return _keyed(frame, name)
