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

Ledgers are packaged, checksummed :class:`~orfdata.reference_table.PackagedTable`
data, built by ``tools/build_reference_ledger.py``, and read as data frames
keyed by entry abbreviation.
"""

import hashlib
from dataclasses import dataclass

import pandas as pd

from .reference_table import PackagedTable, read_table, write_table

SCHEMA = "orfdata.reference_ledger.v2"
COLUMNS = ("index", "abbreviation", "source", "identifier", "state", "length_nt",
           "sequence_md5", "sequence_sha256", "coding_features", "taxon_id", "title")
NUMERIC = ("index", "length_nt", "coding_features", "taxon_id")
# What was measured off the molecule, and so is blank unless the source served one.
MEASURED = ("length_nt", "sequence_md5", "sequence_sha256", "coding_features", "taxon_id")
# What the source turned out to be for an entry.
PUBLISHED = "published"        # It publishes the entry and the record checks out.
UNAVAILABLE = "unavailable"    # It does not publish the entry at all.
INCONSISTENT = "inconsistent"  # It publishes one whose own annotation does not check out.
STATES = (PUBLISHED, UNAVAILABLE, INCONSISTENT)


def ledger_table(resource, sha256):
    """The packaged table one ledger is read from."""
    return PackagedTable(resource=resource, sha256=sha256, schema=SCHEMA, columns=COLUMNS,
                         numeric=NUMERIC)


@dataclass(frozen=True)
class Anchor:
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

    def check(self, ledger):
        """Raise unless a freshly built ledger still reproduces this molecule."""
        if (self.abbreviation not in ledger.index
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
    return frame.set_index("abbreviation")


def read_ledger(table):
    """Read a packaged ledger, keyed by entry abbreviation."""
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
    """A ledger with no rows, for a source that has none packaged."""
    frame = pd.DataFrame({column: pd.Series(dtype="Int64" if column in NUMERIC else "str")
                          for column in COLUMNS})
    return frame.set_index("abbreviation")


def ledger_from_rows(rows, name="reference ledger"):
    """A ledger built from freshly fetched rows, keyed by entry abbreviation."""
    frame = pd.DataFrame(list(rows), columns=list(COLUMNS))
    for column in COLUMNS:
        frame[column] = (frame[column].astype("Int64") if column in NUMERIC
                         else frame[column].fillna("").astype("str"))
    return _keyed(frame, name)
