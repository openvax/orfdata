"""One way to read and write the reference tables ORFdata packages.

ORFdata ships small tables beside its code: the inventory a collection
follows, and what each of its sources publishes for every entry in it. They
are all the same kind of thing — a tab-separated table with a named schema,
pinned by SHA-256 so a packaged file is read only when it is the file that was
built — so they are read and written here once rather than once per table.

A table reads as a :class:`pandas.DataFrame`, with the builder's own
``# key: value`` lines above the header kept in ``frame.attrs``: where the
data came from, when it was checked, and whatever else belongs on the record.
Another source, another panel or another species is a :class:`PackagedTable`
on a registry entry, not another loader.
"""

import hashlib
import io
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

import pandas as pd

SCHEMA_KEY = "schema"


@dataclass(frozen=True)
class PackagedTable:
    """A checksummed table packaged with ORFdata, named by its schema and columns."""

    #: Path of the file inside the ``orfdata`` package.
    resource: str
    #: SHA-256 the packaged bytes must have before a single row is read.
    sha256: str
    #: Schema name the file's own ``# schema:`` line must carry.
    schema: str
    #: Every column the file must have, in order; this is the frame's shape.
    columns: tuple[str, ...]
    #: Which of those columns hold whole numbers, read as nullable ``Int64``
    #: so a blank stays a blank instead of becoming a float or a string.
    numeric: tuple[str, ...] = ()

    @property
    def name(self):
        """The file name this table is packaged under."""
        return self.resource.rsplit("/", 1)[-1]

    def read(self):
        """This table, with its builder's header lines in ``attrs``.

        The frame's columns are exactly :attr:`columns`, in that order. Every
        column reads as ``str`` except those named in :attr:`numeric`, which
        read as nullable ``Int64`` so a blank stays blank. The packaged bytes
        are read and checked once; every caller gets its own frame to do as it
        likes with.
        """
        return _packaged(self).copy()

    def to_dict(self):
        """A JSON-ready description, as the other ORFdata status objects give."""
        return {"resource": self.resource, "sha256": self.sha256, "schema": self.schema,
                "columns": list(self.columns)}


def packaged_bytes(resource, sha256):
    """One packaged file's bytes, read only when they are the bytes that were pinned."""
    raw = files("orfdata").joinpath(resource).read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha256:
        raise ValueError(f"Packaged data checksum mismatch: {resource}")
    return raw


@lru_cache(maxsize=16)
def _packaged(table):
    return read_table(packaged_bytes(table.resource, table.sha256), schema=table.schema,
                      columns=table.columns, numeric=table.numeric, name=table.name)


def read_table(raw, *, schema, columns, numeric=(), name="reference table"):
    """Read a reference table's bytes, checking what it says it is before trusting it.

    Returns a frame whose columns are exactly ``columns``, those in ``numeric``
    as ``Int64`` and the rest as ``str``, with the file's ``# key: value``
    header lines in ``frame.attrs``.
    """
    text = raw.decode("utf-8")
    metadata = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            break
        key, _, value = line[1:].partition(":")
        metadata[key.strip()] = value.strip()
    if metadata.get(SCHEMA_KEY) != schema:
        raise ValueError(f"Unsupported {name} schema {metadata.get(SCHEMA_KEY)!r}; "
                         f"this ORFdata reads {schema!r}")
    body = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    # Everything arrives as text, so a blank is a blank and nothing is guessed
    # at; the declared numeric columns are then converted deliberately.
    frame = pd.read_csv(io.StringIO(body), sep="\t", dtype="str", keep_default_na=False)
    if tuple(frame.columns) != tuple(columns):
        raise ValueError(f"Unsupported {name} columns: {list(frame.columns)}")
    for column in numeric:
        try:
            frame[column] = pd.to_numeric(frame[column].replace("", None)).astype("Int64")
        except (ValueError, TypeError) as error:
            raise ValueError(f"Unreadable {name} column {column!r}: {error}") from error
    frame.attrs.update(metadata)
    return frame


def write_table(frame, *, schema, columns=None, metadata=()):
    """Serialize a reference table, as its builder tool writes it."""
    if columns is not None:
        frame = frame[list(columns)]
    lines = [f"# {SCHEMA_KEY}: {schema}"]
    lines.extend(f"# {key}: {value}" for key, value in dict(metadata).items())
    body = frame.to_csv(sep="\t", index=False, na_rep="", lineterminator="\n")
    return ("\n".join(lines) + "\n" + body).encode("utf-8")
