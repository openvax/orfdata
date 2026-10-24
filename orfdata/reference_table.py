"""One way to read and write the reference tables ORFdata packages.

ORFdata ships small tables beside its code: the inventory a collection
follows, and what each of its sources publishes for every entry in it. They
are all the same kind of thing — a tab-separated table with a named schema,
pinned by SHA-256 so a packaged file is read only when it is the file that was
built — so they are read and written here once rather than once per table.

A table carries its builder's own ``# key: value`` lines above the header:
where the data came from, when it was checked, and whatever else belongs on
the record. Another source, another panel or another species is a
:class:`PackagedTable` on a registry entry, not another loader.
"""

import csv
import hashlib
import io
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

SCHEMA_KEY = "schema"


@dataclass(frozen=True)
class PackagedTable:
    """A checksummed table packaged with ORFdata, named by its schema and columns."""

    resource: str
    sha256: str
    schema: str
    columns: tuple[str, ...]

    @property
    def name(self):
        """The file name this table is packaged under."""
        return self.resource.rsplit("/", 1)[-1]

    def read(self):
        """Every row, as a dict keyed by column name."""
        return _packaged(self)[1]

    def metadata(self):
        """The builder's own header lines, as a mapping."""
        return _packaged(self)[0]

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
    return parse_table(packaged_bytes(table.resource, table.sha256), schema=table.schema,
                       columns=table.columns, name=table.name)


def parse_table(raw, *, schema, columns, name="reference table"):
    """Read a reference table's bytes, checking what it says it is before trusting it."""
    lines = raw.decode("utf-8").splitlines()
    metadata = {}
    for line in lines:
        if not line.startswith("#"):
            break
        key, _, value = line[1:].partition(":")
        metadata[key.strip()] = value.strip()
    if metadata.get(SCHEMA_KEY) != schema:
        raise ValueError(f"Unsupported {name} schema {metadata.get(SCHEMA_KEY)!r}; "
                         f"this ORFdata reads {schema!r}")
    body = io.StringIO("\n".join(line for line in lines if not line.startswith("#")))
    reader = csv.DictReader(body, delimiter="\t")
    if tuple(reader.fieldnames or ()) != tuple(columns):
        raise ValueError(f"Unsupported {name} columns: {reader.fieldnames}")
    return metadata, tuple(reader)


def write_table(rows, *, schema, columns, metadata=()):
    """Serialize a reference table, as its builder tool writes it."""
    lines = [f"# {SCHEMA_KEY}: {schema}"]
    lines.extend(f"# {key}: {value}" for key, value in dict(metadata).items())
    lines.append("\t".join(columns))
    lines.extend("\t".join(str(row[column]) for column in columns) for row in rows)
    return ("\n".join(lines) + "\n").encode("utf-8")
