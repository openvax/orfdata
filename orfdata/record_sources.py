"""Where a reference molecule and its annotation come from.

A reference set names each of its entries in more than one place — an NCBI
accession, a specialist database's reference identifier — and those are not
always the same molecule. This module makes a source an explicit, named thing
that can be asked for one entry's record, so a collection can say which source
each of its genomes came from and store that source's own document beside it.

Two sources are published:

``genbank``
    NCBI's nucleotide records, fetched through E-utilities. Used for any
    reference named by accession.
``pave``
    The Papillomavirus Episteme's reference genomes, which carry spliced
    transcript models (E1^E4, E8^E2, E6*) that the GenBank submissions behind
    the same papillomavirus types usually lack.

Adding a source means adding a :class:`ReferenceSource` with a name, the
inventory column that names an entry at it, a way to fetch one record, and a
way to read that record's coding features. Everything else — storing the
record, building a ledger of what the source publishes, checking a download
against it — is the same for every source. Nothing else in ORFdata needs to
know which source a reference came from, only what the record says.
"""

import gzip
import hashlib
import json
import time
from dataclasses import dataclass, field

from .net import SourceUnavailable, fetch
from .reference_ledger import INCONSISTENT, PUBLISHED, UNAVAILABLE, ledger_from_rows
from .viral_annotation import (
    agrees_with,
    coding_feature,
    extract,
    parse_genbank,
    safe_token,
    translate,
)

NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAVE_API_URL = "https://pave.niaid.nih.gov/api"
GENBANK_SOURCE = "genbank"
PAVE_SOURCE = "pave"
_PAVE_TRANSLATION_TABLE = 1
# Both sources ask that an unauthenticated client stay under three requests a second.
REQUEST_INTERVAL = 0.34


@dataclass(frozen=True)
class FetchedGenome:
    """One fetched genome, its coding features, and the source document verbatim."""

    source: str
    identifier: str
    sequence: str
    coding_features: tuple
    document: bytes = field(repr=False, default=b"")
    suffix: str = ""
    circular: bool = False
    title: str = ""
    taxon_id: int = 0

    @property
    def length(self):
        """Bases in this molecule."""
        return len(self.sequence)

    def checksums(self):
        """The molecule's MD5 and SHA-256, as reference ledgers record them."""
        encoded = self.sequence.encode("ascii")
        return hashlib.md5(encoded).hexdigest(), hashlib.sha256(encoded).hexdigest()

    def document_name(self, label):
        """The file name this record is stored under in a collection."""
        return f"{safe_token(label)}__{self.source}__{safe_token(self.identifier)}{self.suffix}.gz"

    def write_document(self, directory, label):
        """Store the source's own document, exactly as it was served."""
        path = directory / self.document_name(label)
        with gzip.GzipFile(path, "wb", mtime=0) as stream:
            stream.write(self.document)
        return path


class ReferenceSource:
    """A place that publishes reference molecules, asked one entry at a time."""

    name = ""
    description = ""
    url = ""
    # Which inventory column names an entry at this source.
    inventory_field = ""
    batch_size = 1

    def identifier(self, entry):
        """What this source publishes one inventory entry under, or ``""`` if it names none."""
        return getattr(entry, self.inventory_field, None) or ""

    def record(self, identifier, *, timeout):
        """The record this source publishes under ``identifier``."""
        raise NotImplementedError

    def records(self, identifiers, *, timeout):
        """Records for several identifiers; sources that batch override this."""
        return tuple(self.record(identifier, timeout=timeout) for identifier in identifiers)

    def build_ledger(self, inventory, *, timeout, report=None):
        """What this source publishes for every inventory row it names, in inventory order.

        One pass over an inventory produces one source's whole ledger, whatever
        the source is: an identifier it does not publish and a record whose own
        annotation does not check out are both answers to record, not failures
        to stop for.
        """
        named = [(entry, identifier) for entry in inventory.itertuples(index=False)
                 if (identifier := self.identifier(entry))]
        rows = []
        for start in range(0, len(named), self.batch_size):
            batch = named[start:start + self.batch_size]
            found = self._attempt(tuple(identifier for _, identifier in batch), timeout=timeout)
            rows.extend(_ledger_row(entry, identifier, self.name, result, report)
                        for (entry, identifier), result in zip(batch, found, strict=True))
            if start + self.batch_size < len(named):
                time.sleep(REQUEST_INTERVAL)
        return ledger_from_rows(rows)

    def _attempt(self, identifiers, *, timeout):
        """Each identifier's record, or the error standing in its place.

        A batch that fails as a whole is retried one at a time, so a single
        entry the source cannot serve does not take the rest of its batch down.
        """
        try:
            return self.records(identifiers, timeout=timeout)
        except (RuntimeError, ValueError, OSError) as error:
            if len(identifiers) == 1:
                return (error,)
            return tuple(self._attempt((identifier,), timeout=timeout)[0]
                         for identifier in identifiers)


def _ledger_row(entry, identifier, source, result, report=None):
    """One inventory row as a ledger row, whether or not the source could serve it."""
    row = {"index": entry.index, "abbreviation": entry.abbreviation, "source": source,
           "identifier": identifier, "title": ""}
    if isinstance(result, Exception):
        if report is not None:
            report(entry, result)
        # The source not publishing an identifier and the source publishing one
        # whose own annotation does not check out are different answers.
        return {**row, "state": UNAVAILABLE if isinstance(result, SourceUnavailable)
                else INCONSISTENT}
    digest, sha256 = result.checksums()
    return {**row, "identifier": result.identifier, "state": PUBLISHED, "title": result.title,
            "length_nt": result.length, "sequence_md5": digest, "sequence_sha256": sha256,
            "coding_features": len(result.coding_features), "taxon_id": result.taxon_id}


class GenBankSource(ReferenceSource):
    """NCBI nucleotide records, fetched in batches through E-utilities."""

    name = GENBANK_SOURCE
    description = "NCBI GenBank nucleotide records"
    url = NCBI_EFETCH_URL
    inventory_field = "genbank"
    batch_size = 25

    def __init__(self, *, tool="orfdata", email="support@openvax.org"):
        self.tool, self.email = tool, email

    def records(self, identifiers, *, timeout):
        identifiers = tuple(identifiers)
        payload = fetch(self.url, timeout=timeout, description="NCBI EFetch", params={
            "db": "nuccore",
            "id": ",".join(identifiers),
            "rettype": "gbwithparts",
            "retmode": "text",
            "tool": self.tool,
            "email": self.email,
        })
        parsed = parse_genbank(payload)
        if len(parsed) != len(identifiers):
            raise ValueError(f"NCBI returned {len(parsed)} records for {len(identifiers)} accessions")
        return tuple(FetchedGenome(
            source=self.name,
            identifier=record.accession_version,
            sequence=record.sequence,
            coding_features=record.coding_features,
            document=record.source_text.encode("utf-8"),
            suffix=".gb",
            circular=record.circular,
            title=record.definition,
            taxon_id=record.taxon_id,
        ) for record in parsed)

    def record(self, identifier, *, timeout):
        found, = self.records((identifier,), timeout=timeout)
        return found


class PaVESource(ReferenceSource):
    """Papillomavirus Episteme reference genomes and their gene annotations."""

    name = PAVE_SOURCE
    description = "Papillomavirus Episteme (PaVE) reference genomes"
    url = f"{PAVE_API_URL}/genome"
    inventory_field = "pave_id"

    def record(self, identifier, *, timeout):
        document = fetch(f"{self.url}/{identifier}", timeout=timeout,
                         description=f"PaVE genome {identifier}")
        sequence, features, title = self._read(document, identifier)
        return FetchedGenome(
            source=self.name,
            identifier=identifier,
            sequence=sequence,
            coding_features=self._coding_features(identifier, sequence, features),
            document=document,
            suffix=".json",
            title=title,
        )

    @staticmethod
    def _read(document, identifier):
        try:
            genome = json.loads(document)
            sequence = genome["itemSequence"]["sequence"].upper()
            features = genome["features"]
            title = str((genome.get("genome") or {}).get("description") or "")
            if not isinstance(features, list) or not all(isinstance(f, dict) for f in features):
                raise TypeError("features are not a list of objects")
        except (ValueError, KeyError, TypeError, AttributeError) as error:
            raise ValueError(f"Unreadable PaVE genome record for {identifier}: {error}") from error
        return sequence, features, title

    @staticmethod
    def _coding_features(identifier, sequence, features):
        """Convert PaVE's CDS features, checking each against its published protein."""
        coding = []
        for ordinal, feature in enumerate((item for item in features
                                           if item.get("featureType") == "CDS"), 1):
            try:
                exons = tuple((int(block["start"]), int(block["end"]), "+")
                              for block in feature["featureLocations"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"Unreadable PaVE feature locations in {identifier}: {error}") from error
            if any(end < start or start < 1 or end > len(sequence) for start, end, _ in exons):
                raise ValueError(f"PaVE feature {feature.get('featureID', ordinal)} of "
                                 f"{identifier} does not lie on the published sequence")
            nucleotides = extract(sequence, exons)
            protein = feature.get("proteinSequence") or ""
            if not agrees_with(translate(nucleotides, table=_PAVE_TRANSLATION_TABLE), protein):
                raise ValueError(
                    f"PaVE feature {feature.get('featureID', ordinal)} of {identifier} does not "
                    "translate to its own published protein"
                )
            coding.append(coding_feature(
                ordinal, exons, nucleotides,
                {"gene": feature.get("featureName") or "",
                 "product": feature.get("featureDescription") or "",
                 "protein_id": feature.get("proteinID") or "",
                 "translation": protein,
                 "transl_table": _PAVE_TRANSLATION_TABLE},
                annotation_source=PAVE_SOURCE))
        if not coding:
            raise ValueError(f"PaVE genome {identifier} annotates no CDS feature")
        return tuple(coding)


_SOURCES = {source.name: source for source in (GenBankSource(), PaVESource())}


def reference_sources():
    """Every published reference source, by name."""
    return dict(_SOURCES)


def reference_source(name):
    """One published reference source."""
    try:
        return _SOURCES[name]
    except KeyError:
        raise ValueError(
            f"Unknown reference source {name!r}; available: {', '.join(sorted(_SOURCES))}"
        ) from None


def read_document(path):
    """Read back a stored source document, whatever source wrote it."""
    raw = gzip.decompress(path.read_bytes())
    if path.name.endswith(".json.gz"):
        return json.loads(raw)
    return raw.decode("utf-8")
