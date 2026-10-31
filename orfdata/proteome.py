"""Reference proteomes derived from transcript models and selected protein sources.

A reference transcriptome supplies transcript models; a proteome selects
reference CDS products, catalogue ORFs and optionally viral proteins. Each layer keeps its
source and selection, every transcript a protein is encoded on is kept, and
every catalogue ORF that cannot be placed on the genome is reported rather
than dropped silently. A proteome can be saved as a checksummed snapshot and
reloaded offline.
"""

import gzip
import hashlib
import json
import re
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Optional
from urllib.parse import quote

from . import _structure
from .catalog import ORFCatalog, _positive_integer, load_catalog, merge_catalogs
from .errors import NotInstalledError
from .evidence import EVIDENCE_METHODS, SYNONYMS, EvidenceRule, evidence_rule
from .model import TranscriptORF, _contig, _id_forms, _layout_offset, _same_assembly, _versioned
from .reference import ReferenceContext
from .viral import ViralProduct, ViralReference, viral_reference, viral_references

SNAPSHOT_SCHEMA = "orfdata.proteome.v5"
_SNAPSHOT = "proteome.json.gz"
# PyEnsembl's convention: only ENS IDs carry a ".N" version; other dots are part of the ID.
_ENSEMBL_VERSIONED = re.compile(r"(ENS.*)\.(0|[1-9][0-9]*)")
# Evidence labels reported for an included ORF. Snapshots from before 0.28.0 keep "shotgun" and "mhc".
_EVIDENCE_LABELS = ("ribosome", "proteomics", "immunopeptidomics", "ms_unclassified")


def _evidence_labels(evidence):
    """Evidence kinds for a FASTA header; "ms" for MS support whose assay the source doesn't name."""
    labels = [m for m in _EVIDENCE_LABELS if evidence.has_evidence(m)]
    if evidence.has_evidence("ms") and not set(labels) - {"ribosome"}:
        labels.append("ms")
    return tuple(labels)


def _known(value):
    """Spell a three-state boolean without collapsing false into unknown."""
    return "unknown" if value is None else str(value).lower()


def _fasta_description_value(value):
    """Escape one FASTA description value so it cannot add fields or records."""
    return quote(str(value), safe=".,:@+-_")


@dataclass(frozen=True)
class BindingFailure:
    """A catalogue ORF that did not bind to a transcript, and why."""

    orf_id: str
    transcript_id: str
    reason: str


@dataclass(frozen=True)
class Occurrence:
    """One transcript encoding a protein, with the coding path's coordinates.

    ``start`` and ``end`` slice the transcript's spliced sequence; the genomic
    ``coding_blocks`` are 1-based inclusive. ``on_source_transcript`` says
    whether the protein's source defines it on this transcript.
    """

    transcript_id: str
    gene_id: str
    gene_name: Optional[str]
    transcript_biotype: Optional[str]
    contig: str
    strand: str
    start: int
    end: int
    coding_blocks: tuple
    on_source_transcript: bool

    @classmethod
    def from_orf(cls, orf):
        transcript = orf.transcript
        return cls(_versioned(transcript), transcript.gene_id, transcript.gene_name or None,
                   transcript.biotype, transcript.contig, transcript.strand, orf.start, orf.end,
                   tuple(tuple(b) for b in orf.coding_blocks), orf.on_source_transcript)


_OCCURRENCE_FIELDS = ("transcript_id", "gene_id", "gene_name", "transcript_biotype", "contig", "strand",
                      "start", "end", "coding_blocks", "on_source_transcript")


@dataclass(frozen=True)
class DerivedProtein:
    """A supplied sample-specific protein sequence with explicit parentage.

    The caller supplies amino acids; this record does not infer translation or
    establish expression. ``change_ids`` identify the variants or other events
    used to derive it, in the caller's own namespace.
    """

    accession: str
    sequence: str
    parent_accessions: tuple
    sample_id: str
    change_ids: tuple
    source: str
    source_sha256: Optional[str] = None
    parent_transcript_ids: tuple = ()

    def __post_init__(self):
        if any(isinstance(value, str) for value in
               (self.parent_accessions, self.change_ids, self.parent_transcript_ids)):
            raise TypeError("Parent accessions, transcripts and change IDs must be collections")
        object.__setattr__(self, "parent_accessions", tuple(self.parent_accessions))
        object.__setattr__(self, "change_ids", tuple(self.change_ids))
        object.__setattr__(self, "parent_transcript_ids", tuple(self.parent_transcript_ids))
        if not isinstance(self.accession, str) or not self.accession or re.search(r"\s", self.accession):
            raise ValueError("A derived protein needs a whitespace-free accession")
        if not isinstance(self.sequence, str) or re.fullmatch(r"[A-Z*]+", self.sequence) is None:
            raise ValueError("A derived protein needs an uppercase amino-acid sequence")
        if not self.parent_accessions or any(not isinstance(p, str) or not p or re.search(r"[\s,]", p)
                                               for p in self.parent_accessions):
            raise ValueError("A derived protein needs parent protein accessions")
        if self.parent_transcript_ids and len(self.parent_transcript_ids) != len(self.parent_accessions):
            raise ValueError("Specify one transcript per parent protein")
        if any(t is not None and (not isinstance(t, str) or not t or re.search(r"[\s,]", t))
               for t in self.parent_transcript_ids):
            raise ValueError("Parent transcript IDs must be nonempty IDs or None")
        if not isinstance(self.sample_id, str) or not self.sample_id or re.search(r"[\s,]", self.sample_id):
            raise ValueError("A derived protein needs a sample ID")
        if not self.change_ids or any(not isinstance(c, str) or not c or re.search(r"[\s,]", c)
                                      for c in self.change_ids):
            raise ValueError("A derived protein needs change IDs")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("A derived protein needs a source")
        if self.source_sha256 is not None and re.fullmatch(r"[0-9a-f]{64}", self.source_sha256) is None:
            raise ValueError("source_sha256 must be a lowercase SHA-256 hex digest")


@dataclass(frozen=True)
class ProteinDerivation:
    """Parent proteins and sample events behind one supplied protein sequence."""

    parent_accessions: tuple
    parent_sources: tuple
    sample_id: str
    change_ids: tuple
    parent_gene_ids: tuple
    parent_transcript_ids: tuple
    source_sha256: Optional[str] = None

    def __post_init__(self):
        if any(not isinstance(value, tuple) for value in
               (self.parent_accessions, self.parent_sources, self.change_ids,
                self.parent_gene_ids, self.parent_transcript_ids)):
            raise TypeError("Derivation parents, sources, changes and annotations must be tuples")
        if not self.parent_accessions or len(self.parent_accessions) != len(self.parent_sources):
            raise ValueError("A derivation needs one parent source per parent accession")
        if len(set(self.parent_accessions)) != len(self.parent_accessions):
            raise ValueError("A derivation cannot repeat a parent accession")
        if (len(self.parent_gene_ids) != len(self.parent_accessions)
                or len(self.parent_transcript_ids) != len(self.parent_accessions)):
            raise ValueError("A derivation needs gene and transcript annotations for each parent")
        if any(value is not None and (not isinstance(value, str) or not value) for value in
               self.parent_gene_ids + self.parent_transcript_ids):
            raise ValueError("Parent gene and transcript IDs must be strings or absent")
        if (any(not isinstance(value, str) or not value or re.search(r"[\s,]", value) for value in
                self.parent_accessions + self.parent_sources + self.change_ids)
                or not isinstance(self.sample_id, str) or not self.sample_id
                or re.search(r"[\s,]", self.sample_id)
                or not self.change_ids):
            raise ValueError("A derivation needs a sample and change IDs")
        if self.source_sha256 is not None and re.fullmatch(r"[0-9a-f]{64}", self.source_sha256) is None:
            raise ValueError("source_sha256 must be a lowercase SHA-256 hex digest")


@dataclass(frozen=True)
class ProducerInfo:
    """Source identity and declared context for an externally supplied protein.

    The producer supplies the amino acids separately. Gene/transcript/ORF IDs
    here are source claims, not fabricated genomic ``Occurrence`` mappings.
    ``reference=None`` means the producer has not qualified an exact reference.
    """

    origin: str
    method: str
    tool: str
    native_id: str
    tool_version: Optional[str] = None
    provider: Optional[str] = None
    exporter: Optional[str] = field(default=None, repr=False, compare=False)
    sample_id: Optional[str] = None
    run_id: Optional[str] = None
    gene_ids: tuple = ()
    gene_names: tuple = ()
    transcript_ids: tuple = ()
    orf_ids: tuple = ()
    reference: Optional[ReferenceContext] = None
    start_complete: Optional[bool] = None
    stop_complete: Optional[bool] = None
    source_sha256: Optional[str] = None
    host_taxon_ids: tuple = ()
    host_source: Optional[str] = None
    orf_kind: Optional[str] = None
    exporters: tuple[str, ...] = ()

    def __post_init__(self):
        for name in ("origin", "method", "tool", "native_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Producer {name} must be a nonempty string")
        for name in ("tool_version", "provider", "exporter", "sample_id", "run_id",
                     "host_source", "orf_kind"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"Producer {name} must be a nonempty string or None")
        for name in ("gene_ids", "gene_names", "transcript_ids", "orf_ids"):
            value = getattr(self, name)
            if isinstance(value, str):
                raise TypeError(f"Producer {name} must be a collection of IDs")
            value = tuple(value)
            if any(not isinstance(item, str) or not item.strip() for item in value):
                raise ValueError(f"Producer {name} must contain nonempty strings")
            object.__setattr__(self, name, value)
        if isinstance(self.exporters, str):
            raise TypeError("Producer exporters must be a collection of names")
        exporters = tuple(sorted(set((*self.exporters,
                                      *((self.exporter,) if self.exporter is not None else ())))))
        if any(not isinstance(value, str) or not value.strip() for value in exporters):
            raise ValueError("Producer exporters must contain nonempty strings")
        object.__setattr__(self, "exporters", exporters)
        object.__setattr__(self, "exporter", exporters[0] if len(exporters) == 1 else None)
        if self.reference is not None and not isinstance(self.reference, ReferenceContext):
            raise TypeError("Producer reference must be a ReferenceContext or None")
        for name in ("start_complete", "stop_complete"):
            if getattr(self, name) is not None and type(getattr(self, name)) is not bool:
                raise ValueError(f"Producer {name} must be True, False, or None")
        if self.source_sha256 is not None and re.fullmatch(r"[0-9a-f]{64}", self.source_sha256) is None:
            raise ValueError("Producer source_sha256 must be a lowercase SHA-256")
        if isinstance(self.host_taxon_ids, (str, int)):
            raise TypeError("Producer host_taxon_ids must be a collection")
        object.__setattr__(self, "host_taxon_ids", tuple(sorted(set(self.host_taxon_ids))))
        if any(type(taxon) is not int or taxon < 1 for taxon in self.host_taxon_ids):
            raise ValueError("Producer host taxa must be positive NCBI taxon IDs")
        if self.origin == "viral":
            if not self.host_taxon_ids or self.host_source is None:
                raise ValueError("A viral producer entry needs sourced host associations")
        elif self.host_taxon_ids or self.host_source is not None:
            raise ValueError("Only viral producer entries may declare host associations")


def _producer_accession(info, taxon_id):
    """A reversible, source-qualified ID; producer-native IDs remain separate."""
    components = dict(taxon=taxon_id, origin=info.origin, method=info.method,
                      tool=info.tool, version=info.tool_version, provider=info.provider,
                      sample=info.sample_id, run=info.run_id,
                      native=info.native_id, reference=None if info.reference is None else info.reference.key,
                      source_sha256=info.source_sha256)
    return "producer|" + "|".join(f"{name}={quote(str(value), safe='')}"
                                    for name, value in components.items() if value is not None)


def _record_key(record):
    fields = dict(namespace=record.namespace, provider=record.provider,
                  sample=record.sample_id, run=record.run_id,
                  kind=record.kind, native=record.native_id)
    return "record|" + "|".join(f"{name}={quote(value, safe='')}"
                                  for name, value in fields.items() if value is not None)


def _freeze_json(value):
    """Keep caller-owned producer payloads immutable inside a snapshot."""
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("Producer payload object keys must be strings")
        return MappingProxyType({key: _freeze_json(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    if value is None or type(value) in (str, int, float, bool):
        if type(value) is float and not (-float("inf") < value < float("inf")):
            raise ValueError("Producer payload numbers must be finite")
        return value
    raise TypeError(f"Producer payload contains a non-JSON value: {type(value).__name__}")


def _thaw_json(value):
    if isinstance(value, Mapping):
        return {key: _thaw_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


@dataclass(frozen=True)
class ProducerRecord:
    """Versioned producer record, including translations, evidence and events.

    The core retains ``payload`` without interpreting its producer schema.
    ``status`` distinguishes events without a sequence from an absent record.
    Evidence ``support`` states what was observed, independently of links to
    translations or proteins. ``scope`` identifies the observed material.
    """

    namespace: str
    schema_version: str
    kind: str
    native_id: str
    provider: Optional[str] = None
    exporter: Optional[str] = field(default=None, repr=False, compare=False)
    sample_id: Optional[str] = None
    run_id: Optional[str] = None
    status: Optional[str] = None
    support: Optional[str] = None
    scope: Mapping = field(default_factory=dict)
    complete: Optional[bool] = None
    limits: Mapping = field(default_factory=dict)
    payload: Mapping = field(default_factory=dict)
    exporters: tuple[str, ...] = ()

    def __post_init__(self):
        for name in ("namespace", "schema_version", "kind", "native_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Producer record {name} must be a nonempty string")
        for name in ("provider", "exporter", "sample_id", "run_id", "status", "support"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"Producer record {name} must be a nonempty string or None")
        if isinstance(self.exporters, str):
            raise TypeError("Producer record exporters must be a collection of names")
        exporters = tuple(sorted(set((*self.exporters,
                                      *((self.exporter,) if self.exporter is not None else ())))))
        if any(not isinstance(value, str) or not value.strip() for value in exporters):
            raise ValueError("Producer record exporters must contain nonempty strings")
        object.__setattr__(self, "exporters", exporters)
        object.__setattr__(self, "exporter", exporters[0] if len(exporters) == 1 else None)
        if self.complete is not None and type(self.complete) is not bool:
            raise ValueError("Producer record complete must be True, False, or None")
        if self.kind == "evidence" and self.support is None:
            raise ValueError("Producer evidence needs an explicit support type")
        if self.kind == "evidence" and not self.scope:
            raise ValueError("Producer evidence needs an explicit observation scope")
        if self.kind == "event" and self.status is None:
            raise ValueError("Producer event needs an explicit outcome status")
        for name in ("scope", "limits", "payload"):
            value = getattr(self, name)
            if not isinstance(value, Mapping):
                raise TypeError(f"Producer record {name} must be a JSON object")
            object.__setattr__(self, name, _freeze_json(dict(value)))

    @property
    def key(self):
        return _record_key(self)

    def to_dict(self):
        result = {name: getattr(self, name) for name in
                  ("namespace", "schema_version", "kind", "native_id", "provider",
                   "sample_id", "run_id", "status", "support", "complete")}
        result.update(exporters=list(self.exporters), scope=_thaw_json(self.scope),
                      limits=_thaw_json(self.limits),
                      payload=_thaw_json(self.payload))
        return result


@dataclass(frozen=True)
class ProducerLink:
    """A typed relationship from a producer record to a record or protein."""

    source: str
    target: str
    target_kind: str
    relation: str
    payload: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if any(not isinstance(getattr(self, name), str) or not getattr(self, name).strip()
               for name in ("source", "target", "relation")):
            raise ValueError("Producer links need source, target and relation IDs")
        if self.target_kind not in ("protein", "record"):
            raise ValueError("Producer link target_kind must be 'protein' or 'record'")
        if not isinstance(self.payload, Mapping):
            raise TypeError("Producer link payload must be a JSON object")
        object.__setattr__(self, "payload", _freeze_json(dict(self.payload)))

    def to_dict(self):
        return dict(source=self.source, target=self.target, target_kind=self.target_kind,
                    relation=self.relation, payload=_thaw_json(self.payload))


@dataclass(frozen=True)
class ProteinSequence:
    """One protein: its sequence, identity, source and every transcript encoding it.

    ``kind`` is ``"annotated"``, ``"extra_orf"``, ``"protein_set"``, ``"viral"``,
    ``"derived"`` or ``"external"``. ``accession``
    is the versioned protein ID, the source ORF ID or the viral protein ID.
    ``occurrences`` lists every transcript of this genome encoding the protein,
    the source's own transcript first; viral proteins have none. ``gene_id``,
    ``gene_name``, ``transcript_id``, ``contig`` and ``transcript_biotype``
    describe the first occurrence. ``source_binding`` explains why an extra ORF
    could not be placed on its source's own transcript; it is None when it was.
    ``orf_type`` is an extra ORF's source category, kept verbatim, and
    ``evidence`` lists its supporting evidence types; ``evidence_detail`` holds
    the source's evidence records. ``equivalent_ids`` lists other sources'
    ``source:id`` for the same extra ORF when catalogues were merged. Bundled
    viral proteins keep their source's gene and product labels separately.
    Supplied sample-derived proteins keep parent accessions and annotations,
    change IDs and sample provenance in ``derivation``; their ``occurrences``
    remain empty because reference coordinates may no longer describe them.
    """

    accession: str
    sequence: str
    kind: str
    taxon_id: Optional[int]
    source: str
    orf_type: Optional[str] = None
    evidence: tuple = ()
    evidence_detail: Optional[Mapping] = field(default=None, compare=False)
    occurrences: tuple = ()
    source_transcript_id: Optional[str] = None
    source_binding: Optional[str] = None
    equivalent_ids: tuple = ()
    orf: Optional[TranscriptORF] = field(default=None, compare=False, repr=False)
    viral_gene: Optional[str] = None
    viral_product: Optional[str] = None
    derivation: Optional[ProteinDerivation] = None
    producer: Optional[ProducerInfo] = None

    def __post_init__(self):
        if self.taxon_id is not None and (type(self.taxon_id) is not int or self.taxon_id < 1):
            raise ValueError("Protein taxon_id must be a positive NCBI taxon ID or None")
        if (self.kind == "derived") != (self.derivation is not None):
            raise ValueError("Only derived proteins have derivation metadata")
        if (self.kind == "external") != (self.producer is not None):
            raise ValueError("Only external proteins have producer metadata")
        if self.derivation is not None and self.accession in self.derivation.parent_accessions:
            raise ValueError("A derived protein cannot name itself as a parent")
        if self.producer is not None:
            if not isinstance(self.producer, ProducerInfo):
                raise TypeError("External protein producer metadata must be ProducerInfo")
            if (self.accession != _producer_accession(self.producer, self.taxon_id)
                    or self.occurrences or self.orf is not None or self.derivation is not None
                    or self.source_transcript_id is not None):
                raise ValueError("External proteins require a source-qualified ID and no fabricated occurrence")
            if not isinstance(self.sequence, str) or re.fullmatch(r"[A-Z*]+", self.sequence) is None:
                raise ValueError("An external protein needs its supplied uppercase amino-acid sequence")
            if self.producer.reference is not None and self.producer.reference.taxon_id != self.taxon_id:
                raise ValueError("Producer reference and protein taxon disagree")

    @classmethod
    def from_producer(
        cls, sequence: str, *, producer: ProducerInfo, taxon_id: int | None = None
    ) -> "ProteinSequence":
        """Accept supplied amino acids without requiring a parent or mapping.

        An exact producer reference supplies its taxon when ``taxon_id`` is
        omitted. With neither, the source taxon remains unknown.
        """
        if not isinstance(producer, ProducerInfo):
            raise TypeError("producer must be ProducerInfo")
        if taxon_id is None and producer.reference is not None:
            taxon_id = producer.reference.taxon_id
        if taxon_id is not None and (type(taxon_id) is not int or taxon_id < 1):
            raise ValueError("External protein taxon_id must be a positive NCBI taxon ID or None")
        return cls(_producer_accession(producer, taxon_id), sequence, "external", taxon_id,
                   f"producer:{producer.tool}", producer=producer)

    @property
    def _first(self):
        return self.occurrences[0] if self.occurrences else None

    @property
    def gene_id(self):
        if self.producer is not None:
            return self.producer.gene_ids[0] if len(self.producer.gene_ids) == 1 else None
        return self._first and self._first.gene_id

    @property
    def gene_name(self):
        if self.producer is not None:
            return self.producer.gene_names[0] if len(self.producer.gene_names) == 1 else None
        return self._first.gene_name if self._first else self.viral_gene

    @property
    def transcript_id(self):
        if self.producer is not None:
            return self.producer.transcript_ids[0] if len(self.producer.transcript_ids) == 1 else None
        return self._first and self._first.transcript_id

    @property
    def contig(self):
        return self._first and self._first.contig

    @property
    def transcript_biotype(self):
        return self._first and self._first.transcript_biotype

    @property
    def protein_source(self):
        """The protein source, using names independent of annotation status."""
        if self.producer is not None:
            return self.producer.origin
        if self.kind == "derived":
            sources = set(self.derivation.parent_sources)
            return next(iter(sources)) if len(sources) == 1 else "mixed"
        return {"annotated": "reference_cds", "extra_orf": "nc_orf",
                "protein_set": "protein_set", "viral": "viral"}[self.kind]

    def fasta_token(self):
        """Readable, source-qualified first FASTA token, stable offline."""
        producer = self.producer
        if producer is not None:
            names = producer.gene_names or producer.gene_ids
            method, tool, sample, native = (producer.method, producer.tool,
                                            producer.sample_id, producer.native_id)
            orf_kind = producer.orf_kind
            producer_fields = [
                ("version", producer.tool_version),
                ("provider", producer.provider),
                ("exporters", ",".join(producer.exporters) or None),
                ("reference", None if producer.reference is None else producer.reference.key),
                ("start", _known(producer.start_complete)),
                ("stop", _known(producer.stop_complete)),
            ]
        else:
            names = tuple(gene for gene in (self.gene_name or self.gene_id,) if gene)
            if self.derivation is not None:
                names = tuple(dict.fromkeys(gene for gene in
                                            self.derivation.parent_gene_ids if gene)) or names
            method = {"annotated": "annotation", "extra_orf": "catalogue",
                      "protein_set": "import", "viral": "reference",
                      "derived": "derivation"}[self.kind]
            tool = self.source.split(":", 1)[0]
            sample = None if self.derivation is None else self.derivation.sample_id
            native = self.accession
            orf_kind = self.orf_type
            producer_fields = []
        if orf_kind is None:
            orf_kind = "CDS" if self.protein_source == "reference_cds" else "unknown"
        fields = [
            ("origin", self.protein_source), ("orf", orf_kind), ("method", method),
            ("tool", tool), ("gene", ",".join(names) if names else "unknown"),
            ("aa", len(self.sequence)),
            ("taxon", self.taxon_id if self.taxon_id is not None else "unknown"),
            ("source", self.source),
            ("sample", sample), ("run", None if producer is None else producer.run_id),
            ("id", native), ("key", hashlib.sha256(self.accession.encode()).hexdigest()[:16]),
        ]
        fields[4:4] = producer_fields
        if len(names) > 1:
            fields.insert(5, ("gene_count", len(names)))
        return "|".join(f"{name}={quote(str(value), safe='')}" for name, value in fields
                        if value is not None)

    def fasta_header(self, *, token=None):
        first = self._first
        fields = dict(kind=self.kind, protein_source=self.protein_source, taxon=self.taxon_id,
                      source=self.source, gene=self.gene_id,
                      gene_name=self.gene_name, transcript=self.transcript_id, contig=self.contig,
                      transcript_biotype=self.transcript_biotype, orf_type=self.orf_type,
                      product=self.viral_product,
                      sample=None if self.derivation is None else self.derivation.sample_id,
                      producer=None if self.producer is None else self.producer.tool,
                      tool_version=None if self.producer is None else self.producer.tool_version,
                      provider=None if self.producer is None else self.producer.provider,
                      exporters=None if self.producer is None or not self.producer.exporters else
                      ",".join(self.producer.exporters),
                      method=None if self.producer is None else self.producer.method,
                      native_id=None if self.producer is None else self.producer.native_id,
                      producer_reference=None if self.producer is None or
                      self.producer.reference is None else self.producer.reference.key,
                      start_complete=None if self.producer is None else
                      _known(self.producer.start_complete),
                      stop_complete=None if self.producer is None else
                      _known(self.producer.stop_complete),
                      producer_source_sha256=None if self.producer is None else
                      self.producer.source_sha256,
                      producer_sample=None if self.producer is None else self.producer.sample_id,
                      producer_genes=None if self.producer is None or
                      len(self.producer.gene_names) < 2 else ",".join(self.producer.gene_names),
                      derived_from=None if self.derivation is None else
                      ",".join(self.derivation.parent_accessions),
                      parent_genes=None if self.derivation is None or
                      not any(self.derivation.parent_gene_ids) else
                      ",".join(gene or "-" for gene in self.derivation.parent_gene_ids),
                      parent_transcripts=None if self.derivation is None or
                      not any(self.derivation.parent_transcript_ids) else
                      ",".join(transcript or "-" for transcript in self.derivation.parent_transcript_ids),
                      changes=None if self.derivation is None else ",".join(self.derivation.change_ids),
                      source_sha256=None if self.derivation is None else self.derivation.source_sha256,
                      evidence=",".join(self.evidence) or None,
                      source_transcript=self.source_transcript_id
                      if first is not None and not first.on_source_transcript else None,
                      equivalent=",".join(self.equivalent_ids) or None,
                      transcripts=len(self.occurrences) if len(self.occurrences) > 1 else None)
        described = " ".join(f"{key}={_fasta_description_value(value)}"
                             for key, value in fields.items() if value is not None)
        return f">{self.fasta_token() if token is None else token} {described}"

    def to_dict(self):
        return dict(accession=self.accession, sequence=self.sequence, kind=self.kind,
                    taxon_id=self.taxon_id, source=self.source, orf_type=self.orf_type,
                    evidence=list(self.evidence), evidence_detail=self.evidence_detail,
                    occurrences=[[getattr(o, name) for name in _OCCURRENCE_FIELDS] for o in self.occurrences],
                    source_transcript_id=self.source_transcript_id, source_binding=self.source_binding,
                    **({"equivalent_ids": list(self.equivalent_ids)} if self.equivalent_ids else {}),
                    **({"viral_gene": self.viral_gene} if self.viral_gene is not None else {}),
                    **({"viral_product": self.viral_product} if self.viral_product is not None else {}),
                    **({"derivation": asdict(self.derivation)} if self.derivation is not None else {}),
                    **({"producer": {name: value for name, value in asdict(self.producer).items()
                                      if name != "exporter"
                                      and (name != "orf_kind" or value is not None)}}
                       if self.producer is not None else {}))

    @classmethod
    def from_dict(cls, record):
        record = dict(record)
        record["evidence"] = tuple(record.get("evidence", ()))
        record["equivalent_ids"] = tuple(record.get("equivalent_ids", ()))
        if record.get("derivation") is not None:
            derivation = record["derivation"]
            record["derivation"] = ProteinDerivation(
                parent_accessions=tuple(derivation["parent_accessions"]),
                parent_sources=tuple(derivation["parent_sources"]),
                sample_id=derivation["sample_id"], change_ids=tuple(derivation["change_ids"]),
                parent_gene_ids=tuple(derivation["parent_gene_ids"]),
                parent_transcript_ids=tuple(derivation["parent_transcript_ids"]),
                source_sha256=derivation.get("source_sha256"))
        if record.get("producer") is not None:
            producer = dict(record["producer"])
            producer["exporters"] = tuple(producer.get("exporters", ()))
            if producer.get("reference") is not None:
                producer["reference"] = ReferenceContext(**producer["reference"])
            record["producer"] = ProducerInfo(**producer)
        blocks = _OCCURRENCE_FIELDS.index("coding_blocks")
        record["occurrences"] = tuple(
            Occurrence(*values[:blocks], tuple(map(tuple, values[blocks])), *values[blocks + 1 :])
            for values in record.get("occurrences", ()))
        return cls(**record)


def _merge_protein(first, second):
    """Combine exporter paths for one source protein, rejecting biological conflicts."""
    if first.accession != second.accession:
        raise ValueError("Cannot merge proteins with different accessions")
    if first.producer is None or second.producer is None:
        if first.to_dict() != second.to_dict():
            raise ValueError(f"Conflicting protein identity: {first.accession}")
        return first
    first_base = replace(first, producer=replace(first.producer, exporter=None, exporters=()))
    second_base = replace(second, producer=replace(second.producer, exporter=None, exporters=()))
    if first_base.to_dict() != second_base.to_dict():
        raise ValueError(f"Conflicting protein identity: {first.accession}")
    exporters = tuple(sorted(set(first.producer.exporters + second.producer.exporters)))
    return replace(first, producer=replace(first.producer, exporter=None, exporters=exporters))


def _merge_producer_record(first, second):
    """Combine exporter paths for one source record, rejecting content conflicts."""
    if first.key != second.key:
        raise ValueError("Cannot merge producer records with different keys")
    if (replace(first, exporter=None, exporters=()).to_dict()
            != replace(second, exporter=None, exporters=()).to_dict()):
        raise ValueError(f"Conflicting producer record identity: {first.key}")
    exporters = tuple(sorted(set(first.exporters + second.exporters)))
    return replace(first, exporter=None, exporters=exporters)


@dataclass(frozen=True)
class ExtraORFs:
    """Catalogue ORFs on the transcripts of one genome and reference.

    Each ORF appears once for every transcript carrying its exact spliced
    coding path and protein; ``orf.on_source_transcript`` marks the transcript
    its source names. ``source_failures`` lists ORFs whose source transcript
    did not carry them, and ``excluded`` those found on no transcript at all;
    every one is a record of ``catalog``. ``unlifted`` lists ORFs the rule
    selects that a lift to this assembly left out of ``catalog``.
    ``reference`` is kept even when nothing binds. Lookups never return an ORF
    bound to another reference or version.
    """

    orfs: tuple
    excluded: tuple
    catalog: object = field(repr=False)
    reference: ReferenceContext
    evidence: Optional[EvidenceRule] = None
    source_failures: tuple = ()
    genome: object = field(default=None, repr=False, compare=False)
    unlifted: tuple = ()
    min_protein_length: Optional[int] = None
    length_filtered: int = 0

    def __post_init__(self):
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("ExtraORFs requires its resolved ReferenceContext")
        if self.min_protein_length is not None:
            _positive_integer(self.min_protein_length, "min_protein_length")
        if type(self.length_filtered) is not int or self.length_filtered < 0:
            raise ValueError("length_filtered must be a nonnegative integer")
        if self.min_protein_length is None and self.length_filtered:
            raise ValueError("length_filtered requires min_protein_length")
        if any(o.context != self.reference and o.context.key != self.reference.key for o in self.orfs):
            raise ValueError("Every ORF in a collection must be bound to its reference")
        if any(o.kind != "extra_orf" for o in self.orfs):
            raise ValueError("ExtraORFs holds catalogue ORFs; annotated CDSs come from the transcript")
        # Per-contig spans, including initiation flanks, sorted for interval lookup.
        spans = {}
        for orf in self.orfs:
            positions = [p for block in orf.coding_blocks for p in block]
            context = orf.initiation_context
            if context is not None:
                positions.extend(context.genomic_positions)
            spans.setdefault(_contig(orf.transcript.contig), []).append((min(positions), max(positions), orf))
        index = {}
        for contig, items in spans.items():
            items.sort(key=lambda item: (item[0], item[1]))
            index[contig] = ([a for a, _, _ in items], max(b - a for a, b, _ in items), items)
        object.__setattr__(self, "_index", index)

    def __iter__(self):
        return iter(self.orfs)

    def __len__(self):
        return len(self.orfs)

    def for_transcript(self, transcript, *, include_annotated=False):
        """ORFs on one transcript of this reference, in RNA order.

        A string matches a transcript ID exactly. An unversioned Ensembl ID
        means this annotation's version; a versioned one must match it. A
        transcript object must come from this reference and version.
        Incompatible requests raise rather than return another reference's ORF.
        ``include_annotated=True`` adds every annotated CDS on this transcript:
        its own, and any other transcript's CDS whose spliced coding path it
        carries, such as a bicistronic RNA's upstream protein.
        """
        if isinstance(transcript, str):
            found = tuple(o for o in self.orfs
                          if _same_id(transcript, o.transcript.id, _versioned(o.transcript)))
        else:
            self._check_genome(transcript, "Transcript")
            found = tuple(o for o in self.orfs if o.transcript.id == transcript.id)
            if any(_versioned(o.transcript) != _versioned(transcript) for o in found):
                raise ValueError(f"Transcript {_versioned(transcript)} differs from this reference's "
                                 f"{_versioned(found[0].transcript)}")
        if include_annotated:
            host = transcript if not isinstance(transcript, str) else (
                found[0].transcript if found else self._transcript_by_id(transcript))
            found += annotated_orfs(host, context=self.reference)
        return tuple(sorted(found, key=lambda o: (o.start, o.end, o.orf_id)))

    def overlapping(self, contig, start, end, *, genome_build, include_initiation_context=False):
        """ORFs whose coding path overlaps a 1-based inclusive genomic interval.

        ``genome_build`` must be this collection's assembly. Introns between
        coding blocks do not count. With ``include_initiation_context``, an
        ORF's spliced start flanks also count, for example for a Kozak variant.
        """
        if genome_build != self.reference.genome_build:
            raise ValueError(f"Expected genome build {self.reference.genome_build}, got {genome_build}")
        if type(start) is not int or type(end) is not int or start < 1 or end < start:
            raise ValueError("start and end must be positive integers with end >= start")
        starts, widest, items = self._index.get(_contig(contig), ((), 0, ()))
        found = []
        for _, _, orf in items[bisect_left(starts, start - widest) : bisect_right(starts, end)]:
            context = orf.initiation_context if include_initiation_context else None
            if any(a <= end and b >= start for a, b in orf.coding_blocks) or (
                    include_initiation_context and context is not None
                    and any(start <= p <= end for p in context.genomic_positions)):
                found.append(orf)
        return tuple(found)

    def for_orf(self, orf_id, *, source=None):
        """Every transcript placement of one ORF, the source's own transcript first.

        After merging catalogues, any source's ID finds the ORF, bare or as
        ``source:id``; a bare ID two sources use for different ORFs raises.
        Pass ``source`` to select the native ID in one source explicitly.
        """
        if self.catalog is not None:
            try:
                orf_id = self.catalog.get(orf_id, source=source).orf_id
            except KeyError:
                return ()
        elif source is not None:
            raise ValueError("Source-qualified lookup requires the source catalogue")
        return tuple(sorted((o for o in self.orfs if o.orf_id == orf_id),
                            key=lambda o: (not o.on_source_transcript, _versioned(o.transcript))))

    def for_gene(self, gene, *, source=False):
        """ORFs on transcripts of one gene in this genome.

        ``gene`` is a gene object, ID or name, resolved like transcripts. With
        ``source=True``, a string instead matches the source catalogue's own
        gene ID or name, which can differ from this annotation's.
        """
        if source:
            if not isinstance(gene, str):
                raise TypeError("Source gene labels are strings; a gene object selects this genome's genes")
            return tuple(o for o in self.orfs if gene and any(
                gene in (m.gene_id, m.gene_name) for m in (o.record.members or (o.record,))))
        if isinstance(gene, str):
            return tuple(o for o in self.orfs if gene and (
                gene == o.transcript.gene_name
                or _same_id(gene, o.transcript.gene_id, lambda o=o: _versioned(o.transcript.gene))))
        self._check_genome(gene, "Gene")
        found = tuple(o for o in self.orfs if o.transcript.gene_id == gene.id)
        if any(_versioned(o.transcript.gene) != _versioned(gene) for o in found):
            raise ValueError(f"Gene {_versioned(gene)} differs from this reference's "
                             f"{_versioned(found[0].transcript.gene)}")
        return found

    def _transcript_by_id(self, transcript_id):
        if self.genome is None:
            raise ValueError("Pass a transcript object: this collection has no genome to resolve IDs")
        transcript = self.genome.transcript_by_id(transcript_id)
        self._check_genome(transcript, "Transcript")
        return transcript

    def _check_genome(self, feature, label):
        try:
            self.reference.validate_genome(feature.genome)
        except ValueError as error:
            raise ValueError(f"{label} {feature.id} is not from this collection's reference: {error}") from None


def annotated_orfs(transcript, *, context=None):
    """Every annotated CDS on one transcript, as ``TranscriptORF`` objects.

    That is the transcript's own CDS and any other transcript's CDS whose exact
    spliced coding path and sequence this transcript carries, such as another
    gene's protein in a bicistronic RNA. Annotated proteins with an identical
    placement and sequence are reported once: the transcript's own if it is
    among them, with the others' protein IDs in ``equivalent_ids``.
    """
    found = {}
    candidates = transcript.genome.transcripts_at_locus(
        transcript.contig, transcript.start, transcript.end, strand=transcript.strand)
    for other in candidates:
        if not other.protein_id:
            continue
        try:
            orf = TranscriptORF.from_annotation(
                other, on=None if other.id == transcript.id else transcript, context=context)
        except ValueError:
            continue  # This transcript does not carry that CDS.
        found.setdefault((orf.start, orf.end, orf.protein_sequence), []).append(orf)
    placements = []
    for same in found.values():
        same.sort(key=lambda o: (not o.on_source_transcript, o.source_transcript_id))
        first = same[0]
        first._set(equivalent_ids=tuple(o.orf_id for o in same[1:]))
        placements.append(first)
    return tuple(sorted(placements, key=lambda o: (o.start, o.end, o.orf_id)))


@dataclass(frozen=True)
class ProteinSet:
    """Named protein sequences from one source, without inferred placement or use.

    An exact ``reference`` supplies its taxon when ``taxon_id`` is omitted.
    """

    name: str
    proteins: tuple[tuple[str, str], ...]
    source: str
    taxon_id: Optional[int] = None
    reference: Optional[ReferenceContext] = None
    source_sha256: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "proteins", tuple(tuple(p) for p in self.proteins))
        if not all(isinstance(value, str) and value.strip() for value in (self.name, self.source)):
            raise ValueError("A protein set requires a name and source")
        if self.taxon_id is not None and (type(self.taxon_id) is not int or self.taxon_id < 1):
            raise ValueError("Protein set taxon_id must be a positive NCBI taxon ID or None")
        if self.reference is not None and not isinstance(self.reference, ReferenceContext):
            raise TypeError("Protein set reference must be a ReferenceContext or None")
        if self.reference is not None and self.taxon_id is None:
            object.__setattr__(self, "taxon_id", self.reference.taxon_id)
        if self.reference is not None and self.reference.taxon_id != self.taxon_id:
            raise ValueError("Protein set reference and taxon disagree")
        if self.source_sha256 is not None and re.fullmatch(r"[0-9a-f]{64}", self.source_sha256) is None:
            raise ValueError("Protein set source_sha256 must be a lowercase SHA-256")
        if not self.proteins:
            raise ValueError("A protein set requires proteins")
        for accession, sequence in self.proteins:
            if (not isinstance(accession, str) or not accession or re.search(r"\s", accession)
                    or not isinstance(sequence, str)
                    or re.fullmatch(r"[A-Z*]+", sequence) is None):
                raise ValueError(f"Invalid protein record: {accession!r}")
        if len({a for a, _ in self.proteins}) != len(self.proteins):
            raise ValueError("Duplicate protein accessions")

    @classmethod
    def from_fasta(
        cls,
        path: str | Path,
        *,
        name: str,
        taxon_id: int | None = None,
        reference: ReferenceContext | None = None,
        source: str | None = None,
    ) -> "ProteinSet":
        """Read a protein FASTA; a terminal stop character is removed."""
        raw = Path(path).read_bytes()
        text = (gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw).decode("ascii")
        proteins, accession, chunks = [], None, []
        for line in text.splitlines() + [">"]:
            if line.startswith(">"):
                if accession is not None:
                    proteins.append((accession, "".join(chunks).upper().removesuffix("*")))
                accession, chunks = (line[1:].split() or [None])[0], []
            elif line.strip():
                chunks.append(line.strip())
        return cls(name=name, proteins=tuple(proteins), source=source or Path(path).name,
                   taxon_id=taxon_id, reference=reference,
                   source_sha256=hashlib.sha256(raw).hexdigest())

    def to_dict(self) -> dict[str, object]:
        """Return the set's source identity and sequences as plain data."""
        return {
            "name": self.name,
            "proteins": [list(protein) for protein in self.proteins],
            "source": self.source,
            "taxon_id": self.taxon_id,
            "reference": None if self.reference is None else self.reference.to_dict(),
            "source_sha256": self.source_sha256,
        }


@dataclass(frozen=True)
class _ViralInput:
    proteins: ProteinSet
    host_taxon_ids: frozenset[int]
    host_source: str
    products: tuple[ViralProduct, ...] = ()


def _viral_reference_input(reference: ViralReference) -> _ViralInput:
    if not isinstance(reference, ViralReference):
        raise TypeError("Expected a ViralReference")
    proteins = ProteinSet(
        name=f"{reference.virus_type}@{reference.id}",
        proteins=tuple((orf.orf_id, orf.protein_sequence) for orf in reference.orfs),
        source=reference.source_accession,
        taxon_id=reference.taxon_id,
        source_sha256=reference.source_sha256,
    )
    return _ViralInput(proteins, frozenset(reference.host_taxon_ids),
                       reference.host_source, reference.products)


def _canonical_viral_references(references):
    """Collapse coordinate placements of the same proteins, preferring a native source ID."""
    chosen = {}
    for reference in references:
        signature = tuple(sorted(_viral_reference_input(reference).proteins.proteins))
        current = chosen.get(reference.taxon_id)
        if current is not None and current[1] != signature:
            raise ValueError(f"Bundled references for viral taxon {reference.taxon_id} have different proteins; "
                             "choose an exact reference ID with viral= instead of 'all'")
        if current is None or (reference.id == reference.source_accession,
                               reference.id) > (current[0].id == current[0].source_accession,
                                                current[0].id):
            chosen[reference.taxon_id] = reference, signature
    return tuple(chosen[taxon_id][0] for taxon_id in sorted(chosen))


def _viral_inputs(selection, host_taxon_id, *, host_taxon_ids=None, host_source=None):
    """Resolve caller selections to one validated protein input per viral source."""
    if selection is None or selection is False:
        return ()
    if selection is True or selection == "all":
        references = viral_references(host=host_taxon_id)
        if not references:
            raise ValueError(f"No bundled viral references for host taxon {host_taxon_id}")
        selection = _canonical_viral_references(references)
    elif isinstance(selection, (ProteinSet, ViralReference, str)):
        selection = (selection,)
    else:
        try:
            selection = tuple(selection)
        except TypeError as error:
            raise TypeError("viral must be a reference ID, virus type, ViralReference, "
                            "ProteinSet, 'all', or a sequence of them") from error
    if "all" in selection:
        if len(selection) != 1:
            raise ValueError("Choose viral='all' or specific viruses, not both")
        return _viral_inputs("all", host_taxon_id)
    supplied_hosts = None
    if host_taxon_ids is not None or host_source is not None:
        if host_taxon_ids is None or not isinstance(host_source, str) or not host_source.strip():
            raise ValueError("External viral proteins require host_taxon_ids and host_source")
        if isinstance(host_taxon_ids, (str, int)):
            raise TypeError("host_taxon_ids must be a collection of positive NCBI taxon IDs")
        supplied_hosts = frozenset(host_taxon_ids)
        if not supplied_hosts or any(type(taxon) is not int or taxon < 1 for taxon in supplied_hosts):
            raise ValueError("host_taxon_ids must contain positive NCBI taxon IDs")
    resolved = []
    for item in selection:
        if isinstance(item, ProteinSet):
            if supplied_hosts is None:
                raise ValueError("A ProteinSet passed to with_viral requires sourced host taxa")
            resolved.append(_ViralInput(item, supplied_hosts, host_source))
        elif isinstance(item, ViralReference):
            resolved.append(_viral_reference_input(item))
        elif isinstance(item, str):
            try:
                reference = viral_reference(item)
            except ValueError:
                matches = viral_references(host=host_taxon_id, virus_type=item)
                if not matches:
                    raise ValueError(f"No bundled viral reference ID or virus type matches {item!r} "
                                     f"for host taxon {host_taxon_id}") from None
                canonical = _canonical_viral_references(matches)
                if len(canonical) != 1:
                    raise ValueError(f"Virus type {item!r} matches multiple viral taxa; "
                                     "choose an exact reference ID")
                reference = canonical[0]
            resolved.append(_viral_reference_input(reference))
        else:
            raise TypeError("viral entries must be reference IDs, virus types, ViralReference "
                            "or ProteinSet objects")
    if supplied_hosts is not None and any(not isinstance(item, ProteinSet) for item in selection):
        raise ValueError("host_taxon_ids and host_source apply only to ProteinSet inputs")
    return tuple(resolved)


@dataclass(frozen=True)
class ReferenceTranscriptome:
    """One installed PyEnsembl annotation's transcript models, without protein selection.

    The genome supplies transcript structures and sequences. A proteome is a
    projection of selected protein sources onto these models; viral proteins
    can be added as a separate source for the reference's host species.
    """

    genome: object = field(repr=False, compare=False)
    reference: ReferenceContext

    def __post_init__(self):
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("A reference transcriptome requires a ReferenceContext")
        _require_installed(self.genome)
        self.reference.validate_genome(self.genome)

    @classmethod
    def from_genome(cls, genome, *, context=None):
        return cls(genome, ReferenceContext.from_genome(genome, context=context))

    def proteome(self, *, reference_cds=True, isoforms="all", nc_orf_sources=(), viral=(),
                 derived=(), min_protein_length=None, cache=False, cache_root=None):
        """Reference CDS proteins with optional ncORF and viral sources.

        ``isoforms="all"`` keeps every reference protein; ``"longest"`` keeps
        one per gene, chosen by protein length. ``nc_orf_sources`` accepts a
        string or list: ``"all"`` for every compatible ncORF dataset,
        ``"ms"`` or ``"ribosome"`` for evidence-supported records across those
        datasets, or dataset names such as ``"bundled"``. Lists take the union.
        Evidence rules and loaded ORFCatalog objects are also accepted.
        The default empty selection adds no ncORFs. ``viral`` selects
        host-associated viruses by name. Set ``reference_cds=False`` to build
        from only the requested additions. See :meth:`Proteome.build`.
        """
        return Proteome.build(
            self, reference_cds=reference_cds, isoforms=isoforms, nc_orf_sources=nc_orf_sources,
            viral=viral, derived=derived, min_protein_length=min_protein_length,
            cache=cache, cache_root=cache_root)


class Proteome:
    """Selected protein sources for one reference, built in explicit layers.

    Build from :class:`ReferenceTranscriptome`, or start with :meth:`from_genome`
    and add layers with :meth:`with_nc_orfs`, :meth:`with_viral` and
    :meth:`with_derived_proteins`.
    Each method returns a new proteome. ``layers`` records
    every addition with its role (``"host"``, ``"viral"``, ``"component"``,
    ``"external"``, ``"derived"`` or ``"filter"``),
    source and selection; ``excluded`` lists selected ORFs found on no
    transcript. :meth:`save` writes a checksummed snapshot that :meth:`load`
    reads back offline, and ``identity`` changes whenever the content does.
    """

    def __init__(self, reference, proteins: Iterable[ProteinSequence], layers=(), excluded=(), genome=None,
                 reference_policy="strict", records=(), links=()):
        if reference is not None and not isinstance(reference, ReferenceContext):
            raise TypeError("A proteome reference must be a ReferenceContext or None")
        if reference is None and reference_policy != "mixed":
            raise ValueError("An unscoped proteome requires reference_policy='mixed'")
        if reference is not None and reference_policy not in ("strict", "allow_unknown"):
            raise ValueError("A scoped proteome requires 'strict' or 'allow_unknown' reference policy")
        self.reference = reference
        self.reference_policy = reference_policy
        self.proteins: tuple[ProteinSequence, ...] = tuple(proteins)
        self.layers = tuple(MappingProxyType(_json(dict(layer))) for layer in layers)
        self.excluded = tuple(excluded)
        self.genome = genome
        self.records = tuple(records)
        self.links = tuple(links)
        accessions = Counter(e.accession for e in self.proteins)
        duplicated = sorted(a for a, n in accessions.items() if n > 1)
        if duplicated:
            raise ValueError("Duplicate proteome accessions: " + ", ".join(duplicated[:5]))
        viral_layers = [layer for layer in self.layers if layer["role"] == "viral"]
        hosts = {layer["source"]: layer.get("host_taxon_ids", ()) for layer in viral_layers}
        if len(hosts) != len(viral_layers):
            raise ValueError("Duplicate viral protein source")
        component_layers = [layer for layer in self.layers if layer["role"] == "component"]
        components = {layer["source"]: layer for layer in component_layers}
        if len(components) != len(component_layers):
            raise ValueError("Duplicate protein set source")
        preceding = {}
        for entry in self.proteins:
            if entry.kind == "viral":
                if reference is None or reference.taxon_id not in hosts.get(entry.source, ()):
                    raise ValueError(f"Viral protein {entry.accession} has no host layer for this species")
            elif entry.kind == "protein_set":
                component = components.get(entry.source)
                if component is None:
                    raise ValueError(f"Protein {entry.accession} has no component layer")
                if (reference is not None
                        and component.get("source_reference") != reference.key
                        and reference.taxon_id not in component.get("target_taxon_ids", ())):
                    raise ValueError(
                        f"Protein set {component['name']} is not admitted for this species"
                    )
            elif entry.kind == "external":
                source = entry.producer
                if reference is not None:
                    if source.origin == "viral":
                        if reference.taxon_id not in source.host_taxon_ids:
                            raise ValueError(f"Viral producer protein {entry.accession} has no recorded host")
                    elif entry.taxon_id is None:
                        if reference_policy != "allow_unknown":
                            raise ValueError("External protein has unknown source taxon; "
                                             "use compatibility='allow_unknown' explicitly")
                    elif entry.taxon_id != reference.taxon_id:
                        raise ValueError(f"External protein {entry.accession} is from another species")
                    if source.reference is None and reference_policy != "allow_unknown":
                        raise ValueError("External protein has no qualified reference; "
                                         "use compatibility='allow_unknown' explicitly")
                    if (source.reference is not None and source.origin != "viral"
                            and source.reference.key != reference.key):
                        raise ValueError(f"External protein {entry.accession} declares another reference")
            elif reference is None or entry.taxon_id != reference.taxon_id:
                raise ValueError(f"Host protein {entry.accession} is from another species")
            if entry.kind == "derived":
                for accession, source, gene_id, transcript_id in zip(
                        entry.derivation.parent_accessions, entry.derivation.parent_sources,
                        entry.derivation.parent_gene_ids, entry.derivation.parent_transcript_ids):
                    parent = preceding.get(accession)
                    if parent is not None:
                        if parent.protein_source == "viral":
                            raise ValueError("Derived proteins from viral parents are unsupported")
                        occurrence = next((o for o in parent.occurrences
                                           if o.transcript_id == transcript_id), None)
                        if (parent.protein_source != source or
                                (transcript_id is None and (gene_id is not None or parent.occurrences)) or
                                (transcript_id is not None and
                                 (occurrence is None or occurrence.gene_id != gene_id))):
                            raise ValueError(f"Derived protein {entry.accession} has incorrect parent annotations")
            preceding[entry.accession] = entry
        for layer in self.layers:
            if layer.get("reference") not in (None, None if reference is None else reference.key):
                raise ValueError("A layer declares a different reference")
        if any(not isinstance(record, ProducerRecord) for record in self.records):
            raise TypeError("Proteome records must be ProducerRecord objects")
        record_keys = {record.key for record in self.records}
        if len(record_keys) != len(self.records):
            raise ValueError("Duplicate producer record identity")
        if any(not isinstance(link, ProducerLink) for link in self.links):
            raise TypeError("Proteome links must be ProducerLink objects")
        if len({json.dumps(link.to_dict(), sort_keys=True) for link in self.links}) != len(self.links):
            raise ValueError("Duplicate producer links")
        for link in self.links:
            if link.source not in record_keys:
                raise ValueError(f"Producer link has unknown source: {link.source}")
            targets = accessions if link.target_kind == "protein" else record_keys
            if link.target not in targets:
                raise ValueError(f"Producer link has unknown target: {link.target}")

    @classmethod
    def from_proteins(
        cls,
        proteins: ProteinSequence | Iterable[ProteinSequence],
        *,
        reference: ReferenceContext | None = None,
        compatibility: str = "strict",
        records: Iterable[ProducerRecord] = (),
        links: Iterable[ProducerLink] = (),
    ) -> "Proteome":
        """Make a proteome from external proteins, with an explicit reference policy.

        Use ``compatibility='mixed'`` for a standalone collection whose
        reference is unknown or mixed. With a host reference,
        ``compatibility='allow_unknown'`` explicitly admits producer proteins
        lacking an exact source reference; known conflicting references fail.
        """
        result = cls(reference, (), reference_policy=compatibility).with_proteins(proteins)
        return result.with_records(records, links=links)

    def with_proteins(
        self,
        proteins: ProteinSequence | Iterable[ProteinSequence],
        *,
        compatibility: str | None = None,
    ) -> "Proteome":
        """Add producer proteins without inventing parents or genomic occurrences."""
        proteins = (proteins,) if isinstance(proteins, ProteinSequence) else tuple(proteins)
        if any(not isinstance(protein, ProteinSequence) or protein.kind != "external"
               for protein in proteins):
            raise TypeError("with_proteins requires external ProteinSequence records")
        policy = self.reference_policy if compatibility is None else compatibility
        combined = list(self.proteins)
        positions = {protein.accession: i for i, protein in enumerate(combined)}
        added = []
        changed = policy != self.reference_policy
        for protein in proteins:
            position = positions.get(protein.accession)
            if position is not None:
                merged = _merge_protein(combined[position], protein)
                if merged != combined[position]:
                    combined[position] = merged
                    changed = True
                continue
            positions[protein.accession] = len(combined)
            combined.append(protein)
            added.append(protein)
        if not added and not changed:
            return self
        layers = self.layers
        if added:
            layer = dict(kind="external", role="external", source="producer",
                         producers=sorted({entry.producer.tool for entry in added}), added=len(added),
                         reference=None if self.reference is None else self.reference.key)
            layers += (layer,)
        return Proteome(self.reference, combined, layers,
                        self.excluded, self.genome, reference_policy=policy,
                        records=self.records, links=self.links)

    def with_protein_set(
        self,
        proteins: ProteinSet,
        *,
        role: str,
        target_taxon_ids: Iterable[int] | None = None,
        target_source: str | None = None,
    ) -> "Proteome":
        """Add an explicitly selected protein set without inventing genomic placement.

        A set on this proteome's exact reference needs no extra assertion. A
        set with another or unknown reference requires ``target_taxon_ids`` and
        ``target_source`` to say why it belongs in this target proteome. The
        source taxon and reference remain unchanged.
        """
        if not isinstance(proteins, ProteinSet):
            raise TypeError("with_protein_set requires a ProteinSet")
        if not isinstance(role, str) or not role.strip():
            raise ValueError("A protein set role must be a nonempty string")
        if isinstance(target_taxon_ids, (str, int)):
            raise TypeError("target_taxon_ids must be a collection of positive NCBI taxon IDs")
        targets = frozenset(() if target_taxon_ids is None else target_taxon_ids)
        if any(type(taxon) is not int or taxon < 1 for taxon in targets):
            raise ValueError("target_taxon_ids must contain positive NCBI taxon IDs")
        if target_source is not None and (not isinstance(target_source, str) or not target_source.strip()):
            raise ValueError("target_source must be a nonempty string or None")
        if bool(targets) != (target_source is not None):
            raise ValueError("target_taxon_ids and target_source must be supplied together")

        if self.reference is not None:
            same_reference = (proteins.reference is not None
                              and proteins.reference.key == self.reference.key)
            if (targets or not same_reference) and self.reference.taxon_id not in targets:
                raise ValueError(
                    f"Protein set {proteins.name!r} needs an explicit target assertion for "
                    f"{self.reference.species} (taxon {self.reference.taxon_id})"
                )

        source = proteins.name
        added = tuple(ProteinSequence(accession, sequence, "protein_set", proteins.taxon_id, source)
                      for accession, sequence in proteins.proteins)
        layer = dict(
            kind="protein_set", role="component", component_role=role, name=proteins.name,
            source=source, taxon_id=proteins.taxon_id,
            source_reference=None if proteins.reference is None else proteins.reference.key,
            target_taxon_ids=sorted(targets), target_source=target_source,
            provenance=proteins.source, source_sha256=proteins.source_sha256,
            added=len(added),
        )
        existing_layers = [item for item in self.layers
                           if item.get("role") == "component" and item.get("source") == source]
        if existing_layers:
            if len(existing_layers) != 1 or dict(existing_layers[0]) != layer:
                raise ValueError(f"Conflicting protein set: {proteins.name}")
            existing = self.by_accession()
            if any(entry.accession not in existing
                   or existing[entry.accession].to_dict() != entry.to_dict() for entry in added):
                raise ValueError(f"Conflicting protein set: {proteins.name}")
            return self
        collisions = [entry.accession for entry in added
                      if entry.accession in self.by_accession()]
        if collisions:
            raise ValueError("Protein set accessions already exist: " + ", ".join(collisions[:5]))
        return Proteome(
            self.reference, self.proteins + added, self.layers + (layer,),
            self.excluded, self.genome, reference_policy=self.reference_policy,
            records=self.records, links=self.links,
        )

    def with_records(self, records=(), *, links=()):
        """Add portable producer records and relationships to existing proteins."""
        records = (records,) if isinstance(records, ProducerRecord) else tuple(records)
        links = (links,) if isinstance(links, ProducerLink) else tuple(links)
        if any(not isinstance(record, ProducerRecord) for record in records):
            raise TypeError("with_records requires ProducerRecord objects")
        if any(not isinstance(link, ProducerLink) for link in links):
            raise TypeError("with_records requires ProducerLink objects")
        combined = list(self.records)
        positions = {record.key: i for i, record in enumerate(combined)}
        added = []
        changed = False
        for record in records:
            position = positions.get(record.key)
            if position is not None:
                merged = _merge_producer_record(combined[position], record)
                if merged != combined[position]:
                    combined[position] = merged
                    changed = True
                continue
            positions[record.key] = len(combined)
            combined.append(record)
            added.append(record)
        seen_links = {json.dumps(link.to_dict(), sort_keys=True) for link in self.links}
        added_links = []
        for link in links:
            key = json.dumps(link.to_dict(), sort_keys=True)
            if key not in seen_links:
                added_links.append(link)
                seen_links.add(key)
        if not added and not added_links and not changed:
            return self
        return Proteome(self.reference, self.proteins, self.layers, self.excluded, self.genome,
                        reference_policy=self.reference_policy,
                        records=combined, links=self.links + tuple(added_links))

    def merge(self, *others, compatibility=None):
        """Compose compatible proteomes; equal source records are idempotent."""
        result = self
        for other in others:
            if not isinstance(other, Proteome):
                raise TypeError("merge requires Proteome collections")
            if result.reference is not None and other.reference is not None:
                if result.reference.key != other.reference.key:
                    raise ValueError("Cannot merge different reference contexts")
            elif (result.reference is None) != (other.reference is None) and compatibility is None:
                raise ValueError("Merging a mixed collection into a reference needs explicit compatibility")
            reference = result.reference or other.reference
            policy = (compatibility if compatibility is not None else
                      "mixed" if reference is None else
                      "allow_unknown" if "allow_unknown" in (result.reference_policy,
                                                            other.reference_policy) else "strict")
            entries = list(result.proteins)
            positions = {entry.accession: i for i, entry in enumerate(entries)}
            for entry in other.proteins:
                position = positions.get(entry.accession)
                if position is not None:
                    entries[position] = _merge_protein(entries[position], entry)
                    continue
                positions[entry.accession] = len(entries)
                entries.append(entry)
            layers = list(result.layers)
            layer_keys = {_sha256(dict(layer)) for layer in layers}
            for layer in other.layers:
                key = _sha256(dict(layer))
                if key not in layer_keys:
                    layers.append(layer)
                    layer_keys.add(key)
            excluded = list(result.excluded)
            for failure in other.excluded:
                if failure not in excluded:
                    excluded.append(failure)
            ordered_records = list(result.records)
            record_positions = {record.key: i for i, record in enumerate(ordered_records)}
            for record in other.records:
                position = record_positions.get(record.key)
                if position is not None:
                    ordered_records[position] = _merge_producer_record(
                        ordered_records[position], record
                    )
                    continue
                record_positions[record.key] = len(ordered_records)
                ordered_records.append(record)
            links = list(result.links)
            link_keys = {json.dumps(link.to_dict(), sort_keys=True) for link in links}
            for link in other.links:
                key = json.dumps(link.to_dict(), sort_keys=True)
                if key not in link_keys:
                    links.append(link)
                    link_keys.add(key)
            result = Proteome(reference, entries, layers, excluded,
                              result.genome or other.genome, reference_policy=policy,
                              records=ordered_records, links=links)
        return result

    @classmethod
    def from_genome(cls, genome, *, isoforms="all", context=None):
        """Reference CDS proteins with sequences in this PyEnsembl genome.

        ``isoforms="all"`` matches the annotation's protein FASTA, including all transcript
        biotypes with a translation and proteins on patch or haplotype contigs.
        ``isoforms="longest"`` selects one protein per gene by amino-acid length,
        excluding a terminal stop. Ties use accession, then transcript ID in
        lexical order. This is a length-based representative, not an Ensembl
        canonical or MANE designation. The selection is recorded in the layer.
        Use :meth:`filter` to restrict them, for example by ``contig``.
        """
        if isoforms not in ("all", "longest"):
            raise ValueError("isoforms must be 'all' or 'longest'")
        _require_installed(genome)
        reference = ReferenceContext.from_genome(genome, context=context)
        source = f"{reference.annotation_name}:{reference.annotation_version}"
        entries = []
        for structure in _coding_structures(genome):
            accession = structure.protein_id + (
                "" if structure.protein_version is None else f".{structure.protein_version}")
            sequence = genome.protein_sequence(accession)
            if not sequence:
                continue
            offset = _layout_offset(structure.layout, structure.strand, structure.cds_blocks)
            length = sum(b - a + 1 for a, b in structure.cds_blocks)
            occurrence = Occurrence(
                structure.versioned_id, structure.gene_id, structure.gene_name, structure.biotype,
                structure.contig, structure.strand, offset, None if offset is None else offset + length,
                structure.cds_blocks if offset is not None else (), True)
            entries.append(ProteinSequence(
                accession=accession, sequence=sequence, kind="annotated", taxon_id=reference.taxon_id,
                source=source, occurrences=(occurrence,), source_transcript_id=structure.versioned_id))
        available = len(entries)
        if isoforms == "longest":
            by_gene = {}
            for entry in sorted(entries, key=lambda e: (
                    -len(e.sequence.rstrip("*")), e.accession, e.source_transcript_id)):
                if not entry.gene_id:
                    raise ValueError(f"isoforms='longest' requires a gene ID for {entry.accession}")
                by_gene.setdefault(entry.gene_id, entry)
            entries = list(by_gene.values())
        entries.sort(key=lambda e: e.source_transcript_id)
        layer = dict(kind="annotated", role="host", source=source, reference=reference.key,
                     isoforms=isoforms, available=available, added=len(entries))
        return cls(reference, entries, (layer,), genome=genome)

    @classmethod
    def from_refseq(cls, protein_fasta, feature_table, *, reference, expected_sha256):
        """Read pinned RefSeq proteins and NCBI Gene/transcript occurrences offline.

        ``expected_sha256`` requires ``protein`` and ``feature_table`` byte hashes.
        The NCBI genome-download feature table supplies accession joins, not
        spliced CDS coordinates. All FASTA proteins must have source mappings.
        """
        from .refseq import refseq_proteome

        return refseq_proteome(cls, protein_fasta, feature_table,
                              reference=reference, expected_sha256=expected_sha256)

    @classmethod
    def from_ensembl(cls, protein_fasta, gtf, rna_fasta, ncrna_fasta, *, reference, expected_sha256):
        """Join pinned Ensembl GTF, cDNA/ncRNA and peptide files without indexing.

        Keep the complete peptide inventory and every source transcript locus.
        ``expected_sha256`` requires ``protein``, ``gtf``, ``rna`` and ``ncrna``.
        RNA/header identities are checked; RNA abundance and empirical translation
        are not inferred. Original transcript spans remain in the layer manifest;
        spliced CDS coordinates are not reconstructed by this accession join.
        """
        from .ensembl import ensembl_proteome

        return ensembl_proteome(cls, protein_fasta, gtf, rna_fasta, ncrna_fasta,
                               reference=reference, expected_sha256=expected_sha256)

    @classmethod
    def build(cls, genome, *, reference_cds=True, isoforms="all", nc_orf_sources=(),
              viral=(), derived=(), min_protein_length=None, context=None, cache=False, cache_root=None):
        """Reference CDS proteins with optional sources, optionally cached.

        ``genome`` may be a PyEnsembl genome or ReferenceTranscriptome.
        ``reference_cds=True`` includes reference CDS proteins; ``isoforms``
        is ``"all"`` or ``"longest"`` (one protein per gene by length).

        ``nc_orf_sources`` is the sole ncORF selector: an empty list/tuple or
        None adds none; ``"all"`` includes every compatible registered dataset
        including optional ones; ``"ms"``, ``"ribosome"`` and other evidence
        names select supported records across all those datasets. A dataset
        name (``"bundled"`` for the bundled GENCODE data), an EvidenceRule or
        a loaded ORFCatalog is also a source. Lists combine selections by
        union, so ["bundled", "ms"] adds every bundled record and MS-supported
        records from other datasets. Set ``reference_cds=False`` for ncORFs
        alone. Required datasets and liftover chains must already be installed;
        missing resources raise with installation instructions.
        ``min_protein_length`` optionally filters only the ncORF additions; the
        source catalogues and reference CDS layer remain complete.

        ``viral`` accepts host-associated virus names, reference IDs or objects;
        ``"all"`` selects every distinct bundled viral protein set for this host.
        ``derived`` accepts supplied DerivedProtein records linked to parents.
        ``cache=True`` saves a verified snapshot keyed by the reference, source
        content, selection, isoforms and package versions. ``cache_root`` also
        controls where datasets and chains are found. Nothing is downloaded.
        """
        if isinstance(genome, ReferenceTranscriptome):
            if context is not None and context != genome.reference:
                raise ValueError("context differs from the reference transcriptome")
            context, genome = genome.reference, genome.genome
        if type(reference_cds) is not bool:
            raise TypeError("reference_cds must be a boolean")
        if isoforms not in ("all", "longest"):
            raise ValueError("isoforms must be 'all' or 'longest'")
        _require_installed(genome)
        reference = ReferenceContext.from_genome(genome, context=context)
        viral = _viral_inputs(viral, reference.taxon_id)
        derived = (() if derived is None else (derived,) if isinstance(derived, DerivedProtein)
                   else tuple(derived))
        if any(not isinstance(item, DerivedProtein) for item in derived):
            raise TypeError("derived must contain DerivedProtein records")
        catalog, rule = _nc_orf_inputs(reference, nc_orf_sources, cache_root)
        include_nc_orfs = catalog is not None
        if min_protein_length is not None:
            _positive_integer(min_protein_length, "min_protein_length")
            if not include_nc_orfs:
                raise ValueError("min_protein_length requires nc_orf_sources")
        built = []

        def make():
            proteome = (cls.from_genome(genome, isoforms=isoforms, context=context) if reference_cds else
                        cls(reference, (), genome=genome))
            if include_nc_orfs:
                proteome = proteome.with_extra_orfs(
                    rule, catalog=catalog, min_protein_length=min_protein_length)
            for virus in viral:
                proteome = proteome._with_viral_proteins(virus)
            if derived:
                proteome = proteome.with_derived_proteins(derived)
            built.append(proteome)
            return proteome

        if not cache:
            return make()
        from .datasets import _digest
        from .version import __version__

        extra_selection = (dict(evidence=None if rule is None else rule.to_dict(),
                                catalog=_catalog_identity(catalog)) if include_nc_orfs else None)
        if extra_selection is not None and min_protein_length is not None:
            extra_selection["min_protein_length"] = min_protein_length
        key = dict(schema=SNAPSHOT_SCHEMA, orfdata=__version__, pyensembl=_pyensembl_version(),
                   reference=reference.to_dict(), annotation=_annotation_files(genome, include_proteins=reference_cds),
                   annotated=reference_cds, isoforms=isoforms if reference_cds else None,
                   extra_orfs=extra_selection,
                   viral=[_digest(_json(dict(asdict(virus), host_taxon_ids=sorted(virus.host_taxon_ids))))
                          for virus in viral],
                   derived=[asdict(item) for item in derived])
        path = _cached_snapshot(key, make, cache_root)
        if built:
            return built[0]
        loaded = cls.load(path, _verified=True)
        proteome = cls(loaded.reference, loaded.proteins, loaded.layers, loaded.excluded, genome,
                       reference_policy=loaded.reference_policy,
                       records=loaded.records, links=loaded.links)
        proteome.__dict__["_identity"] = loaded.identity
        return proteome

    def with_nc_orfs(self, nc_orf_sources="all", *, min_protein_length=None, cache_root=None):
        """Add ncORFs selected by names, evidence presets/rules or loaded datasets.

        Uses the same selector as :meth:`build`; a list takes the union and
        an empty selection returns this proteome unchanged.
        """
        if nc_orf_sources is None or isinstance(nc_orf_sources, (list, tuple)) and not nc_orf_sources:
            if min_protein_length is not None:
                raise ValueError("min_protein_length requires nc_orf_sources")
            return self
        if self.genome is None:
            raise ValueError("ncORFs bind to a genome; build this proteome with from_genome")
        catalog, rule = _nc_orf_inputs(self.reference, nc_orf_sources, cache_root)
        if catalog is None:
            return self
        return self.with_extra_orfs(
            rule, catalog=catalog, min_protein_length=min_protein_length, cache_root=cache_root)

    def with_extra_orfs(self, evidence="ms", *, catalog=None, require_resolved=False,
                        require_unique=False, require_reviewed=False, min_protein_length=None,
                        cache_root=None, catalogs=None):
        """Add catalogue ORFs with the chosen evidence that this genome encodes.

        ``evidence`` is a rule, an evidence method name, a list of names
        (any of them), or ``None`` (or ``"all"``) for every catalogue ORF. Rules set separate Ribo-seq and MS stringency and combine
        with ``&`` (both) or ``|`` (either), for example
        ``RiboSeq(min_studies=2) & MS(unique=True)``. A name is the simplest rule:
        ``"ms"`` (any supporting MS peptide, the default) is ``MS()``,
        ``"ribosome"`` (or ``"riboseq"``) is ``RiboSeq()``, and ``"proteomics"`` (peptides from
        digested cellular proteins), ``"immunopeptidomics"`` (peptides eluted
        from MHC molecules) or ``"ms_unclassified"`` limit the MS peptide kind. The ``require_*``
        options refine a name as in ``ORFCatalog.query``. The layer records the
        rule. Without ``catalog``, the catalogues covering this species and
        assembly are used, ``catalogs`` choosing among them as in
        :func:`load_catalog`; a species without one raises.

        Each ORF becomes one protein listing every transcript that carries its
        exact spliced coding path. An ORF whose source transcript is gone or
        changed is still added when another transcript carries it, with
        ``source_binding`` saying why; ORFs found on no transcript are excluded.
        ``min_protein_length`` is an explicit selection in amino acids. It does
        not alter the catalogue or the reference CDS proteins already present.
        """
        if self.genome is None:
            raise ValueError("Extra ORFs bind to a genome; build this proteome with from_genome")
        found = extra_orfs(self.genome, evidence, catalog=catalog, context=self.reference,
                           require_resolved=require_resolved, require_unique=require_unique,
                           require_reviewed=require_reviewed, min_protein_length=min_protein_length,
                           cache_root=cache_root, catalogs=catalogs)
        failures = {f.orf_id: f.reason for f in found.source_failures}
        grouped = {}
        for orf in found:
            grouped.setdefault(orf.orf_id, []).append(orf)
        source = found.catalog.source
        annotated = {e.sequence for e in self.proteins if e.kind == "annotated"}
        added = []
        for orf_id, orfs in grouped.items():
            orfs.sort(key=lambda o: (not o.on_source_transcript, _versioned(o.transcript)))
            first = orfs[0]
            added.append(ProteinSequence(
                accession=orf_id, sequence=first.protein_sequence, kind="extra_orf",
                taxon_id=first.model.taxon_id, source=source, orf_type=first.model.orf_type,
                evidence=_evidence_labels(first.evidence),
                evidence_detail=_json(first.evidence.to_dict()),
                occurrences=tuple(Occurrence.from_orf(o) for o in orfs),
                source_transcript_id=first.model.transcript_id,
                source_binding=failures.get(orf_id), equivalent_ids=first.equivalent_ids, orf=first,
            ))
        layer = dict(kind="extra_orf", role="host", source=source, reference=self.reference.key,
                     catalog=_catalog_identity(found.catalog),
                     evidence=None if found.evidence is None else found.evidence.to_dict(),
                     selected=len(added) + len(found.excluded) + len(found.unlifted),
                     added=len(added), excluded=len(found.excluded) + len(found.unlifted),
                     unlifted=len(found.unlifted),
                     not_on_source_transcript=sum(e.source_binding is not None for e in added),
                     identical_to_annotated=sum(e.sequence in annotated for e in added))
        if found.min_protein_length is not None:
            layer["min_protein_length"] = found.min_protein_length
            layer["length_filtered"] = found.length_filtered
        selection = found.catalog._metadata.get("proteome_selection")
        if selection is not None:
            layer["selection"] = selection
        return Proteome(self.reference, self.proteins + tuple(added), self.layers + (layer,),
                        self.excluded + found.excluded + found.unlifted, self.genome,
                        reference_policy=self.reference_policy, records=self.records, links=self.links)

    def with_viral(
        self,
        viral: ProteinSet | ViralReference | str | Iterable[ProteinSet | ViralReference | str],
        *,
        host_taxon_ids: Iterable[int] | None = None,
        host_source: str | None = None,
    ) -> "Proteome":
        """Add curated viral proteins, or a ProteinSet with sourced host taxa."""
        inputs = _viral_inputs(
            viral, self.reference.taxon_id,
            host_taxon_ids=host_taxon_ids, host_source=host_source,
        )
        proteome = self
        for source in inputs:
            proteome = proteome._with_viral_proteins(source)
        return proteome

    def _with_viral_proteins(self, viral: _ViralInput) -> "Proteome":
        proteins = viral.proteins
        if self.reference.taxon_id not in viral.host_taxon_ids:
            raise ValueError(
                f"{proteins.name} (taxon {proteins.taxon_id}) has no recorded host matching "
                f"{self.reference.species} (taxon {self.reference.taxon_id})"
            )
        source = f"viral:{proteins.name}"
        products = {product.protein_id: product for product in viral.products}
        added = tuple(ProteinSequence(accession, sequence, "viral", proteins.taxon_id, source,
                                   viral_gene=products[accession].gene if accession in products else None,
                                   viral_product=products[accession].product if accession in products else None)
                      for accession, sequence in proteins.proteins)
        layer = dict(kind="viral", role="viral", source=source, taxon_id=proteins.taxon_id,
                     host_taxon_ids=sorted(viral.host_taxon_ids), host_source=viral.host_source,
                     provenance=proteins.source, source_sha256=proteins.source_sha256,
                     added=len(added))
        existing_layers = [item for item in self.layers
                           if item.get("role") == "viral" and item.get("source") == source]
        if existing_layers:
            if len(existing_layers) != 1 or dict(existing_layers[0]) != layer:
                raise ValueError(f"Conflicting viral protein set: {proteins.name}")
            existing = self.by_accession()
            if any(entry.accession not in existing
                   or existing[entry.accession].to_dict() != entry.to_dict() for entry in added):
                raise ValueError(f"Conflicting viral protein set: {proteins.name}")
            return self
        collisions = [entry.accession for entry in added
                      if entry.accession in self.by_accession()]
        if collisions:
            raise ValueError("Viral protein accessions already exist: " + ", ".join(collisions[:5]))
        return Proteome(self.reference, self.proteins + added, self.layers + (layer,),
                        self.excluded, self.genome, reference_policy=self.reference_policy,
                        records=self.records, links=self.links)

    def with_derived_proteins(self, proteins):
        """Add supplied sample-specific host sequences linked to parent entries.

        Parent proteins must already be in this proteome. The supplied amino
        acids and change IDs are recorded as provenance; no sequence inference,
        expression evidence or reference measurement is implied.
        """
        proteins = (proteins,) if isinstance(proteins, DerivedProtein) else tuple(proteins)
        if any(not isinstance(item, DerivedProtein) for item in proteins):
            raise TypeError("proteins must contain DerivedProtein records")
        if not proteins:
            return self
        parents = self.by_accession()
        added = []
        for item in proteins:
            missing = [accession for accession in item.parent_accessions if accession not in parents]
            if missing:
                raise ValueError(f"Derived protein {item.accession} has missing parents: {', '.join(missing)}")
            if any(parents[accession].protein_source == "viral" for accession in item.parent_accessions):
                raise ValueError("Derived proteins from viral parents are unsupported")
            selected_occurrences = []
            requested_transcripts = item.parent_transcript_ids or (None,) * len(item.parent_accessions)
            for accession, transcript_id in zip(item.parent_accessions, requested_transcripts):
                occurrences = parents[accession].occurrences
                if transcript_id is None:
                    if len(occurrences) > 1:
                        raise ValueError(f"Derived protein {item.accession} must choose a parent transcript "
                                         f"for {accession}")
                    occurrence = occurrences[0] if occurrences else None
                else:
                    occurrence = next((o for o in occurrences if o.transcript_id == transcript_id), None)
                    if occurrence is None:
                        raise ValueError(f"{transcript_id} is not an occurrence of parent {accession}")
                selected_occurrences.append(occurrence)
            derivation = ProteinDerivation(
                parent_accessions=item.parent_accessions,
                parent_sources=tuple(parents[accession].protein_source
                                     for accession in item.parent_accessions),
                sample_id=item.sample_id, change_ids=item.change_ids,
                parent_gene_ids=tuple(o.gene_id if o is not None else None
                                      for o in selected_occurrences),
                parent_transcript_ids=tuple(o.transcript_id if o is not None else None
                                            for o in selected_occurrences),
                source_sha256=item.source_sha256)
            added.append(ProteinSequence(item.accession, item.sequence, "derived", self.reference.taxon_id,
                                       f"derived:{item.source}", derivation=derivation))
        layer = dict(kind="derived", role="derived", reference=self.reference.key,
                     source="derived", sources=sorted({item.source for item in proteins}),
                     samples=sorted({item.sample_id for item in proteins}), added=len(added))
        return Proteome(self.reference, self.proteins + tuple(added), self.layers + (layer,),
                        self.excluded, self.genome, reference_policy=self.reference_policy,
                        records=self.records, links=self.links)

    def filter(self, predicate, *, description=None):
        """Keep proteins for which ``predicate(protein)`` is true; layers are retained."""
        kept = tuple(e for e in self.proteins if predicate(e))
        layer = dict(kind="filter", role="filter", description=description,
                     removed=len(self.proteins) - len(kept), added=0)
        retained = {entry.accession for entry in kept}
        links = tuple(link for link in self.links
                      if link.target_kind != "protein" or link.target in retained)
        return Proteome(self.reference, kept, self.layers + (layer,), self.excluded, self.genome,
                        reference_policy=self.reference_policy, records=self.records, links=links)

    def host(self):
        """The host's own proteins: annotated and extra ORFs, without viral layers."""
        return self.filter(lambda e: e.protein_source != "viral", description="host proteins")

    def by_sequence(self) -> dict[str, tuple[ProteinSequence, ...]]:
        """Proteins grouped by identical sequence, in first-seen order."""
        groups = {}
        for entry in self.proteins:
            groups.setdefault(entry.sequence, []).append(entry)
        return {sequence: tuple(entries) for sequence, entries in groups.items()}

    def by_accession(self) -> dict[str, ProteinSequence]:
        """Every protein keyed by its unique accession."""
        return {entry.accession: entry for entry in self.proteins}

    def __iter__(self) -> Iterator[ProteinSequence]:
        return iter(self.proteins)

    def __len__(self):
        return len(self.proteins)

    def summary(self):
        """Counts by kind and layer, without sequences."""
        return dict(
            reference=None if self.reference is None else self.reference.to_dict(),
            reference_policy=self.reference_policy, identity=self.identity, proteins=len(self.proteins),
            distinct_sequences=len({e.sequence for e in self.proteins}),
            kinds=dict(Counter(e.kind for e in self.proteins)),
            protein_sources=dict(Counter(e.protein_source for e in self.proteins)),
            derived=sum(e.derivation is not None for e in self.proteins),
            producer_records=len(self.records), producer_links=len(self.links),
            layers=[dict(layer) for layer in self.layers], excluded=len(self.excluded),
        )

    def _content(self):
        content = dict(
            reference=None if self.reference is None else self.reference.to_dict(),
            layers=[dict(layer) for layer in self.layers],
            entries=[e.to_dict() for e in self.proteins],
            excluded=[asdict(e) for e in self.excluded],
        )
        if self.reference_policy != "strict":
            content["reference_policy"] = self.reference_policy
        if self.records or self.links:
            content["producer_graph"] = dict(
                records=[record.to_dict() for record in self.records],
                links=[link.to_dict() for link in self.links])
        return content

    @property
    def identity(self):
        """SHA-256 of the canonical content: sequences, occurrences, layers and exclusions."""
        identity = self.__dict__.get("_identity")
        if identity is None:
            identity = _sha256(self._content())
            self.__dict__["_identity"] = identity
        return identity

    def save(self, path):
        """Write a gzipped JSON snapshot that :meth:`load` reads offline."""
        content = self._content()
        snapshot = dict(schema=SNAPSHOT_SCHEMA, identity=_sha256(content), content=content)
        # No timestamp in the gzip header, so identical proteomes give identical files.
        with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                                                    compresslevel=6) as handle:
            handle.write(json.dumps(snapshot, sort_keys=True, separators=(",", ":")).encode())
        return snapshot["identity"]

    @classmethod
    def load(cls, path, *, _verified=False):
        """Read a snapshot offline; corrupted or inconsistent content raises.

        Entries do not carry ``orf``; use :func:`extra_orfs` for ORF objects.
        Older snapshots are normalized to the current source identities when
        loaded, so a legacy exporter-specific identity may change.
        """
        raw = Path(path).read_bytes()
        snapshot = json.loads(gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw)
        schema = snapshot.get("schema")
        if schema not in ("orfdata.proteome.v1", "orfdata.proteome.v2",
                          "orfdata.proteome.v3", "orfdata.proteome.v4", SNAPSHOT_SCHEMA):
            raise ValueError(f"Unsupported proteome snapshot schema: {schema!r}")
        content = snapshot["content"]
        if _sha256(content) != snapshot.get("identity"):
            raise ValueError("Proteome snapshot content does not match its identity")
        reference = None if content["reference"] is None else ReferenceContext(**content["reference"])
        graph = content.get("producer_graph", {})
        entries = list(content["entries"])
        records = list(graph.get("records", ()))
        links = list(graph.get("links", ()))
        if schema != SNAPSHOT_SCHEMA:
            accession_map = {}
            normalized_entries = []
            for original in entries:
                entry = dict(original)
                if entry.get("producer") is not None:
                    producer = dict(entry["producer"])
                    producer_reference = producer.get("reference")
                    if producer_reference is not None:
                        producer["reference"] = ReferenceContext(**producer_reference)
                    accession = _producer_accession(
                        ProducerInfo(**producer), entry.get("taxon_id")
                    )
                    accession_map[entry["accession"]] = accession
                    entry["accession"] = accession
                normalized_entries.append(entry)
            entries = normalized_entries

            record_map = {}
            normalized_records = []
            for original in records:
                old_fields = dict(
                    namespace=original.get("namespace"), provider=original.get("provider"),
                    exporter=original.get("exporter"), sample=original.get("sample_id"),
                    run=original.get("run_id"), kind=original.get("kind"),
                    native=original.get("native_id"),
                )
                old_key = "record|" + "|".join(
                    f"{name}={quote(str(value), safe='')}"
                    for name, value in old_fields.items() if value is not None
                )
                record = ProducerRecord(**original)
                record_map[old_key] = record.key
                normalized_records.append(record)
            records = normalized_records

            normalized_links = []
            for original in links:
                link = dict(original)
                link["source"] = record_map.get(link["source"], link["source"])
                if link["target_kind"] == "entry":
                    link["target_kind"] = "protein"
                if link["target_kind"] == "protein":
                    link["target"] = accession_map.get(link["target"], link["target"])
                else:
                    link["target"] = record_map.get(link["target"], link["target"])
                normalized_links.append(link)
            links = normalized_links
        proteome = cls(reference,
                       (ProteinSequence.from_dict(e) for e in entries),
                       content["layers"], (BindingFailure(**e) for e in content["excluded"]),
                       reference_policy=content.get("reference_policy", "strict"),
                       records=(record if isinstance(record, ProducerRecord)
                                else ProducerRecord(**record) for record in records),
                       links=(ProducerLink(**link) for link in links))
        # A snapshot from the cache was written by save() and checksummed by DataCache; any other
        # file must also re-serialize to its identity, so non-canonical content is rejected.
        canonical_identity = _sha256(proteome._content())
        if (schema == SNAPSHOT_SCHEMA and not _verified
                and canonical_identity != snapshot["identity"]):
            raise ValueError("Proteome snapshot does not round-trip to its identity")
        proteome.__dict__["_identity"] = (snapshot["identity"] if schema == SNAPSHOT_SCHEMA
                                          else canonical_identity)
        return proteome

    def fasta_index(self, *, deduplicate=False):
        """Map every exported first token to its complete source entry set."""
        if type(deduplicate) is not bool:
            raise ValueError("deduplicate must be True or False")
        groups = (self.by_sequence().values() if deduplicate else
                  ((entry,) for entry in self.proteins))
        index = {}
        for group in groups:
            if deduplicate and len(group) > 1:
                members = sorted(entry.accession for entry in group)
                token = (f"records={len(group)}|aa={len(group[0].sequence)}|"
                         f"key={_sha256(members)[:16]}")
            else:
                token = group[0].fasta_token()
            if token in index:
                raise ValueError(f"Colliding FASTA first token: {token}")
            index[token] = tuple(group)
        return index

    def write_fasta(self, path_or_handle, *, width=60, deduplicate=False):
        """Write one record per source entry, or opt into exact-sequence grouping."""
        if type(width) is not int or width < 1:
            raise ValueError("FASTA width must be a positive integer")
        index = self.fasta_index(deduplicate=deduplicate)
        if hasattr(path_or_handle, "write"):
            return self._write(path_or_handle, width, index)
        with open(path_or_handle, "w") as handle:
            return self._write(handle, width, index)

    def _write(self, handle, width, index):
        for token, group in index.items():
            entry = group[0]
            if len(group) == 1:
                header = entry.fasta_header(token=token)
            else:
                fields = {
                    "records": len(group),
                    "protein_sources": ",".join(sorted({item.protein_source for item in group})),
                    "taxa": ",".join(sorted({"unknown" if item.taxon_id is None
                                               else str(item.taxon_id) for item in group})),
                    "accessions": ",".join(sorted(item.accession for item in group)),
                }
                description = " ".join(
                    f"{name}={_fasta_description_value(value)}" for name, value in fields.items()
                )
                header = f">{token} {description}"
            handle.write(header + "\n")
            for i in range(0, len(entry.sequence), width):
                handle.write(entry.sequence[i : i + width] + "\n")

    def write_annotation(self, directory):
        """Write this proteome's ncORFs as a GTF with matching transcript and protein FASTAs.

        The three files describe exactly what this proteome adds to the
        reference annotation: ``annotation.gtf`` places one coding transcript
        per ncORF on its host gene, ``transcripts.fasta`` holds each ORF's
        spliced coding sequence, and ``proteins.fasta`` its translation. The
        identifiers match, so PyEnsembl reads the three together and varcode,
        isovar, topiary and vaxrank can give a variant inside an ncORF a coding
        consequence. Concatenate the GTF with the assembly's own annotation to
        keep the reference transcripts as well.

        Reference CDS entries are left out, since the reference annotation
        already has them, and viral entries are too: they are not placed on
        this genome. Returns the written paths and counts.
        """
        from pathlib import Path

        from .fasta import write_fasta
        from .gtf import orf_gtf_lines

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        entries = [e for e in self.proteins if e.kind == "extra_orf" and e.orf is not None]
        paths = {name: directory / f"{name}" for name in
                 ("annotation.gtf", "transcripts.fasta", "proteins.fasta")}
        (paths["annotation.gtf"]).write_text(
            f"#!genome-build {self.reference.genome_build if self.reference else 'unknown'}\n"
            f"#!annotation-source orfdata noncanonical ORFs\n"
            + "".join(orf_gtf_lines(entry.orf.record for entry in entries)), encoding="utf-8")
        with paths["transcripts.fasta"].open("w") as transcripts, \
             paths["proteins.fasta"].open("w") as proteins:
            for entry in entries:
                description = (f"gene={entry.orf.record.gene_name or ''} "
                               f"orf_biotype={entry.orf_type or ''} "
                               f"host_transcript={entry.source_transcript_id or ''}").rstrip()
                write_fasta(transcripts, entry.accession, entry.orf.coding_sequence,
                            description=description)
                write_fasta(proteins, f"{entry.accession}.p", entry.sequence,
                            description=description)
        return {"annotation_gtf": str(paths["annotation.gtf"]),
                "transcript_fasta": str(paths["transcripts.fasta"]),
                "protein_fasta": str(paths["proteins.fasta"]),
                "transcript_count": len(entries),
                "protein_count": len(entries),
                "omitted_entries": len(self.proteins) - len(entries)}


def _nc_orf_inputs(reference, sources, cache_root):
    """Resolve one source selection to a catalogue and an optional evidence rule."""
    if sources is None:
        return None, None
    if isinstance(sources, (str, EvidenceRule, ORFCatalog)):
        sources = (sources,)
    if not isinstance(sources, (list, tuple)):
        raise TypeError("nc_orf_sources must be a string, a list, an evidence rule or an ORFCatalog")
    names, rules, supplied = [], [], []
    all_sources = False
    for source in sources:
        if isinstance(source, ORFCatalog):
            supplied.append(source)
        elif isinstance(source, EvidenceRule):
            rules.append(source)
        elif isinstance(source, str) and source:
            if source == "all":
                all_sources = True
            elif source in EVIDENCE_METHODS or source in SYNONYMS:
                rules.append(evidence_rule(source))
            else:
                names.append(source)
        else:
            raise TypeError("Each nc_orf_sources item must be a nonempty string, evidence rule or ORFCatalog")
    if not sources:
        return None, None
    # Resolve explicit names even beside 'all', so a misspelling never disappears.
    named = _species_catalog(reference, None, sorted(set(names)), cache_root) if names else None
    rule = evidence_rule(rules) if rules else None
    registered = (_species_catalog(reference, None, "all", cache_root)
                  if all_sources or rules else named)
    catalogs = ([registered] if registered is not None else []) + supplied
    unique = {}
    for catalog in catalogs:
        unique.setdefault(_catalog_identity(catalog)["sha256"], catalog)
    catalog = _onto(list(unique.values()), reference, cache_root)
    if all_sources or not rules:
        return catalog, None
    if not names and not supplied:
        return catalog, rule

    # A list is a union: named/supplied sources contribute all their records;
    # the rule selects from the remaining sources after evidence is merged.
    from .liftover import unlifted

    complete_sources = set()
    for selected in ([named] if named is not None else []) + supplied:
        records = list(selected.records) + [r for r, _, _ in unlifted(selected.metadata)]
        complete_sources.update(m.source for r in records for m in (r.members or (r,)))

    def accepts(record):
        return (any(m.source in complete_sources for m in (record.members or (record,)))
                or rule.matches(record.evidence))

    metadata = catalog.metadata
    metadata["proteome_selection"] = dict(
        complete_sources=sorted(complete_sources), evidence=rule.to_dict(),
        input=_catalog_identity(catalog))
    _filter_unlifted(metadata, accepts)
    return ORFCatalog([r for r in catalog.records if accepts(r)], metadata), None


def _filter_unlifted(metadata, accepts):
    """Apply a union selection to failed lifts as well as placed records."""
    from .catalog import _record

    if "unresolved" in metadata:
        metadata["unresolved"] = [entry for entry in metadata["unresolved"]
                                  if entry.get("native_record") is None or accepts(_record(entry["native_record"]))]
    lift = metadata.get("liftover")
    if lift:
        _filter_unlifted(lift["input"]["metadata"], accepts)
        lift["unmapped"] = [entry for entry in lift["unmapped"]
                            if entry.get("record") is None or accepts(_record(entry["record"], False))]
    for merged in metadata.get("merged_from", ()):
        _filter_unlifted(merged["metadata"], accepts)


def _species_catalog(reference, catalog, catalogs, cache_root):
    """An explicit catalogue placed onto the reference's assembly, or the species' catalogues ``catalogs`` names."""
    if catalog is not None:
        if catalogs is not None:
            raise ValueError("Pass a catalogue or catalogs, not both")
        return _onto(catalog, reference, cache_root)
    from .catalogs import CatalogRegistry

    return load_catalog(species=reference.species, genome_build=reference.genome_build, catalogs=catalogs,
                        registry=None if cache_root is None else CatalogRegistry(cache_root=cache_root))


def _onto(catalogs, reference, cache_root=None):
    """Explicit catalogues for one genome: any on another assembly of its species lifted to it, then merged."""
    from .liftover import chain_downloaded, chain_routes, lift_to
    from .model import _species

    placed = []
    for catalog in (catalogs,) if isinstance(catalogs, ORFCatalog) else tuple(catalogs):
        if (catalog.genome_build != reference.genome_build and catalog.taxon_id == reference.taxon_id
                and reference.genome_build in chain_routes(_species(catalog.species), catalog.genome_build)):
            route = chain_routes(_species(catalog.species), catalog.genome_build)[reference.genome_build]
            missing = [c.name for c in route if not chain_downloaded(c, cache_root=cache_root)]
            if missing:
                from shlex import join

                command = ["orfdata", "chains"] + ([] if cache_root is None else ["--cache-root", str(cache_root)])
                command += ["download", "--species", _species(catalog.species)]
                command += ["--accept-terms"] if any(c.restricted for c in route) else []
                raise NotInstalledError(f"The {catalog.source} catalogue is on {catalog.genome_build}; lifting it to "
                                        f"{reference.genome_build} needs the chain files {', '.join(missing)}, "
                                        f"which are not installed; download them once with: {join(command)}")
            catalog = lift_to(catalog, reference.genome_build, cache_root=cache_root)
        placed.append(catalog)
    return placed[0] if len(placed) == 1 else merge_catalogs(placed)


def extra_orfs(genome, evidence=None, *, catalog=None, context=None, require_resolved=False,
               require_unique=False, require_reviewed=False, min_protein_length=None,
               cache_root=None, catalogs=None):
    """Place catalogue ORFs on every transcript of a PyEnsembl genome that carries them.

    A transcript carries an ORF when it contains the ORF's exact spliced coding
    path and its sequence translates to the ORF's protein. Without ``catalog``,
    ``load_catalog`` picks the catalogues covering this species and assembly,
    ``catalogs`` choosing among them as there (``"all"`` adds the optional
    ones); a species without one raises. A ``catalog`` (or each of a list) on another
    assembly of the species is lifted to this one through the pinned chains.
    ``cache_root`` is where catalogues and chain files are looked for. ``evidence`` is a rule such as
    ``RiboSeq(min_studies=2) & MS()`` or an evidence method name, as in
    :meth:`Proteome.with_extra_orfs`; ``None`` selects every catalogue ORF.
    ``min_protein_length`` applies a caller-selected amino-acid floor without
    changing the source catalogue.
    """
    _require_installed(genome)
    reference = ReferenceContext.from_genome(genome, context=context)
    catalog = _species_catalog(reference, catalog, catalogs, cache_root)
    if (catalog.taxon_id, catalog.genome_build) != (reference.taxon_id, reference.genome_build):
        lift = ("" if catalog.taxon_id != reference.taxon_id else
                f"; no pinned chain files lift {catalog.genome_build} to {reference.genome_build}, so lift it with "
                f"your own: lift_catalog(catalog, ChainFile(path, source_build={catalog.genome_build!r}, "
                f"target_build={reference.genome_build!r}))")
        raise ValueError(
            f"Catalogue for {catalog.species} {catalog.genome_build} cannot bind to "
            f"{reference.species} {reference.genome_build}{lift}"
        )
    rule = evidence_rule(evidence, require_unique=require_unique, require_reviewed=require_reviewed,
                         require_resolved=require_resolved)
    if min_protein_length is not None:
        _positive_integer(min_protein_length, "min_protein_length")
    selected_records = catalog.query(evidence=rule)
    records = tuple(r for r in selected_records if min_protein_length is None
                    or len(r.protein_sequence) >= min_protein_length)
    spans = {r.orf_id: (r.contig, r.strand, r.coding_blocks[0][0], r.coding_blocks[-1][1]) for r in records}
    try:
        containing = _structure.transcripts_containing(genome, spans)
    except AttributeError:  # No SQLite index to read in bulk: query each locus.
        containing = {}
        for orf_id, (contig, strand, start, end) in spans.items():
            try:
                containing[orf_id] = genome.transcript_ids_at_locus(contig, start, end, strand=strand)
            except ValueError:
                containing[orf_id] = ()
    candidates = {r.orf_id: sorted(set(containing[r.orf_id]).union(
        *(_id_forms(m.transcript_id) for m in r.members or (r,)))) for r in records}
    layouts = _layouts(genome, {i for ids in candidates.values() for i in ids})
    transcripts = {}

    def transcript_by_id(transcript_id):
        if transcript_id not in transcripts:
            try:
                transcripts[transcript_id] = genome.transcript_by_id(transcript_id)
            except ValueError:
                transcripts[transcript_id] = None
        return transcripts[transcript_id]

    bound, failures, excluded = [], [], []
    for record in records:
        model = record.model
        placed, reason, carrying = [], None, None
        source_id = next((i for i in _id_forms(record.transcript_id) if i in layouts), None)
        source = transcript_by_id(source_id) if source_id in layouts else None
        pinned = record.reference_context
        # A record from another release of this assembly is placed by path and protein, even on its own transcript.
        elsewhere = pinned is not None and pinned.key != reference.key and _same_assembly(pinned, reference)
        if source is None:
            reason = "Transcript absent from this annotation"
        else:
            try:
                placed.append(TranscriptORF(record, source, context=reference, any_transcript=elsewhere,
                                            _layout=layouts[source_id]))
                if not placed[-1].on_source_transcript:  # Placed by path on another version of it.
                    reason = "ORF transcript version mismatch"
            except ValueError as error:
                reason = _binding_reason(error)
        for transcript_id in candidates[record.orf_id]:
            layout = layouts.get(transcript_id)
            # The source transcript is tried again by path when strict binding failed, e.g. across releases.
            if layout is None or any(o.transcript.id == transcript_id for o in placed):
                continue
            if _layout_offset(layout, model.strand, model.coding_blocks) is None:
                continue
            transcript = transcript_by_id(transcript_id)
            try:
                placed.append(TranscriptORF(record, transcript, context=reference, any_transcript=True,
                                            _layout=layout))
            except ValueError as error:  # Same path, different sequence or protein.
                carrying = carrying or f"{_binding_reason(error)} on {_versioned(transcript)}, which has the coding path"
        if reason is not None and any(o.on_source_transcript for o in placed):
            reason = None  # Another merged source's own transcript carries it.
        if not placed and carrying and reason == "Transcript absent from this annotation":
            reason = carrying  # Why the transcripts that do carry the path don't encode it.
        if reason is not None:
            failures.append(BindingFailure(record.orf_id, record.transcript_id, reason))
            if not placed:
                excluded.append(failures[-1])
        bound.extend(placed)
    # ORFs a lift left out are reported too, under the same rule, rather than silently missing.
    from .catalog import _merge_key
    from .liftover import unlifted

    left_out, seen = [], set()
    filtered_keys = {_merge_key(record) for record in selected_records if min_protein_length is not None
                     and len(record.protein_sequence) < min_protein_length}
    length_filtered = sum(min_protein_length is not None
                          and len(record.protein_sequence) < min_protein_length
                          for record in selected_records)
    for record, build, reason in unlifted(catalog.metadata):
        if rule is not None and not rule.matches(record.evidence):
            continue
        key = _merge_key(record)
        if min_protein_length is not None and len(record.protein_sequence) < min_protein_length:
            if key not in filtered_keys:
                filtered_keys.add(key)
                length_filtered += 1
            continue
        if key not in seen:
            seen.add(_merge_key(record))  # One ORF two merged catalogues failed to lift counts once.
            left_out.append(BindingFailure(record.orf_id, record.transcript_id, f"Did not lift to {build}: {reason}"))
    left_out = tuple(left_out)
    return ExtraORFs(tuple(bound), tuple(excluded), catalog, reference, rule, tuple(failures), genome, left_out,
                     min_protein_length=min_protein_length, length_filtered=length_filtered)


def _binding_reason(error):
    """Why an ORF doesn't bind, saying so when its start codon reads differently in this sequence."""
    found = re.match(r"Unsupported ORF initiation codon (\w+) for NCBI table (\d+)", str(error))
    if found:  # The catalogue's start is a standard one; this assembly or annotation has another codon there.
        return f"Start codon reads {found.group(1)} (not a start in NCBI table {found.group(2)})"
    return str(error)


def _coding_structures(genome):
    """Exon and CDS structure of every transcript with an annotated protein."""
    try:
        found = [s for s in _structure.structures(genome) if s.protein_id]
    except AttributeError:  # No SQLite index to read in bulk: use transcript objects.
        found = [_structure_from_transcript(t) for t in genome.transcripts() if t.protein_id]
    return found


def _layouts(genome, transcript_ids):
    """Exon layouts by transcript ID, for the IDs this genome has."""
    try:
        return {s.transcript_id: s.layout for s in _structure.structures(genome, transcript_ids=transcript_ids)}
    except AttributeError:
        layouts = {}
        for transcript_id in transcript_ids:
            try:
                layouts[transcript_id] = _structure_from_transcript(genome.transcript_by_id(transcript_id)).layout
            except ValueError:
                continue
        return layouts


def _structure_from_transcript(transcript):
    from .model import _exon_layout, _merge

    exons, offsets, length = _exon_layout(transcript)
    version = getattr(transcript.protein, "protein_version", None) if transcript.protein_id else None
    return _structure.Structure(
        transcript.id, _versioned(transcript), transcript.gene_id, transcript.gene_name or None,
        transcript.biotype, transcript.contig, transcript.strand, tuple(exons), tuple(offsets), length,
        _merge(transcript.coding_sequence_position_ranges) if transcript.protein_id else (),
        transcript.protein_id, version, transcript.contains_start_codon if transcript.protein_id else False,
        transcript.contains_stop_codon if transcript.protein_id else False)


def _cached_snapshot(key, make, cache_root):
    """The snapshot for ``key`` in the DataCache store, built by ``make()`` when absent or invalid."""
    import shutil
    import tempfile
    import time

    from .datasets import _datacache, _digest, get_data_dir

    datacache = _datacache()
    store = get_data_dir(cache_root) / "proteomes" / _digest(key)
    inspected = datacache.inspect_bundle(store)
    if inspected.status == "available":
        return Path(inspected.files[_SNAPSHOT].path)
    if inspected.status not in ("missing", "invalid"):
        raise PermissionError(f"The proteome cache at {store} is {inspected.status}")
    store.parent.mkdir(parents=True, exist_ok=True)
    for stale in store.parent.glob(".partial-*"):  # Left by a build that was killed a day or more ago.
        try:
            if time.time() - stale.stat().st_mtime > 86400:
                shutil.rmtree(stale, ignore_errors=True)
        except OSError:
            pass
    with tempfile.TemporaryDirectory(dir=store.parent, prefix=".partial-") as work:
        files = {_SNAPSHOT: Path(work) / _SNAPSHOT, "recipe.json": Path(work) / "recipe.json"}
        proteome = make()
        proteome.__dict__["_identity"] = proteome.save(files[_SNAPSHOT])
        files["recipe.json"].write_text(json.dumps(dict(key, identity=proteome.identity), indent=1, sort_keys=True))
        assets = {name: dict(path=str(file), size=file.stat().st_size,
                             sha256=hashlib.sha256(file.read_bytes()).hexdigest()) for name, file in files.items()}
        try:
            paths = datacache.install_bundle(store, assets, force=inspected.status == "invalid",
                                             download_options=dict(show_progress=False))
        except datacache.FileValidationError:
            # Another program installed this build meanwhile; snapshots are byte-identical, so use it.
            again = datacache.inspect_bundle(store)
            if again.status != "available":
                raise
            return Path(again.files[_SNAPSHOT].path)
    return Path(paths[_SNAPSHOT])


def _pyensembl_version():
    from pyensembl import __version__

    return __version__


def _annotation_files(genome, *, include_proteins=True):
    """Inputs used by the build: release URLs as given, local files by size and SHA-256."""
    sources = [getattr(genome, "_gtf_path_or_url", None),
               *(getattr(genome, "_transcript_fasta_paths_or_urls", None) or ())]
    if include_proteins:
        sources.extend(getattr(genome, "_protein_fasta_paths_or_urls", None) or ())
    found = []
    for source in (str(s) for s in sources if s is not None):
        path = Path(source)
        if "://" in source or not path.is_file():
            found.append([source])
            continue
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        found.append([path.name, path.stat().st_size, digest.hexdigest()])
    return found


def _json(value):
    """Plain JSON form, so identities do not depend on tuples versus lists."""
    return json.loads(json.dumps(value))


def _sha256(content):
    text = json.dumps(content, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _catalog_identity(catalog):
    identity = getattr(catalog, "_orfdata_identity", None)
    if identity is None:
        identity = dict(source=catalog.source, records=len(catalog.records),
                        sha256=_sha256(catalog.to_dict()))
        catalog._orfdata_identity = identity
    return identity


def _same_id(query, bare, versioned):
    """Exact ID match; a differing Ensembl version is an error, not a fallback.

    ``versioned`` may be a callable, resolved only for a query that could be a
    versioned form of ``bare``.
    """
    if query == bare:
        return True
    if not query.startswith(bare + "."):
        return False
    versioned = versioned() if callable(versioned) else versioned
    if query == versioned:
        return True
    match = _ENSEMBL_VERSIONED.fullmatch(query)
    if match and match.group(1) == bare:
        raise ValueError(f"{query} requested; this reference has {versioned}")
    return False


def _require_installed(genome):
    """Never let a proteome query make PyEnsembl download or index annotation data."""
    inspect = getattr(genome, "inspect_data", None)
    status = None if inspect is None else inspect()["annotation"]
    if status not in (None, "indexed"):
        install = "pyensembl install"
        if getattr(genome, "pinned_name", None) is not None:
            install = f"orfdata genomes install {genome.pinned_name} --force"
        elif getattr(genome, "release", None) is not None and getattr(genome, "species", None) is not None:
            install += f" --release {genome.release} --species {genome.species.latin_name}"
        # Absent, partly downloaded or unindexed annotation is installed by the command; damage is a failure.
        error = NotInstalledError if status in ("missing", "incomplete", "not indexed") else ValueError
        raise error(f"{genome.reference_name} annotation is {status}; install it first with: {install}")


