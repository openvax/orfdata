"""Reference proteomes: annotated proteins plus explicit, scoped additions.

A proteome starts from the proteins annotated in one PyEnsembl genome. Extra
ORFs and viral proteins are separate, opt-in layers. Each layer keeps its
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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Optional

from . import _structure
from .catalog import ORFCatalog, load_catalog, merge_catalogs
from .evidence import EvidenceRule, evidence_rule
from .model import TranscriptORF, _contig, _id_forms, _layout_offset, _same_assembly, _versioned
from .reference import ReferenceContext

SNAPSHOT_SCHEMA = "orfdata.proteome.v1"
_SNAPSHOT = "proteome.json.gz"
# PyEnsembl's convention: only ENS IDs carry a ".N" version; other dots are part of the ID.
_ENSEMBL_VERSIONED = re.compile(r"(ENS.*)\.(0|[1-9][0-9]*)")
# Evidence labels reported for an included ORF; "mhc" includes human HLA.
_EVIDENCE_LABELS = ("ribosome", "shotgun", "mhc", "ms_unclassified")


def _evidence_labels(evidence):
    """Evidence kinds for a FASTA header; "ms" for MS support whose assay the source doesn't name."""
    labels = [m for m in _EVIDENCE_LABELS if evidence.has_evidence(m)]
    if evidence.has_evidence("ms") and not set(labels) - {"ribosome"}:
        labels.append("ms")
    return tuple(labels)


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
class ProteomeEntry:
    """One protein: its sequence, identity, source and every transcript encoding it.

    ``kind`` is ``"annotated"``, ``"extra_orf"`` or ``"viral"``. ``accession``
    is the versioned protein ID, the source ORF ID or the viral protein ID.
    ``occurrences`` lists every transcript of this genome encoding the protein,
    the source's own transcript first; viral proteins have none. ``gene_id``,
    ``gene_name``, ``transcript_id``, ``contig`` and ``transcript_biotype``
    describe the first occurrence. ``source_binding`` explains why an extra ORF
    could not be placed on its source's own transcript; it is None when it was.
    ``orf_type`` is an extra ORF's source category, kept verbatim, and
    ``evidence`` lists its supporting evidence types; ``evidence_detail`` holds
    the source's evidence records. ``equivalent_ids`` lists other sources'
    ``source:id`` for the same extra ORF when catalogues were merged.
    """

    accession: str
    sequence: str
    kind: str
    taxon_id: int
    source: str
    orf_type: Optional[str] = None
    evidence: tuple = ()
    evidence_detail: Optional[Mapping] = field(default=None, compare=False)
    occurrences: tuple = ()
    source_transcript_id: Optional[str] = None
    source_binding: Optional[str] = None
    equivalent_ids: tuple = ()
    orf: Optional[TranscriptORF] = field(default=None, compare=False, repr=False)

    @property
    def _first(self):
        return self.occurrences[0] if self.occurrences else None

    @property
    def gene_id(self):
        return self._first and self._first.gene_id

    @property
    def gene_name(self):
        return self._first and self._first.gene_name

    @property
    def transcript_id(self):
        return self._first and self._first.transcript_id

    @property
    def contig(self):
        return self._first and self._first.contig

    @property
    def transcript_biotype(self):
        return self._first and self._first.transcript_biotype

    def fasta_header(self):
        first = self._first
        fields = dict(kind=self.kind, taxon=self.taxon_id, source=self.source, gene=self.gene_id,
                      gene_name=self.gene_name, transcript=self.transcript_id, contig=self.contig,
                      transcript_biotype=self.transcript_biotype, orf_type=self.orf_type,
                      evidence=",".join(self.evidence) or None,
                      source_transcript=self.source_transcript_id
                      if first is not None and not first.on_source_transcript else None,
                      equivalent=",".join(self.equivalent_ids) or None,
                      transcripts=len(self.occurrences) if len(self.occurrences) > 1 else None)
        described = " ".join(f"{k}={str(v).replace(' ', '_')}" for k, v in fields.items() if v is not None)
        return f">{self.accession} {described}"

    def to_dict(self):
        return dict(accession=self.accession, sequence=self.sequence, kind=self.kind,
                    taxon_id=self.taxon_id, source=self.source, orf_type=self.orf_type,
                    evidence=list(self.evidence), evidence_detail=self.evidence_detail,
                    occurrences=[[getattr(o, name) for name in _OCCURRENCE_FIELDS] for o in self.occurrences],
                    source_transcript_id=self.source_transcript_id, source_binding=self.source_binding,
                    **({"equivalent_ids": list(self.equivalent_ids)} if self.equivalent_ids else {}))

    @classmethod
    def from_dict(cls, record):
        record = dict(record)
        record["evidence"] = tuple(record.get("evidence", ()))
        record["equivalent_ids"] = tuple(record.get("equivalent_ids", ()))
        blocks = _OCCURRENCE_FIELDS.index("coding_blocks")
        record["occurrences"] = tuple(
            Occurrence(*values[:blocks], tuple(map(tuple, values[blocks])), *values[blocks + 1 :])
            for values in record.get("occurrences", ()))
        return cls(**record)


@dataclass(frozen=True)
class ExtraORFs:
    """Catalogue ORFs on the transcripts of one genome and reference.

    Each ORF appears once for every transcript carrying its exact spliced
    coding path and protein; ``orf.on_source_transcript`` marks the transcript
    its source names. ``source_failures`` lists ORFs whose source transcript
    did not carry them, and ``excluded`` those found on no transcript at all.
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

    def __post_init__(self):
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("ExtraORFs requires its resolved ReferenceContext")
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

    def for_orf(self, orf_id):
        """Every transcript placement of one ORF, the source's own transcript first.

        After merging catalogues, any source's ID finds the ORF, bare or as
        ``source:id``; a bare ID two sources use for different ORFs raises.
        """
        if self.catalog is not None:
            try:
                orf_id = self.catalog.get(orf_id).orf_id
            except KeyError:
                return ()
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
class ViralProteome:
    """Proteins of one virus, usable only with its recorded host species.

    ``host_taxon_ids`` are NCBI species taxa, and ``host_source`` says where
    that host range comes from. Adding the proteins to a proteome for any other
    species fails; no host is inferred from orthology or a virus name.
    """

    name: str
    taxon_id: int
    host_taxon_ids: frozenset
    proteins: tuple
    source: str
    host_source: str
    source_sha256: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "host_taxon_ids", frozenset(self.host_taxon_ids))
        object.__setattr__(self, "proteins", tuple(tuple(p) for p in self.proteins))
        if type(self.taxon_id) is not int or self.taxon_id < 1:
            raise ValueError("A viral proteome requires a positive NCBI taxon_id")
        if not self.host_taxon_ids or any(type(t) is not int or t < 1 for t in self.host_taxon_ids):
            raise ValueError("A viral proteome requires explicit NCBI host taxon IDs")
        if not all(isinstance(x, str) and x.strip() for x in (self.name, self.source, self.host_source)):
            raise ValueError("A viral proteome requires a name, source and host-range source")
        if not self.proteins:
            raise ValueError("A viral proteome requires proteins")
        for accession, sequence in self.proteins:
            if not accession or not sequence or not sequence.isalpha() or not sequence.isupper():
                raise ValueError(f"Invalid viral protein record: {accession!r}")
        if len({a for a, _ in self.proteins}) != len(self.proteins):
            raise ValueError("Duplicate viral protein accessions")

    @classmethod
    def from_fasta(cls, path, *, name, taxon_id, host_taxon_ids, host_source, source=None):
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
        return cls(name, taxon_id, host_taxon_ids, tuple(proteins), source or Path(path).name,
                   host_source, hashlib.sha256(raw).hexdigest())


class Proteome:
    """Protein sequences for one species and genome, built in explicit layers.

    Start with :meth:`from_genome`; add layers with :meth:`with_extra_orfs` and
    :meth:`with_viral`. Each method returns a new proteome. ``layers`` records
    every addition with its role (``"host"``, ``"viral"`` or ``"filter"``),
    source and selection; ``excluded`` lists selected ORFs found on no
    transcript. :meth:`save` writes a checksummed snapshot that :meth:`load`
    reads back offline, and ``identity`` changes whenever the content does.
    """

    def __init__(self, reference, entries, layers=(), excluded=(), genome=None):
        if not isinstance(reference, ReferenceContext):
            raise TypeError("A proteome requires a ReferenceContext")
        self.reference = reference
        self.entries = tuple(entries)
        self.layers = tuple(MappingProxyType(_json(dict(layer))) for layer in layers)
        self.excluded = tuple(excluded)
        self.genome = genome
        accessions = Counter(e.accession for e in self.entries)
        duplicated = sorted(a for a, n in accessions.items() if n > 1)
        if duplicated:
            raise ValueError("Duplicate proteome accessions: " + ", ".join(duplicated[:5]))
        hosts = {layer["source"]: layer.get("host_taxon_ids", ()) for layer in self.layers
                 if layer["role"] == "viral"}
        for entry in self.entries:
            if entry.kind == "viral":
                if reference.taxon_id not in hosts.get(entry.source, ()):
                    raise ValueError(f"Viral protein {entry.accession} has no host layer for this species")
            elif entry.taxon_id != reference.taxon_id:
                raise ValueError(f"Host protein {entry.accession} is from another species")
        for layer in self.layers:
            if layer.get("reference") not in (None, reference.key):
                raise ValueError("A layer declares a different reference")

    @classmethod
    def from_genome(cls, genome, *, context=None):
        """Every annotated protein with a sequence in this PyEnsembl genome.

        This matches the annotation's protein FASTA, including all transcript
        biotypes with a translation and proteins on patch or haplotype contigs.
        Use :meth:`filter` to restrict them, for example by ``contig``.
        """
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
            entries.append(ProteomeEntry(
                accession=accession, sequence=sequence, kind="annotated", taxon_id=reference.taxon_id,
                source=source, occurrences=(occurrence,), source_transcript_id=structure.versioned_id))
        entries.sort(key=lambda e: e.source_transcript_id)
        layer = dict(kind="annotated", role="host", source=source, reference=reference.key,
                     added=len(entries))
        return cls(reference, entries, (layer,), genome=genome)

    @classmethod
    def build(cls, genome, evidence="ms", *, catalog=None, viral=(), extra_orfs=True, context=None,
              cache=False, cache_root=None):
        """The usual reference proteome in one call, optionally kept in a cache.

        That is :meth:`from_genome`, then :meth:`with_extra_orfs` with
        ``evidence`` and ``catalog`` unless ``extra_orfs=False``, then
        :meth:`with_viral` for each of ``viral``. With ``cache=True`` the result
        is saved under the OpenVax DataCache root (``ORFDATA_CACHE``,
        ``OPENVAX_DATA_CACHE`` or ``cache_root``), keyed by the reference, the
        evidence rule, the catalogue and viral identities, the snapshot schema
        and the orfdata version. A later identical build loads that verified
        snapshot instead of rebuilding. ``cache_root`` is also where the
        species' registered catalogue is looked for. Nothing is downloaded.
        """
        _require_installed(genome)
        reference = ReferenceContext.from_genome(genome, context=context)
        viral = tuple(viral)
        if not extra_orfs and catalog is not None:
            raise ValueError("A catalogue needs extra_orfs=True")
        rule = None
        if extra_orfs:
            rule = evidence_rule(evidence)
            if catalog is None:
                from .catalogs import CatalogRegistry

                catalog = load_catalog(species=reference.species, genome_build=reference.genome_build,
                                       registry=None if cache_root is None else CatalogRegistry(cache_root=cache_root))
            else:
                catalog = _onto(catalog, reference, cache_root)
        built = []

        def make():
            proteome = cls.from_genome(genome, context=context)
            if extra_orfs:
                proteome = proteome.with_extra_orfs(rule, catalog=catalog)
            for virus in viral:
                proteome = proteome.with_viral(virus)
            built.append(proteome)
            return proteome

        if not cache:
            return make()
        from .datasets import _digest
        from .version import __version__

        key = dict(schema=SNAPSHOT_SCHEMA, orfdata=__version__, pyensembl=_pyensembl_version(),
                   reference=reference.to_dict(), annotation=_annotation_files(genome),
                   extra_orfs=dict(evidence=None if rule is None else rule.to_dict(),
                                   catalog=_catalog_identity(catalog)) if extra_orfs else None,
                   viral=[_digest(_json(dict(asdict(virus), host_taxon_ids=sorted(virus.host_taxon_ids))))
                          for virus in viral])
        path = _cached_snapshot(key, make, cache_root)
        if built:
            return built[0]
        loaded = cls.load(path, _verified=True)
        proteome = cls(loaded.reference, loaded.entries, loaded.layers, loaded.excluded, genome)
        proteome.__dict__["_identity"] = loaded.identity
        return proteome

    def with_extra_orfs(self, evidence="ms", *, catalog=None, require_resolved=False,
                        require_unique=False, require_reviewed=False, cache_root=None):
        """Add catalogue ORFs with the chosen evidence that this genome encodes.

        ``evidence`` is a rule, an evidence method name, or ``None`` for every
        catalogue ORF. Rules set separate Ribo-seq and MS stringency and combine
        with ``&`` (both) or ``|`` (either), for example
        ``RiboSeq(min_studies=2) & MS(unique=True)``. A name is the simplest rule:
        ``"ms"`` (any supporting MS peptide, the default) is ``MS()``,
        ``"ribosome"`` is ``RiboSeq()``, and ``"shotgun"``, ``"mhc"``, ``"hla"``
        or ``"ms_unclassified"`` limit the MS peptide kind. The ``require_*``
        options refine a name as in ``ORFCatalog.query``. The layer records the
        rule. Without ``catalog``, the catalogue covering this species and
        assembly is used; a species without one raises.

        Each ORF becomes one entry listing every transcript that carries its
        exact spliced coding path. An ORF whose source transcript is gone or
        changed is still added when another transcript carries it, with
        ``source_binding`` saying why; ORFs found on no transcript are excluded.
        """
        if self.genome is None:
            raise ValueError("Extra ORFs bind to a genome; build this proteome with from_genome")
        found = extra_orfs(self.genome, evidence, catalog=catalog, context=self.reference,
                           require_resolved=require_resolved, require_unique=require_unique,
                           require_reviewed=require_reviewed, cache_root=cache_root)
        failures = {f.orf_id: f.reason for f in found.source_failures}
        grouped = {}
        for orf in found:
            grouped.setdefault(orf.orf_id, []).append(orf)
        source = found.catalog.source
        annotated = {e.sequence for e in self.entries if e.kind == "annotated"}
        added = []
        for orf_id, orfs in grouped.items():
            orfs.sort(key=lambda o: (not o.on_source_transcript, _versioned(o.transcript)))
            first = orfs[0]
            added.append(ProteomeEntry(
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
                     selected=len(added) + len(found.excluded),
                     added=len(added), excluded=len(found.excluded),
                     not_on_source_transcript=sum(e.source_binding is not None for e in added),
                     identical_to_annotated=sum(e.sequence in annotated for e in added))
        return Proteome(self.reference, self.entries + tuple(added), self.layers + (layer,),
                        self.excluded + found.excluded, self.genome)

    def with_viral(self, viral):
        """Add a viral proteome whose recorded hosts include this species."""
        if not isinstance(viral, ViralProteome):
            raise TypeError("Expected a ViralProteome")
        if self.reference.taxon_id not in viral.host_taxon_ids:
            raise ValueError(
                f"{viral.name} (taxon {viral.taxon_id}) has no recorded host matching "
                f"{self.reference.species} (taxon {self.reference.taxon_id})"
            )
        source = f"viral:{viral.name}"
        added = tuple(ProteomeEntry(accession, sequence, "viral", viral.taxon_id, source)
                      for accession, sequence in viral.proteins)
        layer = dict(kind="viral", role="viral", source=source, taxon_id=viral.taxon_id,
                     host_taxon_ids=sorted(viral.host_taxon_ids), host_source=viral.host_source,
                     provenance=viral.source, source_sha256=viral.source_sha256, added=len(added))
        return Proteome(self.reference, self.entries + added, self.layers + (layer,),
                        self.excluded, self.genome)

    def filter(self, predicate, *, description=None):
        """Keep entries for which ``predicate(entry)`` is true; layers are retained."""
        kept = tuple(e for e in self.entries if predicate(e))
        layer = dict(kind="filter", role="filter", description=description,
                     removed=len(self.entries) - len(kept), added=0)
        return Proteome(self.reference, kept, self.layers + (layer,), self.excluded, self.genome)

    def host(self):
        """The host's own proteins: annotated and extra ORFs, without viral layers."""
        return self.filter(lambda e: e.kind != "viral", description="host proteins")

    def by_sequence(self):
        """Entries grouped by identical protein sequence, in first-seen order."""
        groups = {}
        for entry in self.entries:
            groups.setdefault(entry.sequence, []).append(entry)
        return {sequence: tuple(entries) for sequence, entries in groups.items()}

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def summary(self):
        """Counts by kind and layer, without sequences."""
        return dict(
            reference=self.reference.to_dict(), identity=self.identity, entries=len(self.entries),
            distinct_sequences=len({e.sequence for e in self.entries}),
            kinds=dict(Counter(e.kind for e in self.entries)),
            layers=[dict(layer) for layer in self.layers], excluded=len(self.excluded),
        )

    def _content(self):
        return dict(
            reference=self.reference.to_dict(),
            layers=[dict(layer) for layer in self.layers],
            entries=[e.to_dict() for e in self.entries],
            excluded=[asdict(e) for e in self.excluded],
        )

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
        """
        raw = Path(path).read_bytes()
        snapshot = json.loads(gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw)
        if snapshot.get("schema") != SNAPSHOT_SCHEMA:
            raise ValueError(f"Unsupported proteome snapshot schema: {snapshot.get('schema')!r}")
        content = snapshot["content"]
        if _sha256(content) != snapshot.get("identity"):
            raise ValueError("Proteome snapshot content does not match its identity")
        proteome = cls(ReferenceContext(**content["reference"]),
                       (ProteomeEntry.from_dict(e) for e in content["entries"]),
                       content["layers"], (BindingFailure(**e) for e in content["excluded"]))
        # A snapshot from the cache was written by save() and checksummed by DataCache; any other
        # file must also re-serialize to its identity, so non-canonical content is rejected.
        if not _verified and _sha256(proteome._content()) != snapshot["identity"]:
            raise ValueError("Proteome snapshot does not round-trip to its identity")
        proteome.__dict__["_identity"] = snapshot["identity"]
        return proteome

    def write_fasta(self, path_or_handle, *, width=60):
        """Write every entry as FASTA with a descriptive header."""
        if hasattr(path_or_handle, "write"):
            return self._write(path_or_handle, width)
        with open(path_or_handle, "w") as handle:
            return self._write(handle, width)

    def _write(self, handle, width):
        for entry in self.entries:
            handle.write(entry.fasta_header() + "\n")
            for i in range(0, len(entry.sequence), width):
                handle.write(entry.sequence[i : i + width] + "\n")


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
                accept = " --accept-terms" if any(c.restricted for c in route) else ""
                raise ValueError(f"The {catalog.source} catalogue is on {catalog.genome_build}; lifting it to "
                                 f"{reference.genome_build} needs the chain files {', '.join(missing)}: download "
                                 f"them once with 'orfdata chains download --species "
                                 f"{_species(catalog.species)}{accept}'")
            catalog = lift_to(catalog, reference.genome_build, cache_root=cache_root)
        placed.append(catalog)
    return placed[0] if len(placed) == 1 else merge_catalogs(placed)


def extra_orfs(genome, evidence=None, *, catalog=None, context=None, require_resolved=False,
               require_unique=False, require_reviewed=False, cache_root=None):
    """Place catalogue ORFs on every transcript of a PyEnsembl genome that carries them.

    A transcript carries an ORF when it contains the ORF's exact spliced coding
    path and its sequence translates to the ORF's protein. Without ``catalog``,
    ``load_catalog`` picks the catalogue covering this species and assembly; a
    species without one raises. A ``catalog`` (or each of a list) on another
    assembly of the species is lifted to this one through the pinned chains.
    ``cache_root`` is where catalogues and chain files are looked for. ``evidence`` is a rule such as
    ``RiboSeq(min_studies=2) & MS()`` or an evidence method name, as in
    :meth:`Proteome.with_extra_orfs`; ``None`` selects every catalogue ORF.
    """
    _require_installed(genome)
    reference = ReferenceContext.from_genome(genome, context=context)
    if catalog is None:
        from .catalogs import CatalogRegistry

        catalog = load_catalog(species=reference.species, genome_build=reference.genome_build,
                               registry=None if cache_root is None else CatalogRegistry(cache_root=cache_root))
    else:
        catalog = _onto(catalog, reference, cache_root)
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
    records = catalog.query(evidence=rule)
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
        placed, reason = [], None
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
                reason = str(error)
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
            except ValueError:
                continue  # Same path, different sequence or protein.
        if reason is not None and any(o.on_source_transcript for o in placed):
            reason = None  # Another merged source's own transcript carries it.
        if reason is not None:
            failures.append(BindingFailure(record.orf_id, record.transcript_id, reason))
            if not placed:
                excluded.append(failures[-1])
        bound.extend(placed)
    # ORFs a lift left out are reported too, under the same rule, rather than silently missing.
    from .liftover import unlifted

    for record, build, reason in unlifted(catalog.metadata):
        if rule is None or rule.matches(record.evidence):
            failures.append(BindingFailure(record.orf_id, record.transcript_id, f"Did not lift to {build}: {reason}"))
            excluded.append(failures[-1])
    return ExtraORFs(tuple(bound), tuple(excluded), catalog, reference, rule, tuple(failures), genome)


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
    try:
        from pyensembl import __version__
    except ImportError:
        return None
    return __version__


def _annotation_files(genome):
    """The genome's annotation sources: release URLs as given, local files by size and SHA-256."""
    sources = [getattr(genome, "_gtf_path_or_url", None),
               *(getattr(genome, "_transcript_fasta_paths_or_urls", None) or ()),
               *(getattr(genome, "_protein_fasta_paths_or_urls", None) or ())]
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
        if getattr(genome, "release", None) is not None and getattr(genome, "species", None) is not None:
            install += f" --release {genome.release} --species {genome.species.latin_name}"
        raise ValueError(f"{genome.reference_name} annotation is {status}; install it first with: {install}")


