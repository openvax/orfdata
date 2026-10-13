"""Reference proteomes: annotated proteins plus explicit, scoped additions.

A proteome starts from the proteins annotated in one PyEnsembl genome. Extra
ORFs and viral proteins are separate, opt-in layers. Each layer keeps its
source and selection, and every catalogue ORF that cannot be bound to the
genome is reported rather than dropped silently.
"""

import gzip
import hashlib
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Optional

from .catalog import EVIDENCE_METHODS, load_catalog
from .model import ReferenceORF, _species
from .reference import ReferenceContext

# Evidence labels reported for an included ORF; "mhc" includes human HLA.
_EVIDENCE_LABELS = ("ribosome", "shotgun", "mhc", "ms_unclassified")


@dataclass(frozen=True)
class ProteomeEntry:
    """One protein sequence and where it came from.

    ``kind`` is ``"annotated"``, ``"extra_orf"`` or ``"viral"``. ``biotype`` is
    the transcript biotype for annotated proteins and the source ORF class for
    extra ORFs. ``evidence`` lists the evidence types supporting an extra ORF.
    """

    accession: str
    sequence: str
    kind: str
    taxon_id: int
    source: str
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    transcript_id: Optional[str] = None
    contig: Optional[str] = None
    biotype: Optional[str] = None
    evidence: tuple = ()
    orf: Optional[ReferenceORF] = field(default=None, compare=False, repr=False)

    def fasta_header(self):
        fields = dict(kind=self.kind, taxon=self.taxon_id, source=self.source, gene=self.gene_id,
                      gene_name=self.gene_name, transcript=self.transcript_id, contig=self.contig,
                      biotype=self.biotype, evidence=",".join(self.evidence) or None)
        described = " ".join(f"{k}={str(v).replace(' ', '_')}" for k, v in fields.items() if v is not None)
        return f">{self.accession} {described}"


@dataclass(frozen=True)
class ExcludedORF:
    """A selected catalogue ORF that did not bind to this genome, and why."""

    orf_id: str
    transcript_id: str
    reason: str


@dataclass(frozen=True)
class ExtraORFs:
    """Catalogue ORFs bound to the transcripts of one genome."""

    orfs: tuple
    excluded: tuple
    catalog: object = field(repr=False)
    evidence: Optional[str] = None
    filters: Mapping = field(default_factory=dict)

    def __iter__(self):
        return iter(self.orfs)

    def __len__(self):
        return len(self.orfs)

    def for_transcript(self, transcript):
        """Bound ORFs on one transcript, given as an object or ID."""
        transcript_id = str(getattr(transcript, "id", transcript)).split(".")[0]
        return tuple(o for o in self.orfs if str(o.transcript.id).split(".")[0] == transcript_id)

    def for_gene(self, gene):
        """Bound ORFs whose source gene matches a gene object, ID or name."""
        names = {getattr(gene, "id", gene), getattr(gene, "name", gene)}
        return tuple(o for o in self.orfs if names & {o.model.gene_id, o.model.gene_name})


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
    every addition and ``excluded`` the selected ORFs that did not bind.
    """

    def __init__(self, reference, entries, layers=(), excluded=(), genome=None):
        self.reference = reference
        self.entries = tuple(entries)
        self.layers = tuple(MappingProxyType(dict(layer)) for layer in layers)
        self.excluded = tuple(excluded)
        self.genome = genome
        accessions = Counter(e.accession for e in self.entries)
        duplicated = sorted(a for a, n in accessions.items() if n > 1)
        if duplicated:
            raise ValueError("Duplicate proteome accessions: " + ", ".join(duplicated[:5]))

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
        for transcript in genome.transcripts():
            if not transcript.protein_id:
                continue
            sequence = transcript.protein_sequence
            if not sequence:
                continue
            version = transcript.protein.protein_version
            entries.append(ProteomeEntry(
                accession=transcript.protein_id + ("" if version is None else f".{version}"),
                sequence=sequence, kind="annotated", taxon_id=reference.taxon_id, source=source,
                gene_id=transcript.gene_id, gene_name=transcript.gene_name or None,
                transcript_id=getattr(transcript, "versioned_id", transcript.id),
                contig=transcript.contig, biotype=transcript.biotype,
            ))
        layer = dict(kind="annotated", source=source, reference=reference.key, added=len(entries))
        return cls(reference, entries, (layer,), genome=genome)

    def with_extra_orfs(self, evidence="ms", *, catalog=None, require_resolved=False,
                        require_unique=False, require_reviewed=False):
        """Add catalogue ORFs with the chosen evidence that bind to this genome.

        ``evidence`` is ``"ms"`` (any supporting MS peptide), ``"ribosome"``,
        ``"shotgun"``, ``"mhc"``, ``"hla"``, ``"ms_unclassified"``, or ``None``
        for every catalogue ORF. The ``require_*`` options have the same meaning
        as in ``ORFCatalog.query``. Without ``catalog``, the catalogue covering
        this species and assembly is used; a species without one raises.
        """
        if self.genome is None:
            raise ValueError("Extra ORFs bind to a genome; build this proteome with from_genome")
        found = extra_orfs(self.genome, evidence, catalog=catalog, context=self.reference,
                           require_resolved=require_resolved, require_unique=require_unique,
                           require_reviewed=require_reviewed)
        source = _catalog_source(found.catalog)
        annotated = {e.sequence for e in self.entries if e.kind == "annotated"}
        added = tuple(ProteomeEntry(
            accession=orf.model.orf_id, sequence=orf.protein_sequence, kind="extra_orf",
            taxon_id=orf.model.taxon_id, source=source, gene_id=orf.model.gene_id,
            gene_name=orf.model.gene_name or None,
            transcript_id=getattr(orf.transcript, "versioned_id", orf.transcript.id),
            contig=orf.transcript.contig, biotype=orf.model.orf_type,
            evidence=tuple(m for m in _EVIDENCE_LABELS if orf.evidence.has_evidence(m)), orf=orf,
        ) for orf in found)
        layer = dict(kind="extra_orf", source=source, evidence=evidence, filters=dict(found.filters),
                     selected=len(found.orfs) + len(found.excluded), added=len(added),
                     excluded=len(found.excluded),
                     identical_to_annotated=sum(e.sequence in annotated for e in added))
        return Proteome(self.reference, self.entries + added, self.layers + (layer,),
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
        layer = dict(kind="viral", source=source, taxon_id=viral.taxon_id,
                     host_taxon_ids=sorted(viral.host_taxon_ids), host_source=viral.host_source,
                     provenance=viral.source, source_sha256=viral.source_sha256, added=len(added))
        return Proteome(self.reference, self.entries + added, self.layers + (layer,),
                        self.excluded, self.genome)

    def filter(self, predicate, *, description=None):
        """Keep entries for which ``predicate(entry)`` is true; layers are retained."""
        kept = tuple(e for e in self.entries if predicate(e))
        layer = dict(kind="filter", description=description, removed=len(self.entries) - len(kept),
                     added=0)
        return Proteome(self.reference, kept, self.layers + (layer,), self.excluded, self.genome)

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def summary(self):
        """Counts by kind and layer, without sequences."""
        return dict(
            reference=self.reference.to_dict(), entries=len(self.entries),
            distinct_sequences=len({e.sequence for e in self.entries}),
            kinds=dict(Counter(e.kind for e in self.entries)),
            layers=[dict(layer) for layer in self.layers], excluded=len(self.excluded),
        )

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


def extra_orfs(genome, evidence=None, *, catalog=None, context=None, require_resolved=False,
               require_unique=False, require_reviewed=False):
    """Bind catalogue ORFs to every matching transcript of a PyEnsembl genome.

    Without ``catalog``, the catalogue covering this species and assembly is
    used; a species without one raises. ``evidence=None`` selects every
    catalogue ORF. Selected ORFs whose transcript is absent or disagrees with
    this annotation are returned in ``excluded`` with the reason.
    """
    _require_installed(genome)
    reference = ReferenceContext.from_genome(genome, context=context)
    if catalog is None:
        catalog = _default_catalog(reference)
    if (catalog.taxon_id, catalog.genome_build) != (reference.taxon_id, reference.genome_build):
        raise ValueError(
            f"Catalogue for {catalog.species} {catalog.genome_build} cannot bind to "
            f"{reference.species} {reference.genome_build}"
        )
    if evidence is not None and evidence not in EVIDENCE_METHODS:
        raise ValueError(f"Unknown evidence method: {evidence}")
    filters = {k: v for k, v in dict(require_resolved=require_resolved, require_unique=require_unique,
                                     require_reviewed=require_reviewed).items() if v}
    bound, excluded, transcripts = [], [], {}
    for record in catalog.query(evidence=evidence, **filters):
        transcript_id = record.transcript_id.split(".")[0]
        if transcript_id not in transcripts:
            try:
                transcripts[transcript_id] = genome.transcript_by_id(transcript_id)
            except ValueError:
                transcripts[transcript_id] = None
        transcript = transcripts[transcript_id]
        if transcript is None:
            excluded.append(ExcludedORF(record.orf_id, record.transcript_id,
                                        "Transcript absent from this annotation"))
            continue
        try:
            bound.append(record.bind(transcript, context=reference))
        except ValueError as error:
            excluded.append(ExcludedORF(record.orf_id, record.transcript_id, str(error)))
    return ExtraORFs(tuple(bound), tuple(excluded), catalog, evidence, MappingProxyType(filters))


def _require_installed(genome):
    """Never let a proteome query make PyEnsembl download or index annotation data."""
    inspect = getattr(genome, "inspect_data", None)
    status = None if inspect is None else inspect()["annotation"]
    if status not in (None, "indexed"):
        install = "pyensembl install"
        if getattr(genome, "release", None) is not None and getattr(genome, "species", None) is not None:
            install += f" --release {genome.release} --species {genome.species.latin_name}"
        raise ValueError(f"{genome.reference_name} annotation is {status}; install it first with: {install}")


def _default_catalog(reference):
    if reference.taxon_id == 9606 and reference.genome_build == "GRCh38":
        return load_catalog()
    raise ValueError(
        f"No extra-ORF catalogue is available for {_species(reference.species)} "
        f"{reference.genome_build}; pass catalog= explicitly (see CatalogRegistry)"
    )


def _catalog_source(catalog):
    metadata = catalog.metadata
    manifest = metadata.get("catalog_manifest")
    if manifest:
        return f"{manifest['catalog_id']}:{manifest['version']}"
    return f"orfdata-catalog:{metadata.get('dataset_version', 'unversioned')}"
