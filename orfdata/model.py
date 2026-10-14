"""Reference ORF identity and its binding to an annotated transcript."""

import hashlib
from dataclasses import asdict, dataclass, replace
from typing import Optional

from .initiation_context import annotate_initiation_context
from .reference import ReferenceContext, canonical_species, species_taxon
from .sequence import genetic_code


@dataclass(frozen=True)
class ORFModel:
    """A source-defined ORF; genomic blocks include its terminal stop codon.

    Coordinates are 1-based inclusive; blocks are in ascending genomic order.
    ``orf_type`` describes this reading frame, never the parent RNA's biotype.
    ``translation_table`` is an NCBI genetic-code ID, not inferred from species
    or a contig name. Binding supports tables 1, 2 and 3 (standard, vertebrate
    and yeast mitochondrial); other declared codes fail when bound.
    """

    orf_id: str
    transcript_id: str
    gene_id: str
    gene_name: str
    contig: str
    strand: str
    coding_blocks: tuple
    protein_sequence: str
    orf_type: str
    genome_build: str
    taxon_id: int
    species: str
    translation_table: int = 1
    source: str = ""
    reference_context: Optional[ReferenceContext] = None

    def __post_init__(self):
        if type(self.taxon_id) is not int or self.taxon_id < 1 or not self.species:
            raise ValueError("An explicit species and positive NCBI taxon_id are required")
        known = species_taxon(self.species)
        if known is not None and known != self.taxon_id:
            raise ValueError("ORF species name and taxon_id disagree")
        if self.reference_context is not None:
            if not isinstance(self.reference_context, ReferenceContext):
                raise TypeError("reference_context must be a ReferenceContext")
            if (self.reference_context.taxon_id != self.taxon_id
                    or self.reference_context.genome_build != self.genome_build
                    or canonical_species(self.reference_context.species) != canonical_species(self.species)):
                raise ValueError("ORF and reference context disagree")
        if type(self.translation_table) is not int or self.translation_table < 1:
            raise ValueError("An explicit NCBI translation_table is required")
        if not all((self.orf_id, self.transcript_id, self.genome_build, self.orf_type)):
            raise ValueError("ORF, transcript, assembly and ORF type are required")
        if self.strand not in ("+", "-"):
            raise ValueError("Invalid ORF strand")
        blocks = tuple(tuple(b) for b in self.coding_blocks)
        object.__setattr__(self, "coding_blocks", blocks)
        if (
            not blocks
            or not self.protein_sequence
            or "*" in self.protein_sequence
            or any(type(x) is not int for b in blocks for x in b)
            or any(a < 1 or b < a for a, b in blocks)
            or any(b >= c for (_, b), (c, _) in zip(blocks, blocks[1:]))
            or sum(b - a + 1 for a, b in blocks) != 3 * (len(self.protein_sequence) + 1)
        ):
            raise ValueError("A complete, stopped ORF with valid spliced coding blocks is required")

    @property
    def genomic_positions(self):
        positions = tuple(p for a, b in self.coding_blocks for p in range(a, b + 1))
        return positions[::-1] if self.strand == "-" else positions

    def to_dict(self):
        result = asdict(self)
        if self.reference_context is None:
            result.pop("reference_context")
        return result


def transcript_positions(transcript):
    """Map spliced cDNA offsets to forward-genomic positions on either strand."""
    reverse = transcript.strand == "-"
    exons = sorted(transcript.exons, key=lambda e: e.start, reverse=reverse)
    return tuple(
        p
        for e in exons
        for p in (range(e.end, e.start - 1, -1) if reverse else range(e.start, e.end + 1))
    )


def _exon_layout(transcript):
    """Exon (start, end) pairs in transcript order, their spliced offsets and the RNA length."""
    exons = sorted(((e.start, e.end) for e in transcript.exons), reverse=transcript.strand == "-")
    offsets, total = [], 0
    for start, end in exons:
        offsets.append(total)
        total += end - start + 1
    return exons, offsets, total


def _transcript_match(identifier, transcript):
    """How a source's transcript ID names a transcript.

    ``"id"`` for its ID, ``"versioned"`` for its ID with its version,
    ``"other_version"`` for its ID with another numeric version, else None. IDs
    that contain dots, such as C. elegans ``Y74C9A.3.1``, are compared whole.
    """
    identifier, own = str(identifier), str(transcript.id)
    if identifier == own:
        return "id"
    base, _, version = identifier.rpartition(".")
    if base == own and version.isdigit():
        return "versioned" if identifier == _versioned(transcript) else "other_version"
    return None


def _id_forms(identifier):
    """The transcript IDs a source's ID could be: itself, or without a numeric version."""
    base, _, version = str(identifier).rpartition(".")
    return (str(identifier), base) if base and version.isdigit() else (str(identifier),)


def _layout_position(layout, strand, offset):
    """The forward-genomic position of one spliced offset."""
    exons, offsets, _ = layout
    for (a, b), first in zip(exons, offsets):
        if first <= offset <= first + b - a:
            return a + offset - first if strand == "+" else b - (offset - first)
    raise ValueError(f"Offset {offset} is outside the transcript")


def path_offset(transcript, blocks, *, layout=None):
    """Where a genomic coding path starts in a transcript's spliced sequence.

    ``blocks`` are 1-based inclusive genomic blocks in ascending order. Return
    the 0-based offset into ``transcript.sequence`` at which the path's first
    base lies, or None unless the path occupies consecutive spliced positions:
    each internal junction of the path must be an exon junction of the
    transcript. The coding sequence is then ``sequence[offset:offset + length]``.
    """
    return _layout_offset(layout or _exon_layout(transcript), transcript.strand, blocks)


def _layout_offset(layout, strand, blocks):
    exons, offsets, _ = layout
    forward = strand == "+"
    path = tuple(blocks) if forward else tuple(reversed(blocks))
    if not path:
        return None
    first = path[0][0] if forward else path[0][1]
    index = next((i for i, (a, b) in enumerate(exons) if a <= first <= b), None)
    if index is None:
        return None
    a, b = exons[index]
    offset = offsets[index] + (first - a if forward else b - first)
    for k, (low, high) in enumerate(path):
        a, b = exons[index]
        if not a <= low <= high <= b:
            return None
        if k + 1 == len(path):
            break
        if index + 1 == len(exons):
            return None
        next_low, next_high = path[k + 1]
        following = exons[index + 1]
        # The path may leave an exon only at its splice site, entering the next exon at its start.
        if forward and (high != b or next_low != following[0]):
            return None
        if not forward and (low != a or next_high != following[1]):
            return None
        index += 1
    return offset


def _contig(value):
    """Compare contig names as PyEnsembl normalizes them, e.g. Mito as MITO."""
    value = str(value).removeprefix("chr")
    if value.isalpha():
        value = value.upper()
    return "MT" if value == "M" else value


def _species(value):
    return canonical_species(value)


@dataclass(frozen=True, init=False)
class TranscriptORF:
    """One ORF on one PyEnsembl transcript.

    ``kind`` is ``"annotated"`` for an annotated CDS, or ``"extra_orf"`` for a
    catalogue ORF validated against the transcript's spliced path, sequence
    and protein. ``start`` and ``end`` are Python slice coordinates into the
    transcript's spliced sequence: ``transcript.sequence[start:end]`` is the
    coding sequence, including the stop codon when one is present.

    An ORF can lie on several transcripts with the identical spliced coding
    path. ``source_transcript_id`` names the transcript its source defines it
    on, and ``on_source_transcript`` says whether this is that transcript.
    Another carrying transcript is not evidence that the ORF is translated from
    it. ``gene_id`` and ``gene_name`` describe this transcript's gene in the
    loaded annotation; ``source_gene_id`` and ``source_gene_name`` are the
    source's. ``equivalent_ids`` lists other annotated proteins with the
    identical coding placement and sequence. The transcript, including its
    biotype and annotated protein, is unchanged, and one transcript can carry
    several ORFs. No downloads occur.
    """

    record: object
    transcript: object
    start: int
    end: int
    coding_sequence: str
    context: ReferenceContext
    kind: str = "extra_orf"
    start_complete: bool = True
    stop_complete: bool = True
    annotation_id: Optional[str] = None
    source_transcript_id: Optional[str] = None
    equivalent_ids: tuple = ()

    def __init__(self, record, transcript, *, context=None, any_transcript=False, _layout=None):
        model = record.model
        genome = transcript.genome
        if genome.reference_name != model.genome_build:
            raise ValueError("ORF transcript assembly mismatch")
        species = getattr(genome, "species", None)
        if species is not None and _species(species) != _species(model.species):
            raise ValueError("ORF transcript species mismatch")
        taxon = getattr(genome, "taxon_id", None)
        if taxon is not None and taxon != model.taxon_id:
            raise ValueError("ORF transcript taxon mismatch")
        resolved_context = ReferenceContext.from_genome(genome, context=context)
        if resolved_context.taxon_id != model.taxon_id:
            raise ValueError("ORF transcript taxon mismatch")
        if model.reference_context is not None and model.reference_context.key != resolved_context.key:
            raise ValueError("ORF transcript reference context mismatch (annotation, strain or assembly)")
        if not any_transcript:
            match = _transcript_match(model.transcript_id, transcript)
            if match is None:
                raise ValueError("ORF transcript identity mismatch")
            if match == "other_version":  # The source pins another version of this transcript.
                raise ValueError("ORF transcript version mismatch")
        if _contig(transcript.contig) != _contig(model.contig) or transcript.strand != model.strand:
            raise ValueError("ORF transcript locus mismatch")
        code = genetic_code(model.translation_table)
        layout = _layout or _exon_layout(transcript)
        if transcript.sequence is None:
            raise ValueError("ORF transcript sequence unavailable")
        if len(transcript.sequence) != layout[2]:
            raise ValueError("Transcript sequence length mismatch")
        start = path_offset(transcript, model.coding_blocks, layout=layout)
        if start is None:
            raise ValueError("ORF splice path mismatch")
        length = sum(b - a + 1 for a, b in model.coding_blocks)
        coding = str(transcript.sequence[start : start + length]).upper()
        if any(b not in "ACGT" for b in coding):
            raise ValueError("ORF reference sequence is ambiguous")
        if code.translate(coding) != model.protein_sequence + "*":
            raise ValueError("ORF reference protein/stop mismatch")
        self._set(record=record, transcript=transcript, start=start, end=start + length,
                  coding_sequence=coding, context=resolved_context, kind="extra_orf",
                  start_complete=True, stop_complete=True, annotation_id=None,
                  source_transcript_id=model.transcript_id, equivalent_ids=(), _blocks=model.coding_blocks,
                  _layout=layout)

    @classmethod
    def from_annotation(cls, transcript, *, on=None, context=None):
        """An annotated CDS, on its own transcript or on another one carrying it.

        ``transcript`` is the transcript whose annotated CDS this is. ``on`` is
        a different transcript that contains the identical spliced coding path,
        such as a bicistronic RNA whose "5' UTR" is another gene's CDS; by
        default the CDS is placed on its own transcript. The protein is the
        annotation's, from its protein FASTA. Translation is not re-checked,
        because annotations include selenoproteins and CDSs without start or
        stop codons. ``start_complete`` and ``stop_complete`` report whether
        the annotation marks those codons; the coding sequence includes the
        stop codon when one is marked.
        """
        if not transcript.protein_id:
            raise ValueError(f"Transcript {transcript.id} has no annotated CDS")
        host = transcript if on is None else on
        resolved_context = ReferenceContext.from_genome(host.genome, context=context)
        if on is not None:
            ReferenceContext.from_genome(transcript.genome, context=resolved_context)
        blocks = _merge(transcript.coding_sequence_position_ranges)
        if not blocks:
            raise ValueError(f"Transcript {transcript.id} has no annotated CDS")
        if host.sequence is None:
            raise ValueError("Transcript sequence unavailable")
        layout = _exon_layout(host)
        if len(host.sequence) != layout[2]:
            raise ValueError("Transcript sequence length mismatch")
        start = path_offset(host, blocks, layout=layout) if host.strand == transcript.strand and \
            _contig(host.contig) == _contig(transcript.contig) else None
        if start is None:
            raise ValueError(f"Annotated CDS of {transcript.id} is not one contiguous part of {host.id}")
        length = sum(b - a + 1 for a, b in blocks)
        coding = str(host.sequence[start : start + length]).upper()
        start_complete, stop_complete = transcript.contains_start_codon, transcript.contains_stop_codon
        if on is not None:
            own = path_offset(transcript, blocks)
            if own is None or str(transcript.sequence[own : own + length]).upper() != coding:
                raise ValueError(f"Annotated CDS of {transcript.id} differs in sequence on {host.id}")
            # Without a marked start or stop, the CDS is the same ORF on another RNA only if
            # that RNA continues the same way beyond it; otherwise translation reads on.
            if _versioned(on) != _versioned(transcript) and (
                    (not start_complete and str(host.sequence[:start]) != str(transcript.sequence[:own]))
                    or (not stop_complete and str(host.sequence[start + length :])
                        != str(transcript.sequence[own + length :]))):
                raise ValueError(f"Annotated CDS of {transcript.id} lacks a marked start or stop codon "
                                 f"and {host.id} continues differently beyond it")
        version = transcript.protein.protein_version
        orf = object.__new__(cls)
        orf._set(record=None, transcript=host, start=start, end=start + length, coding_sequence=coding,
                 context=resolved_context, kind="annotated", start_complete=start_complete,
                 stop_complete=stop_complete, _layout=layout,
                 annotation_id=transcript.protein_id + ("" if version is None else f".{version}"),
                 source_transcript_id=_versioned(transcript), equivalent_ids=(), _blocks=blocks,
                 _source_gene=(transcript.gene_id, transcript.gene_name),
                 _protein=transcript.protein_sequence)
        return orf

    def _set(self, **values):
        for name, value in values.items():
            object.__setattr__(self, name, value)

    @property
    def cdna_interval(self):
        """``(start, end)``, the slice of the transcript sequence this ORF reads."""
        return self.start, self.end

    @property
    def transcript_sequence(self):
        return str(self.transcript.sequence).upper()

    @property
    def genomic_positions(self):
        """Forward-genomic position of every base of the transcript, in RNA order."""
        positions = self.__dict__.get("_positions")
        if positions is None:
            positions = transcript_positions(self.transcript)
            object.__setattr__(self, "_positions", positions)
        return positions

    @property
    def initiation_context(self):
        """Up to 9 bases before the start codon and 6 after it, on this transcript.

        An extra ORF on its source's own transcript keeps the source's context.
        Elsewhere, and for an annotated CDS, the flanks are read from this
        transcript, whose 5' sequence can differ from the source's. None when
        the start codon is not complete.
        """
        if "_initiation" not in self.__dict__:
            context = None
            if self.record is not None and self.record.initiation_context is not None and self.on_source_transcript:
                context = self.record.initiation_context
            elif self.start_complete:
                low = max(0, self.start - 9)
                window = str(self.transcript.sequence[low : self.start + 9])
                context = annotate_initiation_context(window, self.start - low)
                layout = self.__dict__.get("_layout") or _exon_layout(self.transcript)
                positions = tuple(_layout_position(layout, self.transcript.strand, offset)
                                  for offset in range(low, low + len(context.upstream + context.start_codon
                                                                     + context.downstream)))
                context = replace(context, genomic_positions=positions, start_offset=self.start,
                                  source_transcript_id=_versioned(self.transcript))
            object.__setattr__(self, "_initiation", context)
        return self._initiation

    @property
    def model(self):
        """The source ORF model; None for an annotated CDS."""
        return None if self.record is None else self.record.model

    @property
    def evidence(self):
        """Source evidence; an annotated CDS carries none from orfdata."""
        if self.record is None:
            from .catalog import ORFEvidence

            return ORFEvidence()
        return self.record.evidence

    @property
    def protein_sequence(self):
        return self._protein if self.record is None else self.model.protein_sequence

    @property
    def genetic_code(self):
        """The ORF's declared NCBI code; modified products must keep using it.

        None for an annotated CDS: the annotation does not declare its code.
        """
        return None if self.record is None else genetic_code(self.model.translation_table)

    @property
    def orf_id(self):
        """The source ORF ID, or the versioned protein ID of an annotated CDS."""
        return self.annotation_id if self.record is None else self.model.orf_id

    id = orf_id

    @property
    def orf_type(self):
        """Source ORF category, kept verbatim; None for an annotated CDS."""
        return None if self.record is None else self.model.orf_type

    @property
    def coding_blocks(self):
        """1-based inclusive genomic blocks of the coding path, in ascending order."""
        return self._blocks

    @property
    def on_source_transcript(self):
        """Whether this is the transcript the ORF's source defines it on."""
        return _transcript_match(self.source_transcript_id, self.transcript) in ("id", "versioned")

    @property
    def transcript_id(self):
        return self.transcript.id

    @property
    def gene(self):
        return self.transcript.gene

    @property
    def gene_id(self):
        return self.transcript.gene_id

    @property
    def gene_name(self):
        return self.transcript.gene_name

    @property
    def source_gene_id(self):
        return self._source_gene[0] if self.record is None else self.model.gene_id

    @property
    def source_gene_name(self):
        return self._source_gene[1] if self.record is None else self.model.gene_name

    def to_dict(self):
        transcript = self.transcript
        genome = transcript.genome
        return dict(
            kind=self.kind,
            orf_id=self.orf_id,
            equivalent_ids=self.equivalent_ids,
            source_transcript_id=self.source_transcript_id,
            on_source_transcript=self.on_source_transcript,
            model=None if self.model is None else self.model.to_dict(),
            evidence=self.evidence.to_dict(),
            reference_context=self.context.to_dict(),
            reference_key=self.context.key,
            sequence_fingerprints={
                name: hashlib.sha256(sequence.encode("ascii")).hexdigest()
                for name, sequence in (
                    ("transcript_sha256", self.transcript_sequence),
                    ("coding_sha256", self.coding_sequence),
                    ("protein_sha256", self.protein_sequence),
                )
            },
            resolved_transcript=dict(
                id=transcript.id,
                version=getattr(transcript, "version", None),
                release=getattr(genome, "release", None),
                genome_build=genome.reference_name,
                species=_species(genome.species) if getattr(genome, "species", None) else None,
                biotype=getattr(transcript, "biotype", None),
            ),
            start=self.start,
            end=self.end,
            cdna_interval=self.cdna_interval,
            coding_blocks=self.coding_blocks,
            start_complete=self.start_complete,
            stop_complete=self.stop_complete,
            coordinates="0-based-half-open-cdna;1-based-inclusive-genome",
        )


def _merge(ranges):
    """Sorted, merged 1-based inclusive genomic ranges."""
    merged = []
    for a, b in sorted(ranges):
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return tuple(tuple(x) for x in merged)


def _versioned(feature):
    return str(getattr(feature, "versioned_id", None) or feature.id)


# The original name; Varcode and Isovar import it.
ReferenceORF = TranscriptORF
