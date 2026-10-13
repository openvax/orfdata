"""Reference ORF identity and its binding to an annotated transcript."""

import hashlib
from dataclasses import asdict, dataclass
from typing import Optional

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
    """One ORF placement on one PyEnsembl transcript.

    ``kind`` is ``"annotated"`` for the transcript's own annotated CDS, or
    ``"extra_orf"`` for a source catalogue ORF validated against the
    transcript's spliced path, sequence and protein. ``placement`` is
    ``"source"`` on the transcript the source declared, or ``"compatible"``
    on another transcript carrying the identical spliced coding path; a
    compatible placement is not evidence of translation on that transcript.
    The transcript, including its biotype and annotated protein, is unchanged,
    and one transcript can carry several ORFs. No downloads occur.
    """

    record: object
    transcript: object
    cdna_interval: tuple
    coding_sequence: str
    transcript_sequence: str
    genomic_positions: tuple
    context: ReferenceContext
    kind: str = "extra_orf"
    placement: str = "source"
    start_complete: bool = True
    stop_complete: bool = True
    annotation_id: Optional[str] = None

    def __init__(self, record, transcript, *, context=None, placement="source"):
        if placement not in ("source", "compatible"):
            raise ValueError("placement must be 'source' or 'compatible'")
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
        if placement == "source":
            if str(transcript.id).split(".")[0] != model.transcript_id.split(".")[0]:
                raise ValueError("ORF transcript identity mismatch")
            # If the source explicitly pins a transcript version, respect it.
            if "." in model.transcript_id:
                resolved = getattr(transcript, "versioned_id", transcript.id)
                if str(resolved) != model.transcript_id:
                    raise ValueError("ORF transcript version mismatch")
        if _contig(transcript.contig) != _contig(model.contig) or transcript.strand != model.strand:
            raise ValueError("ORF transcript locus mismatch")
        code = genetic_code(model.translation_table)
        positions = transcript_positions(transcript)
        if transcript.sequence is None:
            raise ValueError("ORF transcript sequence unavailable")
        sequence = str(transcript.sequence).upper()
        if len(sequence) != len(positions):
            raise ValueError("Transcript sequence length mismatch")
        wanted = model.genomic_positions
        try:
            start = positions.index(wanted[0])
        except ValueError as error:
            raise ValueError("ORF splice path mismatch") from error
        if positions[start : start + len(wanted)] != wanted:
            raise ValueError("ORF splice path mismatch")
        coding = sequence[start : start + len(wanted)]
        if any(b not in "ACGT" for b in coding):
            raise ValueError("ORF reference sequence is ambiguous")
        if code.translate(coding) != model.protein_sequence + "*":
            raise ValueError("ORF reference protein/stop mismatch")
        self._set(record=record, transcript=transcript, cdna_interval=(start, start + len(wanted)),
                  coding_sequence=coding, transcript_sequence=sequence, genomic_positions=positions,
                  context=resolved_context, kind="extra_orf", placement=placement,
                  start_complete=True, stop_complete=True, annotation_id=None)

    @classmethod
    def from_annotation(cls, transcript, *, context=None):
        """The transcript's own annotated CDS, as the annotation defines it.

        The protein is the annotation's, from its protein FASTA; translation is
        not re-checked, because annotations include selenoproteins and CDSs
        without start or stop codons. ``start_complete`` and ``stop_complete``
        report whether the annotation marks those codons; the coding sequence
        includes the stop codon when one is marked.
        """
        if not transcript.protein_id:
            raise ValueError(f"Transcript {transcript.id} has no annotated CDS")
        resolved_context = ReferenceContext.from_genome(transcript.genome, context=context)
        positions = transcript_positions(transcript)
        if transcript.sequence is None:
            raise ValueError("Transcript sequence unavailable")
        sequence = str(transcript.sequence).upper()
        if len(sequence) != len(positions):
            raise ValueError("Transcript sequence length mismatch")
        coding = sorted((p for a, b in transcript.coding_sequence_position_ranges for p in range(a, b + 1)),
                        reverse=transcript.strand == "-")
        start = positions.index(coding[0]) if coding and coding[0] in positions else None
        if start is None or positions[start : start + len(coding)] != tuple(coding):
            raise ValueError(f"Annotated CDS of {transcript.id} is not one contiguous part of its RNA")
        version = transcript.protein.protein_version
        orf = object.__new__(cls)
        orf._set(record=None, transcript=transcript, cdna_interval=(start, start + len(coding)),
                 coding_sequence=sequence[start : start + len(coding)], transcript_sequence=sequence,
                 genomic_positions=positions, context=resolved_context, kind="annotated",
                 placement="source", start_complete=transcript.contains_start_codon,
                 stop_complete=transcript.contains_stop_codon,
                 annotation_id=transcript.protein_id + ("" if version is None else f".{version}"))
        return orf

    def _set(self, **values):
        for name, value in values.items():
            object.__setattr__(self, name, value)

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
        return self.transcript.protein_sequence if self.record is None else self.model.protein_sequence

    @property
    def genetic_code(self):
        """The ORF's declared NCBI code; modified products must keep using it.

        None for an annotated CDS: the annotation does not declare its code.
        """
        return None if self.record is None else genetic_code(self.model.translation_table)

    # Identity in the bound annotation; the source's own labels stay on ``model``.
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
        if self.record is not None:
            return self.model.coding_blocks
        start, end = self.cdna_interval
        blocks = []
        for position in sorted(self.genomic_positions[start:end]):
            if blocks and position == blocks[-1][1] + 1:
                blocks[-1][1] = position
            else:
                blocks.append([position, position])
        return tuple(tuple(b) for b in blocks)

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

    def to_dict(self):
        transcript = self.transcript
        genome = transcript.genome
        return dict(
            kind=self.kind,
            placement=self.placement,
            orf_id=self.orf_id,
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
            cdna_interval=self.cdna_interval,
            start_complete=self.start_complete,
            stop_complete=self.stop_complete,
            coordinates="0-based-half-open-cdna;1-based-inclusive-genome",
        )


# The original name; Varcode and Isovar import it.
ReferenceORF = TranscriptORF
