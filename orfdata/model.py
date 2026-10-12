"""Reference ORF identity and its binding to an annotated transcript."""

from dataclasses import asdict, dataclass

from .sequence import CODONS


@dataclass(frozen=True)
class ORFModel:
    """A source-defined ORF; genomic blocks include its terminal stop codon.

    Coordinates are 1-based inclusive; blocks are in ascending genomic order.
    ``orf_type`` describes this reading frame, never the parent RNA's biotype.
    ``translation_table`` is an NCBI genetic-code ID, not inferred from species
    or a contig name. The initial sequence binder supports table 1 explicitly.
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

    def __post_init__(self):
        if type(self.taxon_id) is not int or self.taxon_id < 1 or not self.species:
            raise ValueError("An explicit species and positive NCBI taxon_id are required")
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
        return asdict(self)


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
    value = str(value).removeprefix("chr")
    return "MT" if value == "M" else value


def _species(value):
    return str(getattr(value, "latin_name", value)).lower().replace(" ", "_")


@dataclass(frozen=True, init=False)
class ReferenceORF:
    """An ORF record validated against a specific PyEnsembl transcript.

    The original transcript remains accessible as ``transcript`` with its
    canonical protein/biotype unchanged. Multiple ReferenceORFs can bind the
    same transcript, including one annotated as noncoding. No downloads occur.
    """

    record: object
    transcript: object
    cdna_interval: tuple
    coding_sequence: str
    transcript_sequence: str
    genomic_positions: tuple

    def __init__(self, record, transcript):
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
        if str(transcript.id).split(".")[0] != model.transcript_id.split(".")[0]:
            raise ValueError("ORF transcript identity mismatch")
        # If the source explicitly pins a transcript version, respect it.
        if "." in model.transcript_id:
            resolved = getattr(transcript, "versioned_id", transcript.id)
            if str(resolved) != model.transcript_id:
                raise ValueError("ORF transcript version mismatch")
        if _contig(transcript.contig) != _contig(model.contig) or transcript.strand != model.strand:
            raise ValueError("ORF transcript locus mismatch")
        if model.translation_table != 1:
            raise ValueError(
                "Unsupported ORF translation table %s; no implicit standard-code fallback"
                % model.translation_table
            )
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
        translated = "".join(CODONS[coding[i : i + 3]] for i in range(0, len(coding), 3))
        if coding[:3] in ("ATG", "CTG", "TTG"):
            translated = "M" + translated[1:]
        if translated != model.protein_sequence + "*":
            raise ValueError("ORF reference protein/stop mismatch")
        for name, value in dict(
            record=record,
            transcript=transcript,
            cdna_interval=(start, start + len(wanted)),
            coding_sequence=coding,
            transcript_sequence=sequence,
            genomic_positions=positions,
        ).items():
            object.__setattr__(self, name, value)

    @property
    def model(self):
        return self.record.model

    @property
    def evidence(self):
        return self.record.evidence

    @property
    def protein_sequence(self):
        return self.model.protein_sequence

    def to_dict(self):
        transcript = self.transcript
        genome = transcript.genome
        return dict(
            model=self.model.to_dict(),
            evidence=self.evidence.to_dict(),
            resolved_transcript=dict(
                id=transcript.id,
                version=getattr(transcript, "version", None),
                release=getattr(genome, "release", None),
                genome_build=genome.reference_name,
                species=_species(genome.species) if getattr(genome, "species", None) else None,
                biotype=getattr(transcript, "biotype", None),
            ),
            cdna_interval=self.cdna_interval,
            coordinates="0-based-half-open-cdna;1-based-inclusive-genome",
        )
