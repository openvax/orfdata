"""Source-defined RNAs, RNA features such as IRESs, and processed protein products.

All RNA coordinates are Python slices into the RNA's spliced sequence, and all
protein coordinates are Python slices into the protein sequence. Genomic blocks
are 1-based and inclusive.
"""

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Optional

from .reference import ReferenceContext

IRES = ("internal_ribosome_entry_site", "SO:0000243")
_TOPOLOGIES = ("linear", "circular")


@dataclass(frozen=True)
class SourceRNA:
    """An RNA defined by a source rather than by the loaded annotation.

    Use it for a circular RNA, a novel or retired transcript, or a sequence-only
    reference. ``exons`` are 1-based inclusive genomic blocks in RNA order that
    cover the sequence exactly once; a circular RNA's junction joins its last
    base to its first. ``related_transcript_ids`` names annotated transcripts it
    derives from; no Ensembl transcript identity is invented for it.
    """

    id: str
    sequence: str
    reference: ReferenceContext
    source: str
    topology: str = "linear"
    contig: Optional[str] = None
    strand: Optional[str] = None
    exons: tuple = ()
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    biotype: Optional[str] = None
    related_transcript_ids: tuple = ()

    def __post_init__(self):
        object.__setattr__(self, "sequence", str(self.sequence).upper())
        object.__setattr__(self, "exons", tuple(tuple(e) for e in self.exons))
        object.__setattr__(self, "related_transcript_ids", tuple(self.related_transcript_ids))
        if not self.id or not self.source or not self.sequence or set(self.sequence) - set("ACGTN"):
            raise ValueError("A source RNA needs an ID, a source and an A/C/G/T/N sequence")
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("A source RNA needs a ReferenceContext")
        if self.topology not in _TOPOLOGIES:
            raise ValueError(f"topology must be one of {_TOPOLOGIES}")
        if self.exons:
            if self.strand not in ("+", "-") or not self.contig:
                raise ValueError("Genomic exons need a contig and strand")
            if any(a < 1 or b < a for a, b in self.exons):
                raise ValueError("Exons are 1-based inclusive blocks")
            ascending = sorted(self.exons)
            if list(self.exons) != (ascending[::-1] if self.strand == "-" else ascending) or any(
                    b >= c for (_, b), (c, _) in zip(ascending, ascending[1:])):
                raise ValueError("Exons must be in RNA order, descending on the minus strand, without overlaps")
            if sum(b - a + 1 for a, b in self.exons) != len(self.sequence):
                raise ValueError("Exons must cover the RNA sequence exactly once")

    @property
    def versioned_id(self):
        return self.id

    @property
    def name(self):
        return self.id

    @property
    def transcript_id(self):
        return self.id

    @property
    def genome(self):
        return None

    @property
    def gene(self):
        return None

    @property
    def protein_id(self):
        return None

    @property
    def layout(self):
        """Exons in RNA order, their offsets in the sequence, and its length."""
        offsets, total = [], 0
        for a, b in self.exons:
            offsets.append(total)
            total += b - a + 1
        return self.exons, tuple(offsets), len(self.sequence)

    @classmethod
    def from_transcript(cls, transcript, start, end, *, id, source, topology="linear", context=None,
                        biotype=None):
        """The slice ``transcript.sequence[start:end]`` as a source RNA.

        For example a circular RNA formed by back-splicing exons that are
        consecutive in an annotated transcript. Its genomic exons follow from
        the transcript's, and the transcript is recorded as related.
        """
        from .model import _exon_layout, _layout_blocks, _versioned

        if not 0 <= start < end <= len(transcript.sequence):
            raise ValueError("The slice must lie within the transcript sequence")
        reference = ReferenceContext.from_genome(transcript.genome, context=context)
        blocks = _layout_blocks(_exon_layout(transcript), transcript.strand, range(start, end), ordered=True)
        return cls(id, str(transcript.sequence)[start:end], reference, source, topology, transcript.contig,
                   transcript.strand, blocks, transcript.gene_id, transcript.gene_name or None, biotype,
                   (_versioned(transcript),))


@dataclass(frozen=True)
class RNAFeature:
    """A source-defined element of an RNA, such as an IRES.

    ``start`` and ``end`` slice the RNA's sequence. ``kind`` and ``so_term``
    name the Sequence Ontology type, for example ``IRES``. ``evidence`` says
    what the source establishes; a sequence-family match does not show activity
    in every transcript or condition. ``linked_orf_ids`` lists ORFs the source
    reports this element initiating; none is inferred.
    """

    feature_id: str
    kind: str
    so_term: Optional[str]
    rna: object
    start: int
    end: int
    source: str
    evidence: str = ""
    linked_orf_ids: tuple = ()

    def __post_init__(self):
        if not 0 <= self.start < self.end <= len(self.rna.sequence):
            raise ValueError("A feature must be a non-empty slice of its RNA's sequence")
        object.__setattr__(self, "linked_orf_ids", tuple(self.linked_orf_ids))

    @classmethod
    def find(cls, rna, element_sequence, *, feature_id, kind, so_term, source, evidence="", linked_orf_ids=()):
        """Place an element at its one exact occurrence in the RNA's sequence.

        Raises when the sequence occurs nowhere or more than once, rather than
        guessing a location.
        """
        target = str(element_sequence).upper().replace("U", "T")
        sequence = str(rna.sequence).upper()
        start = sequence.find(target)
        if start < 0:
            raise ValueError(f"{feature_id} does not occur in {rna.id}")
        if sequence.find(target, start + 1) >= 0:
            raise ValueError(f"{feature_id} occurs more than once in {rna.id}")
        return cls(feature_id, kind, so_term, rna, start, start + len(target), source, evidence, linked_orf_ids)

    @property
    def sequence(self):
        return str(self.rna.sequence)[self.start : self.end].upper()

    @property
    def genomic_blocks(self):
        """1-based inclusive genomic blocks of the element, ascending; empty if unmapped."""
        from .model import _layout_blocks, _rna_layout

        layout = _rna_layout(self.rna)
        return _layout_blocks(layout, self.rna.strand, range(self.start, self.end)) if layout[0] else ()

    def overlaps(self, contig, start, end):
        """Whether the element overlaps a 1-based inclusive genomic interval on this contig."""
        from .model import _contig

        return _contig(contig) == _contig(self.rna.contig) and any(
            a <= end and b >= start for a, b in self.genomic_blocks)

    def to_dict(self):
        return dict(feature_id=self.feature_id, kind=self.kind, so_term=self.so_term,
                    rna_id=self.rna.versioned_id if hasattr(self.rna, "versioned_id") else self.rna.id,
                    start=self.start, end=self.end, sequence=self.sequence,
                    genomic_blocks=[list(b) for b in self.genomic_blocks], source=self.source,
                    evidence=self.evidence, linked_orf_ids=list(self.linked_orf_ids))


@dataclass(frozen=True)
class ProteinProduct:
    """A processed product of a precursor protein, as a slice of its sequence.

    ``kind`` is the source's feature type, such as a signal peptide, propeptide,
    chain or peptide. Products can nest or overlap. A product is part of its
    precursor's translation, not a separately initiated ORF.
    """

    product_id: str
    name: str
    kind: str
    start: int
    end: int
    source: str
    evidence: tuple = ()

    def __post_init__(self):
        if type(self.start) is not int or type(self.end) is not int or not 0 <= self.start < self.end:
            raise ValueError("A product is a non-empty slice of its precursor protein")
        object.__setattr__(self, "evidence", tuple(self.evidence))

    def sequence(self, orf):
        """The product's residues in a precursor ``TranscriptORF``."""
        protein = orf.protein_sequence
        if self.end > len(protein):
            raise ValueError(f"{self.product_id} extends beyond the precursor protein")
        return protein[self.start : self.end]

    def coding_slice(self, orf):
        """``(start, end)`` slicing the RNA sequence to this product's codons.

        For an annotated CDS that begins with a partial codon, residues are
        counted from Ensembl's leading ``X``, which has no complete codon.
        """
        self.sequence(orf)
        if getattr(orf, "wraps_junction", False):
            raise ValueError("Products of an ORF crossing a circular junction have no single slice")
        origin = orf._origin
        if origin is None:
            raise ValueError(f"The reading frame of {orf.orf_id} is unknown")
        if origin + 3 * self.start < 0:
            raise ValueError(f"{self.product_id} starts at a partial codon")
        return orf.start + origin + 3 * self.start, orf.start + origin + 3 * self.end

    def genomic_blocks(self, orf):
        """1-based inclusive genomic blocks of this product's codons, ascending."""
        from .model import _layout_blocks, _rna_layout

        start, end = self.coding_slice(orf)
        layout = _rna_layout(orf.transcript)
        return _layout_blocks(layout, orf.transcript.strand, range(start, end)) if layout[0] else ()


@lru_cache(maxsize=1)
def _rfam_ires_data():
    return json.loads(files("orfdata").joinpath("data/ires/rfam-15.1-human-ires.json").read_text())


def rfam_ires(rna):
    """Human cellular IRES sequences from Rfam 15.1 that occur exactly in an RNA.

    Each match is an :class:`RNAFeature` of kind ``internal_ribosome_entry_site``
    (SO:0000243). The bundled seeds come from the 16 Rfam IRES families with a
    human seed found exactly in Ensembl 115 cDNA. Rfam curates these families
    from published reports, several of them disputed; a match locates a
    reported IRES sequence and does not establish activity in this transcript,
    cell type or condition.
    """
    data = _rfam_ires_data()
    release = data["source"]["release"]
    sequence = str(rna.sequence).upper()
    found = []
    for family in data["families"]:
        pmids = ", ".join(r["pmid"] for r in family["references"])
        matches = {}
        for seed in family["seeds"]:
            start = sequence.find(seed["sequence"])
            while start >= 0:
                matches.setdefault((start, start + len(seed["sequence"])), []).extend(seed["seed_names"])
                start = sequence.find(seed["sequence"], start + 1)
        # Seeds of one family can lie inside each other: report each element once, at its widest.
        spans = sorted(matches, key=lambda span: (span[0], -span[1]))
        for start, end in spans:
            if any(s <= start and end <= e and (s, e) != (start, end) for s, e in spans):
                continue
            names = list(dict.fromkeys(n for s, e in spans if start <= s and e <= end for n in matches[s, e]))
            found.append(RNAFeature(
                f"{family['accession']}:{names[0]}", *IRES, rna, start, end,
                f"Rfam {release} {family['accession']} {family['id']}",
                f"Reported IRES: Rfam seed {', '.join(names)}; PMID {pmids}. "
                "Activity in this transcript or condition is not established."))
    return tuple(sorted(found, key=lambda f: (f.start, f.end, f.feature_id)))


_UNIPROT_PRODUCTS = {"Signal": "signal_peptide", "Propeptide": "propeptide", "Chain": "chain",
                     "Peptide": "peptide", "Transit peptide": "transit_peptide"}


def products_from_uniprot(entry, *, source=None):
    """Processed products from a UniProt JSON entry's features.

    Signal, transit and propeptide, chain and peptide features become
    ``ProteinProduct`` slices of the canonical sequence; UniProt's 1-based
    inclusive positions are converted. Features with uncertain boundaries are
    skipped rather than guessed. Evidence codes are kept.
    """
    accession = entry["primaryAccession"]
    source = source or f"UniProt {accession} entry version {entry['entryAudit']['entryVersion']}"
    products = []
    for number, feature in enumerate(entry.get("features", ()), 1):
        kind = _UNIPROT_PRODUCTS.get(feature.get("type"))
        location = feature.get("location", {})
        begin, end = location.get("start", {}), location.get("end", {})
        if kind is None or begin.get("modifier") != "EXACT" or end.get("modifier") != "EXACT":
            continue
        products.append(ProteinProduct(
            feature.get("featureId") or f"{accession}:{kind}:{number}", feature.get("description") or kind,
            kind, begin["value"] - 1, end["value"], source,
            tuple(e.get("evidenceCode") for e in feature.get("evidences", ()))))
    return tuple(products)
