"""Versioned reference ORFs and source-specific translation evidence.

Genomic coordinates are 1-based inclusive on the explicitly named assembly.
These are reference observations, not evidence for a patient's mutant protein.
"""

import gzip
import hashlib
import json
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from importlib.resources import files
from pathlib import Path
from typing import NamedTuple, Optional

from .evidence_summary import EvidenceSummary, unique_summaries
from .initiation_context import KOZAK_PREFERENCES, InitiationContext
from .model import ORFModel, TranscriptORF, _contig, _id_forms, _species, _transcript_match
from .reference import ReferenceContext
from .sequence import genetic_code, reverse_complement_dna

EVIDENCE_METHODS = ("ribosome", "shotgun", "mhc", "hla", "ms_unclassified", "ms")


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class PeptideEvidence:
    """A reported MS peptide, with every compatible protein interval.

    ``protein_intervals`` are 1-based inclusive. ``mapping_count`` and
    spectrum counts have the scope of the named source, not all human proteins
    or independent biological specimens. Missing counts remain unknown.
    """

    method: str
    sequence: str
    protein_intervals: tuple
    source: str
    accession: str
    mapping_count: Optional[int] = None
    spectrum_count: Optional[int] = None
    datasets: tuple = ()
    sequence_match: str = "exact"
    attribution: str = "source_orf"
    quality: str = "reported"
    review_ratings: tuple = ()
    spectrum_ids: tuple = ()

    @property
    def unique_in_source(self):
        """Whether the source reports one protein mapping (not global novelty)."""
        return self.mapping_count == 1

    @property
    def is_supporting(self):
        """Reported or convincingly reviewed, excluding low/mixed reviews."""
        return self.quality in ("reported", "reviewed_support")


@dataclass(frozen=True)
class ORFEvidence:
    """Observations attached to one exact ORF, not a binary protein-existence call."""

    ribosome_studies: tuple = ()
    peptides: tuple = ()
    summaries: tuple = ()

    def __post_init__(self):
        object.__setattr__(self, "summaries", unique_summaries(self.summaries))

    def to_dict(self):
        result = asdict(self)
        if not self.summaries:
            result.pop("summaries")
        return result

    def has_evidence(self, method, *, require_unique=False, require_reviewed=False, require_resolved=False):
        """Select reported support, or require named studies/resolved peptides.

        Aggregate counts never satisfy uniqueness or spectrum-review filters.
        No selection establishes evidence for a modified protein product.
        """
        if method not in EVIDENCE_METHODS:
            raise ValueError(f"Unknown evidence method: {method}")
        if method == "ribosome":
            if require_unique or require_reviewed:
                raise ValueError("MS filters do not apply to ribosome evidence")
            return bool(self.ribosome_studies) or (not require_resolved and any(
                s.method == "ribosome" and s.is_supporting for s in self.summaries))
        resolved = any(
            (method == "ms" or p.method == method or (method == "mhc" and p.method == "hla"))
            and p.is_supporting
            and (not require_unique or p.unique_in_source)
            and (not require_reviewed or p.quality == "reviewed_support")
            for p in self.peptides
        )
        return resolved or (not (require_resolved or require_unique or require_reviewed) and any(
            (method == "ms" and s.method != "ribosome"
             or s.method == method or method == "mhc" and s.method == "hla")
            and s.is_supporting for s in self.summaries))


@dataclass(frozen=True)
class ORFRecord:
    """One reported, spliced reference ORF, not a gene-wide evidence label."""

    orf_id: str
    gene_id: str
    gene_name: str
    transcript_id: str
    contig: str
    strand: str
    coding_blocks: tuple
    protein_sequence: str
    biotype: str
    ribosome_studies: tuple
    peptides: tuple
    genome_build: str = "GRCh38"
    initiation_context: Optional[InitiationContext] = None
    taxon_id: Optional[int] = None
    species: Optional[str] = None
    translation_table: int = 1
    source: str = ""
    reference_context: Optional[ReferenceContext] = None
    summaries: tuple = ()
    start_codon: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "summaries", unique_summaries(self.summaries))
        if self.summaries:
            model = self.model
            for summary in self.summaries:
                summary.validate_orf(model)

    @property
    def model(self):
        """Source-defined ORF identity, independent of its supporting evidence."""
        return ORFModel(
            orf_id=self.orf_id,
            transcript_id=self.transcript_id,
            gene_id=self.gene_id,
            gene_name=self.gene_name,
            contig=self.contig,
            strand=self.strand,
            coding_blocks=self.coding_blocks,
            protein_sequence=self.protein_sequence,
            orf_type=self.biotype,
            genome_build=self.genome_build,
            taxon_id=self.taxon_id,
            species=self.species,
            translation_table=self.translation_table,
            source=self.source,
            reference_context=self.reference_context,
            start_codon=self.start_codon,
        )

    @property
    def evidence(self):
        return ORFEvidence(self.ribosome_studies, self.peptides, self.summaries)

    @classmethod
    def from_model(cls, model, evidence=None, *, initiation_context=None):
        """Construct a record from the generic model and independent evidence."""
        evidence = ORFEvidence() if evidence is None else evidence
        fields = model.to_dict()
        if model.reference_context is not None:
            fields["reference_context"] = model.reference_context
        fields["biotype"] = fields.pop("orf_type")
        return cls(
            **fields,
            ribosome_studies=evidence.ribosome_studies,
            peptides=evidence.peptides,
            summaries=evidence.summaries,
            initiation_context=initiation_context,
        )

    def bind(self, transcript, *, context=None, any_transcript=False):
        """Validate and bind this ORF to a PyEnsembl transcript without relabelling it.

        By default the transcript must be the one the source names. With
        ``any_transcript=True``, any transcript carrying the identical spliced
        coding path and sequence is accepted.
        """
        return TranscriptORF(self, transcript, context=context, any_transcript=any_transcript)

    def to_dict(self):
        """Return a serializable copy of this record."""
        result = asdict(self)
        if not self.summaries:
            result.pop("summaries")
        if self.reference_context is None:
            result.pop("reference_context")
        if self.start_codon is None:
            result.pop("start_codon")
        if self.initiation_context is not None:
            result["initiation_context"] = self.initiation_context.to_dict()
        return result

    def map_initiation_position(self, position, *, genome_build):
        """Map a 1-based genomic base into the spliced initiation context.

        Kozak numbering has A of ATG at +1 and no zero. Return None for
        unavailable context or a base outside the retained window. Upstream
        context, start and coding-body roles remain separate; +4 has both
        initiation-context and coding consequences.
        """
        _positive_integer(position, "position")
        if genome_build != self.genome_build:
            raise ValueError(f"Expected genome build {self.genome_build}, got {genome_build}")
        context = self.initiation_context
        if context is None or position not in context.genomic_positions:
            return None
        index = context.genomic_positions.index(position)
        relative = context.relative_positions[index]
        sequence = context.upstream + context.start_codon + context.downstream
        coding = self.map_genomic_position(position, genome_build=genome_build)
        roles = ["initiation_context"]
        if relative in KOZAK_PREFERENCES:
            roles.append("kozak_preference_position")
        if coding is not None:
            roles.append(coding["role"])
        return {
            "orf_id": self.orf_id,
            "genome_build": self.genome_build,
            "contig": self.contig,
            "position": position,
            "context_position": relative,
            "roles": roles,
            "transcript_base": sequence[index],
            "genomic_ref": sequence[index]
            if self.strand == "+"
            else reverse_complement_dna(sequence[index]),
            "coding_annotation": coding,
            "initiation_assessment": context.assessment,
            "mutant_translation_evidence": "not_assessed",
        }

    def annotate_initiation_snv(self, position, ref, alt, *, genome_build):
        """Compare REF-validated reference/alternate initiation contexts.

        Alleles are single A/C/G/T bases on the forward genomic strand.
        REF mismatch and unavailable/outside context raise ValueError. This
        deliberately does not reconstruct indels, splice events or fusions.
        Directional preference changes are sequence priors, not measured
        initiation changes or proof of a translated mutant protein.
        """
        if ref not in ("A", "C", "G", "T") or alt not in ("A", "C", "G", "T") or ref == alt:
            raise ValueError("Provide distinct single-base forward-genomic REF and ALT")
        annotation = self.map_initiation_position(position, genome_build=genome_build)
        if annotation is None:
            raise ValueError("Initiation context is unavailable or does not cover position")
        if ref != annotation["genomic_ref"]:
            raise ValueError(
                "REF mismatch: expected {}, got {}".format(annotation["genomic_ref"], ref)
            )
        context = self.initiation_context
        index = context.genomic_positions.index(position)
        transcript_alt = alt if self.strand == "+" else reverse_complement_dna(alt)
        sequence = context.upstream + context.start_codon + context.downstream
        altered = sequence[:index] + transcript_alt + sequence[index + 1 :]
        start = len(context.upstream)
        mutant = replace(
            context,
            upstream=altered[:start],
            start_codon=altered[start : start + 3],
            downstream=altered[start + 3 :],
        )
        changes = []
        relative = annotation["context_position"]
        if 1 <= relative <= 3:
            changes.append(
                "ATG_lost"
                if context.start_codon == "ATG"
                else "ATG_gained"
                if mutant.start_codon == "ATG"
                else "start_codon_changed"
            )
        if relative in KOZAK_PREFERENCES:
            preferred = KOZAK_PREFERENCES[relative]
            before, after = sequence[index] in preferred, transcript_alt in preferred
            changes.append(
                "kozak_preference_gained"
                if not before and after
                else "kozak_preference_lost"
                if before and not after
                else "kozak_preference_unchanged"
            )
        codons = None
        if relative > 0:
            offset = start + (relative - 1) // 3 * 3
            before, after = sequence[offset : offset + 3], altered[offset : offset + 3]
            if len(before) == 3 and all(b in "ACGT" for b in before + after):
                table = genetic_code(self.translation_table).codons
                codons = {
                    "reference_codon": before,
                    "alternate_codon": after,
                    "reference_codon_translation": table[before],
                    "alternate_codon_translation": table[after],
                    "interpretation": "codon_translation_if_frame_retained; initiation_not_demonstrated",
                }
        return dict(
            **annotation,
            variant={"ref": ref, "alt": alt, "transcript_alt": transcript_alt},
            sequence_changes=changes,
            reference_context=context.to_dict(),
            alternate_context=mutant.to_dict(),
            codon_consequence=codons,
            translation_change="not_measured",
        )

    def has_evidence(self, method, *, require_unique=False, require_reviewed=False, require_resolved=False):
        return self.evidence.has_evidence(
            method, require_unique=require_unique, require_reviewed=require_reviewed,
            require_resolved=require_resolved,
        )

    def _coding_offset(self, position):
        """Map a placed base without constructing peptide evidence exports."""
        offset = 0
        blocks = self.coding_blocks if self.strand == "+" else self.coding_blocks[::-1]
        for start, end in blocks:
            if start <= position <= end:
                return offset + (position - start if self.strand == "+" else end - position)
            offset += end - start + 1
        return None

    def map_genomic_position(self, position, *, genome_build):
        """Map one genomic base into the spliced ORF and reference MS coverage.

        Parameters
        ----------
        position : int
            A 1-based genomic position on this record's contig.
        genome_build : str
            Must equal the record's assembly. No implicit liftover is done.

        Returns
        -------
        dict or None
            Coding nucleotide offset (0-based), residue position (1-based),
            reference residue, start/stop/body role and covering reference
            peptides. An intronic or outside position returns None. This does
            not compute a mutant sequence or demonstrate mutant translation.
        """
        _positive_integer(position, "position")
        if genome_build != self.genome_build:
            raise ValueError(f"Expected genome build {self.genome_build}, got {genome_build}")
        offset = self._coding_offset(position)
        if offset is None:
            return None
        aa = offset // 3 + 1
        stop = aa > len(self.protein_sequence)
        covered = [p for p in self.peptides if any(a <= aa <= b for a, b in p.protein_intervals)]
        return {
            "orf_id": self.orf_id,
            "genome_build": self.genome_build,
            "coding_nucleotide_offset": offset,
            "amino_acid_position": aa,
            "reference_amino_acid": "*" if stop else self.protein_sequence[aa - 1],
            "role": "stop_codon" if stop else "start_codon" if offset < 3 else "body",
            "reference_peptides": [asdict(p) for p in covered],
            "reference_ms_support": any(p.is_supporting for p in covered),
            "mutant_translation_evidence": "not_assessed",
        }


class ORFCatalog:
    """Offline ORF records for one explicitly identified species and assembly."""

    def __init__(self, records, metadata):
        self.records = tuple(records)
        self._metadata = deepcopy(metadata)
        self.genome_build = metadata["genome_build"]
        self.taxon_id = metadata["taxon_id"]
        self.species = metadata["species"]
        self._by_id = {r.orf_id: r for r in self.records}
        if len(self._by_id) != len(self.records):
            raise ValueError("Duplicate ORF identifiers")
        for r in self.records:
            if r.genome_build != self.genome_build or r.strand not in ("+", "-"):
                raise ValueError(f"Invalid assembly or strand for {r.orf_id}")
            r.model  # Validate identity and geometry, including non-upstream ORFs.
            if (r.taxon_id, _species(r.species)) != (self.taxon_id, _species(self.species)):
                raise ValueError("Catalogue and ORF species disagree")
            if (
                not r.coding_blocks
                or not r.protein_sequence
                or any(a < 1 or b < a for a, b in r.coding_blocks)
                or any(b >= c for (_, b), (c, _) in zip(r.coding_blocks, r.coding_blocks[1:]))
                or sum(b - a + 1 for a, b in r.coding_blocks) != 3 * (len(r.protein_sequence) + 1)
            ):
                raise ValueError(f"Invalid spliced coding blocks for {r.orf_id}")
            context = r.initiation_context
            if (
                context is not None
                and context.source_transcript_id is not None
                and not set(_id_forms(context.source_transcript_id)) & set(_id_forms(r.transcript_id))
            ):
                raise ValueError("Initiation-context transcript identity disagrees with ORF")
            if context is not None and context.genomic_positions:
                positions = context.genomic_positions
                if any(
                    (b <= a if r.strand == "+" else b >= a)
                    for a, b in zip(positions, positions[1:])
                ):
                    raise ValueError("Initiation-context placements disagree with strand")
                for relative, position in zip(context.relative_positions, positions):
                    mapped = r._coding_offset(position)
                    if (relative < 0 and mapped is not None) or (
                        relative > 0 and mapped != relative - 1
                    ):
                        raise ValueError("Initiation-context placements disagree with coding path")
            for p in r.peptides:
                if p.method == "hla" and r.taxon_id != 9606:
                    raise ValueError("Use mhc for non-human immunopeptidomics evidence")
                if (
                    p.method not in ("shotgun", "mhc", "hla", "ms_unclassified")
                    or not p.protein_intervals
                ):
                    raise ValueError(f"Invalid peptide evidence for {r.orf_id}")
                if (
                    p.quality
                    not in (
                        "reported",
                        "reviewed_support",
                        "reviewed_low_quality",
                        "reviewed_mixed",
                    )
                    or p.attribution not in ("source_orf", "sequence_compatible")
                    or p.sequence_match not in ("exact", "I_L_equivalent")
                    or not p.source
                    or not p.accession
                ):
                    raise ValueError(f"Invalid peptide provenance for {r.orf_id}")
                for count in (p.mapping_count, p.spectrum_count):
                    if count is not None:
                        _positive_integer(count, "evidence count")
                for curator, rating in p.review_ratings:
                    if not curator or type(rating) is not int or not 1 <= rating <= 5:
                        raise ValueError("Invalid spectrum review rating")
                if p.quality != "reported" and not p.review_ratings:
                    raise ValueError("Reviewed evidence requires review ratings")
                for a, b in p.protein_intervals:
                    if not 1 <= a <= b <= len(r.protein_sequence):
                        raise ValueError(f"Peptide interval outside {r.orf_id}")
                    if r.protein_sequence[a - 1 : b].replace("I", "L") != p.sequence.replace(
                        "I", "L"
                    ):
                        raise ValueError(f"Peptide sequence does not match {r.orf_id}")

        self._evidence_summaries = unique_summaries(s for r in self.records for s in r.summaries)

    @property
    def evidence_summaries(self):
        """Distinct source counts, without multiplying transcript assignments."""
        return self._evidence_summaries

    @property
    def metadata(self):
        """Return a copy of dataset version, sources, scope and provenance."""
        return deepcopy(self._metadata)

    def get(self, orf_id):
        """Return an exact source ORF ID; unknown IDs raise KeyError."""
        return self._by_id[orf_id]

    def for_transcript(self, transcript, *, context=None, **filters):
        """Return validated ORFs on this transcript, with optional evidence filters.

        A source model that disagrees with the selected annotation raises rather
        than silently transferring its evidence to a different protein.
        """
        if transcript.genome.reference_name != self.genome_build:
            raise ValueError("Catalogue transcript assembly mismatch")
        species = getattr(transcript.genome, "species", None)
        if species is not None and _species(species) != _species(self.species):
            raise ValueError("Catalogue transcript species mismatch")
        resolved = ReferenceContext.from_genome(transcript.genome, context=context)
        if resolved.taxon_id != self.taxon_id:
            raise ValueError("Catalogue transcript taxon mismatch")
        return tuple(
            record.bind(transcript, context=resolved)
            for record in self.query(**filters)
            if _transcript_match(record.transcript_id, transcript) is not None
        )

    def to_dict(self):
        """Portable catalogue with explicit species and source provenance."""
        return dict(
            schema="orfdata.catalog.v1",
            metadata=self.metadata,
            records=[r.to_dict() for r in self.records],
        )

    def query(
        self,
        *,
        gene=None,
        transcript_id=None,
        contig=None,
        start=None,
        end=None,
        genome_build=None,
        evidence=None,
        require_unique=False,
        require_reviewed=False,
        require_resolved=False,
        biotype=None,
        include_initiation_context=False,
    ):
        """Filter records; genomic intervals are 1-based inclusive and spliced.

        A locus query requires an explicit assembly. Gene and transcript lookup
        is a discovery aid, not evidence transfer between transcript isoforms.
        Transcript identifiers match exactly as supplied in the source,
        including versions when present. Empty results do not prove no translation.
        ``require_unique`` means unique within the evidence source's mapping
        scope, not unique in every normal human proteome.
        ``require_resolved`` excludes aggregate-only support: MS needs a peptide
        observation and ribosome evidence needs a named study.
        ``include_initiation_context`` also intersects placed spliced start
        flanks, for example an upstream Kozak SNV. The default remains coding
        blocks only. It does not include introns between context bases.
        """
        locus = any(x is not None for x in (contig, start, end))
        if locus:
            if contig is None or start is None or end is None or genome_build is None:
                raise ValueError("Locus queries require contig, start, end and genome_build")
            _positive_integer(start, "start")
            _positive_integer(end, "end")
            if end < start:
                raise ValueError("end must be >= start")
            contig = _contig(contig)
        if genome_build is not None and genome_build != self.genome_build:
            raise ValueError(f"Expected genome build {self.genome_build}, got {genome_build}")
        if evidence is not None and evidence not in EVIDENCE_METHODS:
            raise ValueError(f"Unknown evidence method: {evidence}")
        if require_unique and evidence not in EVIDENCE_METHODS[1:]:
            raise ValueError("require_unique applies to MS evidence only")
        if require_reviewed and evidence not in EVIDENCE_METHODS[1:]:
            raise ValueError("require_reviewed applies to MS evidence only")
        if require_resolved and evidence is None:
            raise ValueError("require_resolved requires an evidence method")

        return tuple(
            r
            for r in self.records
            if (gene is None or gene in (r.gene_name, r.gene_id))
            and (transcript_id is None or r.transcript_id == transcript_id)
            and (
                not locus
                or (
                    _contig(r.contig) == contig
                    and (
                        any(a <= end and b >= start for a, b in r.coding_blocks)
                        or (
                            include_initiation_context
                            and r.initiation_context is not None
                            and any(
                                start <= p <= end for p in r.initiation_context.genomic_positions
                            )
                        )
                    )
                )
            )
            and (biotype is None or r.biotype == biotype)
            and (
                evidence is None
                or r.has_evidence(
                    evidence, require_unique=require_unique, require_reviewed=require_reviewed,
                    require_resolved=require_resolved,
                )
            )
        )


def _load_catalog(path=None, *, expected_sha256=None, resource="orf-evidence"):
    """Load the packaged reference offline, or a compatible JSON/JSON.gz file.

    No genome downloads, spreadsheet dependencies or LLM are required. The
    catalogue stores protein sequences, not genomic REF alleles; callers must
    validate alleles against their own assembly before reconstructing mutants.
    """
    raw = (
        files("orfdata").joinpath("data/" + resource + "/catalog.json.gz").read_bytes()
        if path is None
        else Path(path).read_bytes()
    )
    if path is None:
        manifest = json.loads(
            files("orfdata").joinpath("data/" + resource + "/manifest.json").read_text()
        )
        if hashlib.sha256(raw).hexdigest() != manifest["catalog_sha256"]:
            raise ValueError("Packaged uORF catalogue checksum mismatch")
    if expected_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ValueError("uORF catalogue checksum mismatch")
    if raw.startswith(b"\x1f\x8b"):
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    if data.get("schema") not in ("isovar.uorf_evidence.v1", "orfdata.catalog.v1"):
        raise ValueError("Unsupported uORF evidence schema")
    if data["metadata"].get("coordinates") != "1-based-inclusive":
        raise ValueError("Unsupported uORF coordinate convention")
    legacy = data["schema"] == "isovar.uorf_evidence.v1"
    metadata = data["metadata"]
    if legacy:
        if metadata.get("genome_build") != "GRCh38":
            raise ValueError("Legacy Isovar data are explicitly human GRCh38 only")
        metadata = dict(metadata, taxon_id=9606, species="homo_sapiens")
    elif not metadata.get("sources"):
        raise ValueError("Generic catalogues require source provenance")
    records = []
    for row in data["records"]:
        row = dict(row)
        if legacy:
            row.update(
                taxon_id=9606, species="homo_sapiens", translation_table=1, source="gencode_phase1"
            )
        elif not all(
            k in row for k in ("taxon_id", "species", "translation_table", "genome_build", "source")
        ):
            raise ValueError(
                "Generic ORF records require explicit species, assembly, genetic code and source"
            )
        row["coding_blocks"] = tuple(tuple(x) for x in row["coding_blocks"])
        row["ribosome_studies"] = tuple(row["ribosome_studies"])
        if row.get("reference_context") is not None:
            row["reference_context"] = ReferenceContext(**row["reference_context"])
        if row.get("initiation_context") is not None:
            context = dict(row["initiation_context"])
            assessment = context.pop("assessment", None)
            context["genomic_positions"] = tuple(context.get("genomic_positions", ()))
            row["initiation_context"] = InitiationContext(**context)
            if assessment is not None and assessment != row["initiation_context"].assessment:
                raise ValueError("Initiation-context assessment disagrees with sequence")
        peptides = []
        for p in row["peptides"]:
            p = dict(p)
            p["protein_intervals"] = tuple(tuple(x) for x in p["protein_intervals"])
            p["datasets"] = tuple(p.get("datasets", ()))
            p["review_ratings"] = tuple(tuple(x) for x in p.get("review_ratings", ()))
            p["spectrum_ids"] = tuple(p.get("spectrum_ids", ()))
            peptides.append(PeptideEvidence(**p))
        row["peptides"] = tuple(peptides)
        row["summaries"] = tuple(EvidenceSummary.from_dict(s) for s in row.get("summaries", ()))
        records.append(ORFRecord(**row))
    return ORFCatalog(records, metadata)


# The released Isovar API was explicitly human/GRCh38 and upstream-only. Keep
# that compatibility surface separate from the generic constructors and loader.
@dataclass(frozen=True)
class UORFRecord(ORFRecord):
    taxon_id: int = 9606
    species: str = "homo_sapiens"
    source: str = "gencode_phase1"


class UORFCatalog(ORFCatalog):
    def __init__(self, records, metadata):
        records = tuple(records)
        if any(r.biotype not in ("uORF", "uoORF") for r in records):
            raise ValueError("Invalid upstream ORF/evidence")
        super().__init__(records, {"taxon_id": 9606, "species": "homo_sapiens", **metadata})

    def query(self, *, biotype=None, **kwargs):
        if biotype is not None and biotype not in ("uORF", "uoORF"):
            raise ValueError(f"Unknown upstream ORF biotype: {biotype}")
        return super().query(biotype=biotype, **kwargs)


UORFPeptide = PeptideEvidence


def load_uorf_catalog(path=None, *, expected_sha256=None):
    """Load the released human upstream-ORF API; use load_catalog for all ORFs."""
    catalog = _load_catalog(path, expected_sha256=expected_sha256, resource="uorf-evidence")
    records = tuple(r for r in catalog.records if r.biotype in ("uORF", "uoORF"))
    return UORFCatalog(records, catalog.metadata)


# Catalogues shipped in the package, by species and assembly.
_BUNDLED = {("homo_sapiens", "GRCh38"): "orf-evidence"}


def load_catalog(path=None, *, species=None, genome_build=None, registry=None, expected_sha256=None):
    """Load the ORF catalogue for a species and assembly, or a local catalogue file.

    ``species`` is a name such as ``"mouse"`` or ``"mus_musculus"``; as in
    PyEnsembl, it defaults to human. ``genome_build`` is needed only when
    catalogues cover several assemblies of the species. A bundled catalogue
    is preferred, so the default never changes when another is downloaded.
    Otherwise the one catalogue in ``registry`` (by default the packaged
    :class:`CatalogRegistry`) covering the species is opened offline, at its
    default release unless only another release covers the assembly; it must
    already be downloaded and is never downloaded here. A species without a
    catalogue raises, naming those that exist. Absence from a catalogue is
    never evidence of absence.

    With ``path``, that file is loaded, and ``species`` and ``genome_build``,
    when given, must match it.
    """
    if path is not None:
        if registry is not None:
            raise ValueError("Pass a catalogue path or a registry, not both")
        catalog = _load_catalog(path, expected_sha256=expected_sha256)
        if species is not None and _species(species) != _species(catalog.species):
            raise ValueError(f"{path} is a {_species(catalog.species)} catalogue, not {_species(species)}")
        if genome_build is not None and genome_build != catalog.genome_build:
            raise ValueError(f"{path} is a {catalog.genome_build} catalogue, not {genome_build}")
        return catalog
    species = _species("human" if species is None else species)
    choices = _catalog_choices(registry)
    found = [c for c in choices if c.species == species]
    builds = sorted({c.genome_build for c in found})
    if genome_build is None and len(builds) > 1:
        raise ValueError(f"ORF catalogues cover several {species} assemblies ({', '.join(builds)}); "
                         "pass genome_build")
    found = [c for c in found if genome_build in (None, c.genome_build)]
    if not found:
        scope = species if genome_build is None else f"{species} {genome_build}"
        known = "; ".join(f"{c.species} {c.genome_build} ({c.label})" for c in choices) or "none"
        raise ValueError(f"No ORF catalogue covers {scope}. Catalogues exist for: {known}. "
                         "Load a local file with load_catalog(path) or add one to a CatalogRegistry.")
    bundled = [c for c in found if c.registry is None]
    if bundled:
        return _load_catalog(expected_sha256=expected_sha256, resource=bundled[0].resource)
    found = _releases(found)
    if len(found) > 1:
        raise ValueError(f"Several ORF catalogues cover {species} {found[0].genome_build} "
                         f"({', '.join(c.label for c in found)}); open one with "
                         "CatalogRegistry().open(name, version)")
    return _open_registered(found[0], expected_sha256)


class _Choice(NamedTuple):
    species: str
    genome_build: str
    name: str
    registry: object = None
    resource: Optional[str] = None
    version: Optional[str] = None

    @property
    def label(self):
        return self.name if self.version is None else f"{self.name} {self.version}"


def _catalog_choices(registry):
    """Every catalogue load_catalog can open: the bundled ones, then each registry release."""
    from .catalogs import CatalogRegistry

    registry = CatalogRegistry() if registry is None else registry
    return [_Choice(s, b, "bundled", resource=r) for (s, b), r in _BUNDLED.items()] + [
        _Choice(_species(m.references[0].species), m.genome_build, m.catalog_id, registry, version=m.version)
        for m in registry.list()]


def _releases(found):
    """One release per registered catalogue: its default version if it qualifies, else its only one."""
    chosen = []
    for name in dict.fromkeys(c.name for c in found):
        releases = [c for c in found if c.name == name]
        default = releases[0].registry.info(name).version
        chosen.extend([c for c in releases if c.version == default] or releases)
    return chosen


def _open_registered(choice, expected_sha256):
    registry, name, version = choice.registry, choice.name, choice.version
    if registry.inspect(name, version, verify_files=False).status == "missing":
        raise ValueError(f"The {name} {version} ORF catalogue is not installed; download it first with "
                         f"CatalogRegistry().download({name!r}, {version!r}) or "
                         f"'orfdata catalogs download {name} --version {version}'")
    if expected_sha256 is not None:
        pinned = registry.download_metadata(name, version)["assets"]["catalog.json.gz"]["sha256"]
        if pinned != expected_sha256:
            raise ValueError("ORF catalogue checksum mismatch")
    return registry.open(name, version)
