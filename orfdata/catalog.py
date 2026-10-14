"""Versioned reference ORFs and source-specific translation evidence.

Genomic coordinates are 1-based inclusive on the explicitly named assembly.
These are reference observations, not evidence for a patient's mutant protein.
"""

import gzip
import hashlib
import json
from collections import Counter
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from importlib.resources import files
from pathlib import Path
from typing import NamedTuple, Optional

from .evidence import EVIDENCE_METHODS as EVIDENCE_METHODS  # Re-exported for existing imports.
from .evidence import evidence_rule
from .evidence_summary import EvidenceSummary, unique_summaries
from .initiation_context import KOZAK_PREFERENCES, InitiationContext
from .liftover import Liftover
from .model import ORFModel, TranscriptORF, _contig, _id_forms, _species, _transcript_match
from .reference import ReferenceContext
from .sequence import genetic_code, reverse_complement_dna


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
        """Whether this evidence satisfies an evidence rule or method name.

        ``method`` is a rule such as ``RiboSeq(min_studies=2) & MS()``, or a
        method name with ``require_*`` options (see :func:`evidence_rule`).
        Aggregate counts never satisfy uniqueness or spectrum-review filters.
        No selection establishes evidence for a modified protein product.
        """
        if method is None:
            raise ValueError("Unknown evidence method: None")
        return evidence_rule(method, require_unique=require_unique, require_reviewed=require_reviewed,
                             require_resolved=require_resolved).matches(self)


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
    members: tuple = ()
    lifted_from: Optional["ORFRecord"] = None
    liftover: Optional[Liftover] = None

    def __post_init__(self):
        object.__setattr__(self, "summaries", unique_summaries(self.summaries))
        object.__setattr__(self, "members", tuple(self.members))
        if (self.lifted_from is None) != (self.liftover is None):
            raise ValueError("A lifted record names both its original record and the chain")
        original = self.lifted_from
        if original is not None and any(getattr(self, name) != getattr(original, name) for name in (
                "orf_id", "protein_sequence", "translation_table", "start_codon", "source",
                "ribosome_studies", "peptides", "summaries")):
            raise ValueError(f"A lifted record keeps its original's ID, protein and evidence: {self.orf_id}")
        if self.members:
            # Each member validated its own summaries against its own identity.
            key = _merge_key(self)
            if any(m.members or _merge_key(m) != key for m in self.members):
                raise ValueError(f"Members of {self.orf_id} must be source records of one ORF")
            union = _union(self.members)
            # Studies may be fewer: one study labelled differently by two sources is listed once.
            studies = union.pop("ribosome_studies")
            if (any(getattr(self, name) != value for name, value in union.items())
                    or not set(self.ribosome_studies) <= set(studies)
                    or bool(self.ribosome_studies) != bool(studies)
                    or self.initiation_context not in (None, *(m.initiation_context for m in self.members))):
                raise ValueError(f"The evidence of {self.orf_id} must be its members' evidence")
        elif self.summaries and self.lifted_from is None:  # A lifted record's original validated them.
            model = self.model
            for summary in self.summaries:
                summary.validate_orf(model)

    @property
    def source_ids(self):
        """``source:id`` for each source record of this ORF, in catalogue order."""
        return tuple(f"{m.source}:{m.orf_id}" for m in self.members or (self,))

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
        if self.members:
            result["members"] = [m.to_dict() for m in self.members]
        else:
            result.pop("members")
        if self.lifted_from is not None:
            result["lifted_from"] = self.lifted_from.to_dict()
        else:
            result.pop("lifted_from")
            result.pop("liftover")
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
            Must equal the record's assembly. No implicit liftover is done; see ``lift_catalog``.

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
        self._aliases = {}
        for r in self.records:
            names = {f"{m.source}:{m.orf_id}" for m in r.members or (r,)} | {m.orf_id for m in r.members}
            # A source accession split by coding path, "IP_1@ENST...", still answers to "IP_1".
            names |= {m.orf_id.split("@")[0] for m in r.members or (r,) if "@" in m.orf_id}
            for alias in names - {r.orf_id}:
                self._aliases.setdefault(alias, []).append(r)
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
        """Return the ORF with this ID, bare or as ``source:id``; unknown IDs raise KeyError.

        In a merged catalogue every source's ID finds its ORF. A bare ID that
        two sources use for different ORFs raises rather than choosing one.
        """
        found = list(dict.fromkeys(
            ([self._by_id[orf_id]] if orf_id in self._by_id else []) + list(self._aliases.get(orf_id, ()))))
        if len(found) > 1:
            raise ValueError(f"{orf_id} names several ORFs: {', '.join(r.orf_id for r in found)}")
        if not found:
            raise KeyError(orf_id)
        return found[0]

    @property
    def source(self):
        """The catalogue's release label, as proteome layers and entries record it."""
        metadata = self._metadata
        if metadata.get("liftover"):
            return f"{metadata['liftover']['source']}@{metadata['genome_build']}"
        if metadata.get("merged_from"):
            return "+".join(m["source"] for m in metadata["merged_from"])
        manifest = metadata.get("catalog_manifest")
        if manifest:
            return f"{manifest['catalog_id']}:{manifest['version']}"
        return f"orfdata-catalog:{metadata.get('dataset_version', 'unversioned')}"

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
        found = []
        for record in self.query(**filters):
            if _transcript_match(record.transcript_id, transcript) is not None:
                found.append(record.bind(transcript, context=resolved))
            elif any(_transcript_match(m.transcript_id, transcript) is not None for m in record.members):
                # Another source's transcript: it must carry the exact path and protein.
                found.append(record.bind(transcript, context=resolved, any_transcript=True))
        return tuple(found)

    def to_dict(self):
        """Portable catalogue with explicit species and source provenance."""
        return dict(
            schema="orfdata.catalog.v1",
            metadata=self.metadata,
            records=[r.to_dict() for r in self.records],
        )

    def save(self, path):
        """Write the catalogue as JSON, gzipped when ``path`` ends in ``.gz``; ``load_catalog(path)`` reads it.

        The same catalogue always writes the same bytes, so a file's checksum identifies it.
        """
        raw = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        with open(path, "wb") as handle:
            if str(path).endswith(".gz"):
                with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
                    compressed.write(raw)
            else:
                handle.write(raw)

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
        ``evidence`` is a method name or a rule such as
        ``RiboSeq(min_studies=2) & MS(unique=True)``; see :mod:`orfdata.evidence`.
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
        rule = evidence_rule(evidence, require_unique=require_unique, require_reviewed=require_reviewed,
                             require_resolved=require_resolved)

        return tuple(
            r
            for r in self.records
            if (gene is None or any(gene in (m.gene_name, m.gene_id) for m in r.members or (r,)))
            and (transcript_id is None or any(m.transcript_id == transcript_id for m in r.members or (r,)))
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
            and (rule is None or rule.matches(r.evidence))
        )


def _merge_key(record):
    """What makes two source records one ORF: path, protein and translation.

    Declaring a start codon the genetic code already starts at changes nothing.
    """
    start = record.start_codon
    try:
        if start is not None and start in genetic_code(record.translation_table).start_codons:
            start = None
    except ValueError:  # An unsupported code keeps its declaration as given.
        pass
    return (_contig(record.contig), record.strand, record.coding_blocks, record.protein_sequence,
            record.translation_table, start, record.genome_build, record.taxon_id)


def _path_key(record):
    return _merge_key(record)[:3] + _merge_key(record)[6:]


def _strain(record):
    context = record.reference_context
    return (None, None) if context is None else (context.strain, context.strain_taxon_id)


def _union(members, study=None):
    """One ORF's evidence from several sources: identical observations listed once.

    ``study(member, label)`` gives a study's identity, such as its PMID, so a
    study two sources label differently is listed once, under its first label.
    """
    studies = {}
    for m in members:
        for label in m.ribosome_studies:
            studies.setdefault(label if study is None else study(m, label), label)
    return dict(ribosome_studies=tuple(studies.values()),
                peptides=tuple(dict.fromkeys(x for m in members for x in m.peptides)),
                summaries=unique_summaries(x for m in members for x in m.summaries))


def _shared_context(primary, members):
    """The first source initiation context that describes the primary record's transcript."""
    return next((m.initiation_context for m in members if m.initiation_context is not None and (
        m.initiation_context.source_transcript_id is None
        or set(_id_forms(m.initiation_context.source_transcript_id)) & set(_id_forms(primary.transcript_id)))),
        None)


def merge_catalogs(catalogs):
    """Merge ORF catalogues for one species, assembly and strain into one catalogue.

    Records from different catalogues with the same spliced coding path,
    protein, genetic code and start codon are one ORF; records within one
    catalogue stay as it gave them. A record joins an ORF found on its own
    transcript before any other. The ORF's record
    keeps every source record in ``members``, in catalogue order, and the
    union of their evidence; each peptide keeps its source, and identical
    observations are listed once. Different starts, stops, proteins or codes
    stay separate ORFs, and ``metadata["conflicts"]`` lists the ones that share
    a coding path. Each ORF keeps the first source's ID; when two sources use
    one ID for different ORFs, both become ``source:id``. Every source's ID
    still finds the ORF through ``get``. Merged catalogues can be merged again.

    Catalogues must share species, assembly and strain, an unknown strain
    included, and must not repeat a source, such as two releases of one
    catalogue. A Ribo-seq study that two catalogues label differently is
    listed, and counted by ``RiboSeq(min_studies=...)``, once when their
    ``ribosome_study_pmids`` give it the same PMID; otherwise it counts twice.
    Queries by transcript or gene match any source's labels.
    """
    catalogs = tuple(catalogs)
    if not catalogs or any(not isinstance(c, ORFCatalog) for c in catalogs):
        raise TypeError("merge_catalogs needs one or more ORFCatalog objects")
    if len(catalogs) == 1:
        return catalogs[0]
    first = catalogs[0]
    for catalog in catalogs[1:]:
        if (catalog.taxon_id, catalog.genome_build) != (first.taxon_id, first.genome_build):
            raise ValueError(f"Cannot merge a {_species(catalog.species)} {catalog.genome_build} catalogue "
                             f"with a {_species(first.species)} {first.genome_build} one")
    # A merged input is merged again from its source records, each source its own input.
    sources = [[m for r in c.records for m in (r.members or (r,))] for c in catalogs]
    # Records join across inputs, never within one; a merged input's sources were separate inputs.
    origin = {id(m): (i, m.source if c.metadata.get("merged_from") else None)
              for i, c in enumerate(catalogs) for r in c.records for m in (r.members or (r,))}
    pmids = {}
    for c in catalogs:
        for label, pmid in (c.metadata.get("ribosome_study_pmids") or {}).items():
            pmids.setdefault(label, pmid)
    sources_pmids = {m.source: c.metadata.get("ribosome_study_pmids") or {}
                     for c in catalogs for r in c.records for m in (r.members or (r,))}

    def study(member, label):
        pmid = sources_pmids.get(member.source, {}).get(label)
        return label if pmid is None else f"PMID:{pmid}"
    strains = {_strain(m) for records in sources for m in records}
    if len(strains) > 1:
        raise ValueError("Cannot merge catalogues for different or unknown strains: "
                         + ", ".join(sorted(str(s[0]) for s in strains)))
    seen = set()
    for records in sources:
        names = {m.source for m in records}
        if names & seen:
            raise ValueError(f"Two catalogues contain {min(names & seen)} records; merge one release of each source")
        seen |= names
    groups, by_key = [], {}
    for records in sources:
        # Records on the transcript of an ORF found so far join it first, then the rest.
        later = []
        for passing in (0, 1):
            for record in (records if passing == 0 else later):
                key, forms = _merge_key(record), set(_id_forms(record.transcript_id))
                eligible = [g for g in by_key.get(key, ()) if origin[id(record)] not in g[0]]
                group = next((g for g in eligible if any(forms & set(_id_forms(m.transcript_id)) for m in g[1])),
                             None)
                if group is None and passing == 0 and eligible:
                    later.append(record)
                    continue
                if group is None:
                    group = eligible[0] if eligible else None
                if group is None:
                    group = (set(), [])
                    by_key.setdefault(key, []).append(group)
                    groups.append(group)
                group[0].add(origin[id(record)])
                group[1].append(record)
    uses = Counter(orf_id for _, members in groups for orf_id in {m.orf_id for m in members})
    records = []
    for _, members in groups:
        primary = members[0]
        orf_id = primary.orf_id if uses[primary.orf_id] == 1 else f"{primary.source}:{primary.orf_id}"
        if len(members) == 1 and orf_id == primary.orf_id:
            records.append(primary)
            continue
        pinned = {None if m.reference_context is None else m.reference_context.key for m in members}
        records.append(replace(
            # The merged ORF is no one source's lift; each member keeps its own lift provenance.
            primary, orf_id=orf_id, members=tuple(members), lifted_from=None, liftover=None,
            initiation_context=_shared_context(primary, members),
            # Sources pinned to different annotations of one assembly leave the ORF unpinned.
            reference_context=primary.reference_context if len(pinned) == 1 else None, **_union(members, study)))
    paths = {}
    for record in records:
        paths.setdefault(_path_key(record), set()).add((_merge_key(record), record.orf_id))
    conflicts = [sorted(i for _, i in found) for found in paths.values() if len({k for k, _ in found}) > 1]
    metadata = dict(
        species=first.species, taxon_id=first.taxon_id, genome_build=first.genome_build,
        coordinates="1-based-inclusive",
        sources=[s for c in catalogs for s in c.metadata.get("sources", ())],
        merged_from=[dict(source=c.source, records=len(c.records),
                          sha256=hashlib.sha256(json.dumps(c.to_dict(), sort_keys=True).encode()).hexdigest(),
                          metadata={k: v for k, v in c.metadata.items() if k != "sources"})
                     for c in catalogs],
        conflicts=sorted(conflicts))
    if pmids:
        metadata["ribosome_study_pmids"] = pmids
    return ORFCatalog(records, metadata)


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
    records = [_record(row, legacy) for row in data["records"]]
    return ORFCatalog(records, metadata)


def _record(row, legacy=False):
    """One ORFRecord from its JSON form, including a merged record's members."""
    row = dict(row)
    if legacy:
        row.update(taxon_id=9606, species="homo_sapiens", translation_table=1, source="gencode_phase1")
    elif not all(k in row for k in ("taxon_id", "species", "translation_table", "genome_build", "source")):
        raise ValueError("Generic ORF records require explicit species, assembly, genetic code and source")
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
    row["members"] = tuple(_record(m) for m in row.get("members", ()))
    if row.get("lifted_from") is not None:
        row["lifted_from"] = _record(row["lifted_from"])
        row["liftover"] = Liftover(**row["liftover"])
    return ORFRecord(**row)


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
    catalogues cover several assemblies of the species.

    Each assembly of the species contributes one catalogue: a bundled one if
    there is one, so the default never changes when another is downloaded,
    otherwise the one catalogue in ``registry`` (by default the packaged
    :class:`CatalogRegistry`), at its default release unless only another
    release covers that assembly. One on another assembly is lifted to
    ``genome_build`` through the pinned chain files, in as many hops as
    needed, and several are merged. Catalogues and chain files must already
    be installed (``orfdata catalogs download``/``build``, ``orfdata chains
    download``); nothing is downloaded here, and a missing one raises with
    the command that installs it. Which catalogues and chains are used
    depends only on the registry and the pins, never on what is installed. A
    species without a catalogue raises, naming those that exist. Absence
    from a catalogue is never evidence of absence.

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
    from .liftover import chain_downloaded, chain_routes, lift_to

    species = _species("human" if species is None else species)
    choices = _catalog_choices(registry)
    found = [c for c in choices if c.species == species]
    builds = sorted({c.genome_build for c in found})
    known = "; ".join(f"{c.species} {c.genome_build} ({c.label})" for c in choices) or "none"
    if not found:
        raise ValueError(f"No ORF catalogue covers {species}{'' if genome_build is None else ' ' + genome_build}. "
                         f"Catalogues exist for: {known}. Load a local file with load_catalog(path) or add one "
                         "to a CatalogRegistry.")
    if genome_build is None:
        if len(builds) > 1:
            raise ValueError(f"ORF catalogues cover several {species} assemblies ({', '.join(builds)}); "
                             "pass genome_build")
        genome_build = builds[0]
    # One catalogue per assembly: on this one as is, on the others lifted here through pinned chains.
    sources = [(c, () if c.genome_build == genome_build else chain_routes(species, c.genome_build).get(genome_build))
               for c in _sources(found, genome_build)]
    sources = [(choice, route) for choice, route in sources if route is not None]
    if not sources:
        raise ValueError(f"No ORF catalogue covers {species} {genome_build}, and no pinned chain files lift one "
                         f"there. Catalogues exist for: {known}. Lift one with your own chain file: 'orfdata lift "
                         f"--species {species} --genome-build {builds[0]} --to {genome_build} --chain FILE --output "
                         f"FILE', or lift_catalog(load_catalog(species={species!r}, genome_build={builds[0]!r}), "
                         f"ChainFile(path, source_build={builds[0]!r}, target_build={genome_build!r})).")
    if expected_sha256 is not None and (len(sources) > 1 or sources[0][1]):
        raise ValueError(f"expected_sha256 checks one catalogue file, but {species} {genome_build} ORFs are "
                         "lifted or merged from several")
    chains = list(dict.fromkeys(c for _, route in sources for c in route))
    try:
        missing, step = [c for c in chains if not chain_downloaded(c)], ""
    except ImportError:  # Without DataCache no chain can have been downloaded.
        missing, step = chains, "install orfdata[data], then "
    if missing:
        accept = " --accept-terms" if any(c.restricted for c in missing) else ""
        raise ValueError(f"{species} ORF catalogues reach {genome_build} by liftover, which needs the chain files "
                         f"{', '.join(c.name for c in missing)}; {step}download them once with 'orfdata chains "
                         f"download --species {species}{accept}' or download_chains(species={species!r}"
                         f"{', accept_terms=True' if accept else ''})")
    catalogs = [lift_to(_open_choice(choice, expected_sha256), genome_build) if route
                else _open_choice(choice, expected_sha256) for choice, route in sources]
    return catalogs[0] if len(catalogs) == 1 else merge_catalogs(catalogs)


def _sources(found, genome_build):
    """The catalogues load_catalog combines for ``genome_build``, the requested assembly's first.

    One release per catalogue, so two releases of one source never count
    twice: its release on ``genome_build`` if it has one, else its default.
    Then one catalogue per assembly: a bundled one if there is one, else the
    one registered catalogue; two are never chosen between.
    """
    chosen = [c for c in found if c.registry is None]
    for name in dict.fromkeys(c.name for c in found if c.registry is not None):
        releases = [c for c in found if c.name == name and c.registry is not None]
        here = [c for c in releases if c.genome_build == genome_build]
        default = releases[0].registry.info(name).version
        chosen.extend(_releases(here) if here else [c for c in releases if c.version == default] or releases[:1])
    result = []
    for build in sorted({c.genome_build for c in chosen}, key=lambda b: (b != genome_build, b)):
        here = [c for c in chosen if c.genome_build == build]
        here = [c for c in here if c.registry is None][:1] or here
        if len(here) > 1:
            raise ValueError(f"Several ORF catalogues cover {here[0].species} {build} "
                             f"({', '.join(c.label for c in here)}); open one with "
                             "CatalogRegistry().open(name, version)")
        result.append(here[0])
    return result


def _open_choice(choice, expected_sha256):
    if choice.registry is None:
        return _load_catalog(expected_sha256=expected_sha256, resource=choice.resource)
    return _open_registered(choice, expected_sha256)


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


def _installed(choice):
    """How a registered release is installed: "downloaded", "built", "damaged" or None. Offline, unhashed."""
    registry, name, version = choice.registry, choice.name, choice.version
    try:
        damaged = False
        if registry.published(name, version):
            status = registry.inspect(name, version, verify_files=False).status
            if status == "available":
                return "downloaded"
            damaged = status != "missing"
        status = registry.inspect_build(name, version, verify_files=False).status
    except ImportError:  # Without DataCache nothing can have been downloaded or built.
        return None
    if status == "available":
        return "built"
    return "damaged" if damaged or status != "missing" else None


def _open_registered(choice, expected_sha256):
    """A published download if installed, else a local build from the pinned sources."""
    registry, name, version = choice.registry, choice.name, choice.version
    installed = _installed(choice)
    if installed in (None, "damaged"):
        verb = "download" if registry.published(name, version) else "build"
        when = " again" if installed else (" first" if verb == "download" else " from its pinned sources")
        how = (f"{verb} it{when} with CatalogRegistry().{verb}({name!r}, {version!r}"
               f"{', force=True' if installed else ''}) or 'orfdata catalogs {verb} {name} --version {version}"
               f"{' --force' if installed else ''}'")
        try:
            from .datasets import _datacache

            _datacache()
        except ImportError:
            how = "install orfdata[data], then " + how
        state = "damaged" if installed else "not installed"
        raise ValueError(f"The {name} {version} ORF catalogue is {state}; {how}")
    if installed == "built":
        inspected = registry.inspect_build(name, version)
        if inspected.status != "available":
            return registry.open_build(name, version)  # Raises with what is wrong.
        if expected_sha256 is not None and hashlib.sha256(
                Path(inspected.files["catalog.json.gz"].path).read_bytes()).hexdigest() != expected_sha256:
            raise ValueError("ORF catalogue checksum mismatch")
        return registry.open_build(name, version)
    if expected_sha256 is not None:
        pinned = registry.download_metadata(name, version)["assets"]["catalog.json.gz"]["sha256"]
        if pinned != expected_sha256:
            raise ValueError("ORF catalogue checksum mismatch")
    return registry.open(name, version)


@dataclass(frozen=True)
class CatalogStatus:
    """One ORF catalogue load_catalog can open, and whether it is ready offline.

    ``selected`` marks the one ``load_catalog(species=..., genome_build=...)``
    uses for its species and assembly. ``installed`` is ``"bundled"``,
    ``"downloaded"``, ``"built"``, ``"damaged"`` or None; ``acquire`` names the explicit step
    that installs a missing one, ``"download"`` or ``"build"``. ``lifts_to``
    lists the other assemblies the pinned chains reach, each with its chains
    and whether they are all downloaded.
    """

    species: str
    genome_build: str
    name: str
    version: Optional[str]
    title: str
    records: int
    selected: bool
    installed: Optional[str]
    acquire: Optional[str]
    lifts_to: tuple = ()


def _lifts(species, genome_build):
    from .liftover import chain_downloaded, chain_routes

    def downloaded(chain):
        try:
            return chain_downloaded(chain)
        except ImportError:  # Without DataCache nothing can be downloaded.
            return False

    return tuple(dict(genome_build=build, chains=[c.name for c in route], downloaded=all(map(downloaded, route)))
                 for build, route in sorted(chain_routes(species, genome_build).items()))


def catalog_status(*, registry=None):
    """Every ORF catalogue for every species, without loading or downloading any.

    Lists the bundled catalogue and each registry release with its species,
    assembly, declared record count and installation, so a missing catalogue
    can be told apart from an empty query.
    """
    choices = _catalog_choices(registry)
    selected = set()
    for key in dict.fromkeys((c.species, c.genome_build) for c in choices):
        found = [c for c in choices if (c.species, c.genome_build) == key]
        bundled = [c for c in found if c.registry is None]
        chosen = bundled or _releases(found)
        if len(chosen) == 1:  # Two candidate catalogues are never chosen between.
            selected.add(chosen[0])
    rows = []
    for c in choices:
        if c.registry is None:
            manifest = json.loads(files("orfdata").joinpath(f"data/{c.resource}/manifest.json").read_text())
            rows.append(CatalogStatus(c.species, c.genome_build, c.name, manifest["dataset_version"],
                                      "Bundled ORF catalogue", manifest["statistics"]["records"], c in selected,
                                      "bundled", None, _lifts(c.species, c.genome_build)))
            continue
        info = c.registry.info(c.name, c.version)
        rows.append(CatalogStatus(c.species, c.genome_build, c.name, c.version, info.title, info.coverage["records"],
                                  c in selected, _installed(c),
                                  "download" if c.registry.published(c.name, c.version) else "build",
                                  _lifts(c.species, c.genome_build)))
    return tuple(rows)
