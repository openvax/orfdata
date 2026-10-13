"""Explicit taxonomy and annotation identities, independent of local names."""

import hashlib
import json
import re
from dataclasses import asdict, dataclass
from typing import Optional

# NCBI Taxonomy IDs; source names are retained on ReferenceContext. This mapping
# is deliberately versioned: changing aliases must never silently change keys.
TAXONOMY_VERSION = "ncbi-2026-10-12"
_SPECIES = {
    9606: ("homo_sapiens", "human"),
    10090: ("mus_musculus", "mouse"),
    9615: ("canis_lupus_familiaris", "canis_familiaris", "dog"),
    4932: ("saccharomyces_cerevisiae", "budding_yeast", "yeast"),
    9685: ("felis_catus", "cat"),
    9913: ("bos_taurus", "cattle", "cow"),
    10116: ("rattus_norvegicus", "rat"),
}
_ALIASES = {name: taxon for taxon, names in _SPECIES.items() for name in names}
_STRAINS = {559292: (4932, "saccharomyces_cerevisiae_s288c")}


def species_name(value):
    """Normalize a display name or PyEnsembl Species without guessing a taxon."""
    return str(getattr(value, "latin_name", value)).strip().lower().replace(" ", "_")


def species_taxon(value):
    """Return a known species taxon, or None for an unregistered name.

    User-defined taxa remain supported through an explicit ReferenceContext;
    an unknown name is never assigned the human taxon or a similar species.
    """
    return _ALIASES.get(species_name(value))


def canonical_species(value):
    """Resolve versioned aliases while retaining unknown scientific names."""
    name = species_name(value)
    taxon = species_taxon(name)
    return _SPECIES[taxon][0] if taxon is not None else name


def _positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive NCBI taxonomy ID")


@dataclass(frozen=True)
class ReferenceContext:
    """One species/strain, assembly and annotation, with a stable scoped key.

    ``annotation_version=None`` means unknown, never the current release.
    Sequence-only or incomplete source references may also leave the assembly
    or annotation provider None, with a pinned sequence/source SHA-256. They
    remain queryable but cannot be used to bind a genomic transcript.
    An assembly accession must include its version. No assembly-name aliases,
    coordinate conversion, orthology or evidence transfer occur here.
    """

    taxon_id: int
    species: str
    genome_build: Optional[str]
    annotation_name: Optional[str]
    annotation_version: Optional[str] = None
    assembly_accession: Optional[str] = None
    strain: Optional[str] = None
    strain_taxon_id: Optional[int] = None
    sequence_sha256: Optional[str] = None
    taxonomy_version: str = TAXONOMY_VERSION
    source_sha256: Optional[str] = None

    def __post_init__(self):
        _positive(self.taxon_id, "taxon_id")
        if not isinstance(self.species, str) or not self.species.strip():
            raise ValueError("Species is required")
        for name in ("genome_build", "annotation_name"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{name} must be a nonempty source label or None")
        if (self.genome_build is None or self.annotation_name is None) and not (
                self.sequence_sha256 or self.source_sha256):
            raise ValueError("An incomplete source reference requires a pinned sequence/source SHA-256")
        if self.taxonomy_version != TAXONOMY_VERSION:
            raise ValueError("Unsupported taxonomy alias version")
        known = species_taxon(self.species)
        if known is not None and known != self.taxon_id:
            raise ValueError("Species name and taxon_id disagree")
        if self.taxon_id in _SPECIES and known != self.taxon_id:
            raise ValueError("Unrecognized name for a registered species taxon")
        if self.taxon_id in _STRAINS:
            raise ValueError("Use species taxon_id and a separate strain_taxon_id")
        if self.strain_taxon_id is not None:
            _positive(self.strain_taxon_id, "strain_taxon_id")
            if self.strain_taxon_id in _SPECIES:
                raise ValueError("Strain taxon cannot be a registered species taxon")
            if self.strain_taxon_id in _STRAINS:
                parent, _ = _STRAINS[self.strain_taxon_id]
                if parent != self.taxon_id:
                    raise ValueError("Strain and species taxonomy disagree")
        if self.annotation_version is not None:
            if isinstance(self.annotation_version, bool) or not str(self.annotation_version):
                raise ValueError("Annotation version must be explicit or None")
            object.__setattr__(self, "annotation_version", str(self.annotation_version))
        if self.assembly_accession is not None and not re.fullmatch(
            r"GC[AF]_[0-9]+\.[0-9]+", self.assembly_accession
        ):
            raise ValueError("Assembly accession must be a versioned GCA/GCF accession")
        if self.sequence_sha256 is not None and not re.fullmatch(
            r"[0-9a-f]{64}", self.sequence_sha256
        ):
            raise ValueError("Sequence-set fingerprint must be a lowercase SHA-256")
        if self.source_sha256 is not None and not re.fullmatch(r"[0-9a-f]{64}", self.source_sha256):
            raise ValueError("Source fingerprint must be a lowercase SHA-256")
        if self.strain is not None and (not isinstance(self.strain, str) or not self.strain):
            raise ValueError("Strain must be a nonempty source label or None")

    @property
    def key(self):
        """Content identity, including strain and reference pins, excluding aliases."""
        values = self.to_dict()
        values.pop("species")
        # Preserve every already-published genomic reference key from 0.2-0.4.
        if values["source_sha256"] is None:
            values.pop("source_sha256")
        payload = json.dumps(values, sort_keys=True, separators=(",", ":")).encode()
        return "ref:" + hashlib.sha256(payload).hexdigest()

    def to_dict(self):
        return asdict(self)

    @property
    def is_genomic(self):
        """Whether an assembly and annotation provider are explicitly declared."""
        return self.genome_build is not None and self.annotation_name is not None

    @classmethod
    def from_genome(cls, genome, *, context=None):
        """Resolve PyEnsembl metadata offline, or validate an explicit custom context.

        An explicit context supplies metadata absent from a custom Genome. Any
        metadata the Genome does contain must still agree. No reference is
        inferred merely from a colliding assembly or transcript identifier.
        """
        species = getattr(genome, "species", None)
        taxon = getattr(genome, "taxon_id", None)
        named_taxon = species_taxon(species) if species is not None else None
        if taxon is not None and named_taxon is not None and taxon != named_taxon:
            raise ValueError("Genome species and taxonomy disagree")
        taxon = taxon if taxon is not None else named_taxon
        if context is None:
            if (species is None or taxon is None or not getattr(genome, "annotation_name", None)
                    or not getattr(genome, "reference_name", None)):
                raise ValueError("Genome species/taxonomy/assembly unavailable; provide an explicit ReferenceContext")
            return cls(
                taxon_id=taxon,
                species=species_name(species),
                genome_build=genome.reference_name,
                annotation_name=genome.annotation_name,
                annotation_version=getattr(genome, "annotation_version", None),
                strain=getattr(genome, "strain", None),
                strain_taxon_id=getattr(genome, "strain_taxon_id", None),
                assembly_accession=getattr(genome, "assembly_accession", None),
                sequence_sha256=getattr(genome, "sequence_sha256", None),
                source_sha256=getattr(genome, "source_sha256", None),
            )
        if not isinstance(context, cls):
            raise TypeError("context must be a ReferenceContext")
        context.validate_genome(genome)
        return context

    def validate_genome(self, genome, *, include_annotation=True):
        """Check available Genome metadata against this explicit declaration.

        Set ``include_annotation=False`` only for genomic-coordinate objects
        such as variants: an allele does not belong to one gene annotation
        release. Taxonomy, strain and assembly checks still apply. This method
        does not claim that missing metadata were observed on the Genome.
        """
        if self.genome_build is None or (include_annotation and self.annotation_name is None):
            raise ValueError("Incomplete source reference cannot bind a Genome; use an explicit reference mapping")
        species = getattr(genome, "species", None)
        taxon = getattr(genome, "taxon_id", None)
        named_taxon = species_taxon(species) if species is not None else None
        if taxon is not None and named_taxon is not None and taxon != named_taxon:
            raise ValueError("Genome species and taxonomy disagree")
        taxon = taxon if taxon is not None else named_taxon
        if taxon is not None and taxon != self.taxon_id:
            raise ValueError("Genome taxon mismatch")
        if species is not None and canonical_species(species) != canonical_species(self.species):
            raise ValueError("Genome species mismatch")
        fields = [
            ("reference_name", self.genome_build),
            ("assembly_accession", self.assembly_accession),
            ("strain", self.strain),
            ("strain_taxon_id", self.strain_taxon_id),
            ("sequence_sha256", self.sequence_sha256),
            ("source_sha256", self.source_sha256),
        ]
        if include_annotation:
            fields.extend((("annotation_name", self.annotation_name),
                           ("annotation_version", self.annotation_version)))
        for attribute, expected in fields:
            actual = getattr(genome, attribute, None)
            if actual is not None and (expected is None or str(actual) != str(expected)):
                raise ValueError(f"Genome {attribute} mismatch")


@dataclass(frozen=True)
class ScopedTarget:
    """A source identifier whose reference, kind and dataset cannot be discarded.

    ``region`` names a source-defined quantification window (for example a UTR
    or intron). It does not assert resolved coordinates or a complete ORF.
    """

    reference: ReferenceContext
    kind: str
    identifier: str
    namespace: str

    def __post_init__(self):
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Target requires a ReferenceContext")
        if self.kind not in ("gene", "transcript", "orf", "protein", "protein_group", "peptide", "region"):
            raise ValueError("Unknown target kind")
        if not all(isinstance(x, str) and x.strip() for x in (self.identifier, self.namespace)):
            raise ValueError("Target identifier and source namespace are required")

    @property
    def key(self):
        return (self.reference.key, self.namespace, self.kind, self.identifier)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_orf(cls, reference_orf, *, namespace):
        """Identify a bound reference ORF in an explicit catalogue/release namespace."""
        from .model import ReferenceORF

        if not isinstance(reference_orf, ReferenceORF):
            raise TypeError("Expected a bound ReferenceORF, not a modified variant/RNA product")
        return cls(reference_orf.context, "orf", reference_orf.model.orf_id, namespace)

    @classmethod
    def from_pyensembl(cls, feature, *, kind, namespace, context=None, versioned=False):
        """Scope a PyEnsembl gene/transcript without implicit release conversion.

        Choose versioned=True only for sources that use versioned IDs. A missing
        versioned identifier is an error; versions are never silently stripped.
        The exact namespace should match the quantitative source's mapping.
        """
        if kind not in ("gene", "transcript"):
            raise ValueError("PyEnsembl features must be gene or transcript targets")
        identifier = getattr(feature, "versioned_id" if versioned else "id", None)
        if not identifier:
            raise ValueError("Feature lacks the requested identifier")
        return cls(ReferenceContext.from_genome(feature.genome, context=context), kind,
                   identifier, namespace)
