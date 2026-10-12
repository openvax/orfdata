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
    An assembly accession must include its version. No assembly-name aliases,
    coordinate conversion, orthology or evidence transfer occur here.
    """

    taxon_id: int
    species: str
    genome_build: str
    annotation_name: str
    annotation_version: Optional[str] = None
    assembly_accession: Optional[str] = None
    strain: Optional[str] = None
    strain_taxon_id: Optional[int] = None
    sequence_sha256: Optional[str] = None
    taxonomy_version: str = TAXONOMY_VERSION

    def __post_init__(self):
        _positive(self.taxon_id, "taxon_id")
        if not all(isinstance(x, str) and x.strip() for x in
                   (self.species, self.genome_build, self.annotation_name)):
            raise ValueError("Species, assembly and annotation provider are required")
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
        if self.strain is not None and (not isinstance(self.strain, str) or not self.strain):
            raise ValueError("Strain must be a nonempty source label or None")

    @property
    def key(self):
        """Content identity, including strain and reference pins, excluding aliases."""
        values = self.to_dict()
        values.pop("species")
        payload = json.dumps(values, sort_keys=True, separators=(",", ":")).encode()
        return "ref:" + hashlib.sha256(payload).hexdigest()

    def to_dict(self):
        return asdict(self)

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
            if species is None or taxon is None or not getattr(genome, "annotation_name", None):
                raise ValueError("Genome species/taxonomy unavailable; provide an explicit ReferenceContext")
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
            )
        if not isinstance(context, cls):
            raise TypeError("context must be a ReferenceContext")
        if taxon is not None and taxon != context.taxon_id:
            raise ValueError("Genome taxon mismatch")
        if species is not None and canonical_species(species) != canonical_species(context.species):
            raise ValueError("Genome species mismatch")
        for attribute, expected in (
            ("reference_name", context.genome_build),
            ("annotation_name", context.annotation_name),
            ("annotation_version", context.annotation_version),
            ("assembly_accession", context.assembly_accession),
            ("strain", context.strain),
            ("strain_taxon_id", context.strain_taxon_id),
            ("sequence_sha256", context.sequence_sha256),
        ):
            actual = getattr(genome, attribute, None)
            if actual is not None and (expected is None or str(actual) != str(expected)):
                raise ValueError(f"Genome {attribute} mismatch")
        return context


@dataclass(frozen=True)
class ScopedTarget:
    """A source identifier whose reference, kind and dataset cannot be discarded."""

    reference: ReferenceContext
    kind: str
    identifier: str
    namespace: str

    def __post_init__(self):
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Target requires a ReferenceContext")
        if self.kind not in ("gene", "transcript", "orf", "protein", "protein_group", "peptide"):
            raise ValueError("Unknown target kind")
        if not all(isinstance(x, str) and x.strip() for x in (self.identifier, self.namespace)):
            raise ValueError("Target identifier and source namespace are required")

    @property
    def key(self):
        return (self.reference.key, self.namespace, self.kind, self.identifier)

    def to_dict(self):
        return asdict(self)
