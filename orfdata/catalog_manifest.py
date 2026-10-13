"""Pinned provenance and coverage for an optional ORF catalogue."""

import re
from collections.abc import Mapping
from dataclasses import dataclass

from .measurements import SourceAsset, _identifier, _metadata, _Record
from .reference import ReferenceContext


@dataclass(frozen=True)
class CatalogManifest(_Record):
    """One species/assembly catalogue release, independent of abundance data."""

    catalog_id: str
    version: str
    title: str
    references: tuple
    sources: tuple
    citations: tuple
    coverage: Mapping
    processing: Mapping
    schema: str = "orfdata.catalog-manifest.v1"

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["references"] = tuple(ReferenceContext(**c) for c in values["references"])
        values["sources"] = tuple(SourceAsset.from_dict(s) for s in values["sources"])
        return cls(**values)

    def __post_init__(self):
        for name in ("catalog_id", "version", "title"):
            _identifier(getattr(self, name), name)
        for name in ("catalog_id", "version"):
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", getattr(self, name)):
                raise ValueError(f"{name} must be a portable catalogue/path identifier")
        if self.schema != "orfdata.catalog-manifest.v1":
            raise ValueError("Unsupported catalogue manifest schema")
        for name in ("references", "sources", "citations"):
            values = getattr(self, name)
            if isinstance(values, str):
                raise TypeError(f"Catalogue {name} must be a sequence of entries")
            values = tuple(values)
            if not values:
                raise ValueError(f"Catalogue {name} cannot be empty")
            object.__setattr__(self, name, values)
        if any(not isinstance(r, ReferenceContext) for r in self.references):
            raise TypeError("Catalogue references must be ReferenceContexts")
        if len({r.key for r in self.references}) != len(self.references):
            raise ValueError("Duplicate catalogue reference contexts")
        if len({r.taxon_id for r in self.references}) != 1:
            raise ValueError("A catalogue must have one species")
        builds = {r.genome_build for r in self.references if r.genome_build is not None}
        if len(builds) != 1 or not any(r.is_genomic for r in self.references):
            raise ValueError("A catalogue must have one assembly and a complete genomic reference")
        if any(not isinstance(s, SourceAsset) for s in self.sources):
            raise TypeError("Catalogue sources must be SourceAssets")
        if len({s.asset_id for s in self.sources}) != len(self.sources):
            raise ValueError("Duplicate catalogue source assets")
        for citation in self.citations:
            _identifier(citation, "citation")
        _metadata(self, "coverage", "processing")
        count = self.coverage.get("records")
        if type(count) is not int or count < 0:
            raise ValueError("Catalogue coverage must declare a nonnegative records count")
        for key in ("distinct_protein_sequences", "evidence_summaries"):
            if key in self.coverage and (type(self.coverage[key]) is not int or self.coverage[key] < 0):
                raise ValueError(f"Catalogue {key} must be a nonnegative integer")
        if not self.processing:
            raise ValueError("Catalogue processing provenance is required")

    @property
    def taxon_id(self):
        return self.references[0].taxon_id

    @property
    def genome_build(self):
        return next(r.genome_build for r in self.references if r.is_genomic)

    def validate_catalog(self, catalog):
        """Check content identity, accepted records and declared reference scope."""
        if catalog.metadata.get("catalog_manifest") != self.to_dict():
            raise ValueError("Catalogue manifest disagrees with the registered release")
        if catalog.taxon_id != self.taxon_id or catalog.genome_build != self.genome_build:
            raise ValueError("Catalogue species or assembly disagrees with manifest")
        if len(catalog.records) != self.coverage["records"]:
            raise ValueError("Catalogue record count disagrees with manifest")
        references = {r.key for r in self.references}
        for record in catalog.records:
            reference = record.reference_context
            if reference is None or not reference.is_genomic or reference.key not in references:
                raise ValueError("Managed ORF records require a declared complete genomic reference")
        for summary in catalog.evidence_summaries:
            if summary.target.reference.key not in references:
                raise ValueError("Evidence summary reference is absent from catalogue manifest")
        counts = dict(distinct_protein_sequences=len({r.protein_sequence for r in catalog.records}),
                      evidence_summaries=len(catalog.evidence_summaries))
        for key, actual in counts.items():
            if key in self.coverage and self.coverage[key] != actual:
                raise ValueError(f"Catalogue {key} count disagrees with manifest")
