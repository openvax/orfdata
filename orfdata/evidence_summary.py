"""Source-reported support counts, separate from resolved observations."""

import hashlib
import re
from dataclasses import asdict, dataclass, replace
from typing import Optional
from urllib.parse import urlparse

from .reference import ReferenceContext, ScopedTarget


@dataclass(frozen=True)
class EvidenceSummary:
    """One native count of experimental support for an exact source entity.

    ``value`` is a nonnegative integer or unknown (None), never a confidence
    score or abundance. ``definition`` preserves the source's counting rules;
    a sum of unique peptides per study is not a globally unique peptide count.
    A protein-sequence summary can support several transcript assignments,
    which retain the same ``key`` instead of multiplying independent evidence.
    """

    source: str
    release: str
    target: ScopedTarget
    protein_sha256: str
    method: str
    metric: str
    value: Optional[int]
    definition: str
    scope: str
    source_url: str
    source_sha256: str

    def __post_init__(self):
        if not isinstance(self.target, ScopedTarget):
            raise TypeError("Evidence summaries require a ScopedTarget")
        if self.method not in ("ribosome", "shotgun", "mhc", "hla", "ms_unclassified", "ms"):
            raise ValueError("Unknown summary evidence method")
        if self.method == "hla" and self.target.reference.taxon_id != 9606:
            raise ValueError("Use mhc for non-human immunopeptidomics evidence")
        if self.scope not in ("protein_sequence", "orf") or self.target.kind != {
            "protein_sequence": "protein", "orf": "orf"
        }.get(self.scope):
            raise ValueError("Evidence summary scope and target kind disagree")
        if any(not isinstance(x, str) or not x.strip() for x in (self.source, self.release, self.metric, self.definition)):
            raise ValueError("Source, release, native metric and counting definition are required")
        for value in (self.protein_sha256, self.source_sha256):
            if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
                raise ValueError("Evidence summaries require lowercase SHA-256 pins")
        if not isinstance(self.source_url, str) or urlparse(self.source_url).scheme not in ("https", "http") or not urlparse(self.source_url).netloc:
            raise ValueError("Evidence summaries require a source URL")
        if self.value is not None and (type(self.value) is not int or self.value < 0):
            raise ValueError("A support count must be a nonnegative integer or None")

    @property
    def key(self):
        """Shared source metric identity, independent of transcript assignment."""
        return (self.source, self.release, self.target.key, self.scope, self.method, self.metric)

    @property
    def is_supporting(self):
        """Positive reported support; zero and unknown both fail selection."""
        return self.value is not None and self.value > 0

    def validate_orf(self, model):
        """Reject attachment to another species, protein sequence or native ORF."""
        if self.target.reference.taxon_id != model.taxon_id:
            raise ValueError("Summary and ORF species disagree")
        reference = self.target.reference
        if reference.genome_build is not None and reference.genome_build != model.genome_build:
            raise ValueError("Summary and ORF assemblies disagree")
        for field in ("strain", "strain_taxon_id"):
            source = getattr(reference, field)
            destination = getattr(model.reference_context, field, None)
            if source != destination:
                raise ValueError("Summary and ORF strain contexts disagree")
        if self.protein_sha256 != hashlib.sha256(model.protein_sequence.encode("ascii")).hexdigest():
            raise ValueError("Summary and ORF protein sequences disagree")
        if self.scope == "orf" and (
            self.target.identifier != model.orf_id
            or reference.genome_build != model.genome_build
            or (model.reference_context is not None and reference.key != model.reference_context.key)
        ):
            raise ValueError("Summary and native ORF identities disagree")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, record):
        record = dict(record)
        target = dict(record["target"])
        target["reference"] = ReferenceContext(**target["reference"])
        record["target"] = ScopedTarget(**target)
        return cls(**record)


def unique_summaries(summaries):
    """Preserve first occurrence; conflicting facts for one identity are errors."""
    result = {}
    for summary in summaries:
        if not isinstance(summary, EvidenceSummary):
            raise TypeError("Expected EvidenceSummary")
        previous = result.setdefault(summary.key, summary)
        if previous != replace(summary, target=previous.target):
            raise ValueError("Conflicting evidence summaries for one source metric")
    return tuple(result.values())
