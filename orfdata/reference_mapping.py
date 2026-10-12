"""Explicit mappings between native source and annotated reference targets."""

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Optional

from .measurements import Measurement, _identifier, _metadata, _Record
from .reference import ReferenceContext, ScopedTarget


def _target(data):
    values = dict(data)
    values["reference"] = ReferenceContext(**values["reference"])
    return ScopedTarget(**values)


@dataclass(frozen=True)
class ReferenceMapping(_Record):
    """All candidate assignments for a native target in another reference.

    Unique identifier resolution does not imply sequence validation. A shared
    group can identify candidate proteins while its measurement remains a group
    observation. Parent RNA relationships belong in TargetRelation, not here.
    """

    mapping_id: str
    dataset_id: str
    dataset_version: str
    source: ScopedTarget
    target_reference: ReferenceContext
    targets: tuple
    status: str
    method: str
    validation: Mapping
    reason: Optional[str] = None
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        for name in ("mapping_id", "dataset_id", "dataset_version", "method"):
            _identifier(getattr(self, name), name)
        if not isinstance(self.source, ScopedTarget) or not isinstance(self.target_reference, ReferenceContext):
            raise TypeError("Reference mappings require explicit source and destination contexts")
        if self.source.reference.taxon_id != self.target_reference.taxon_id:
            raise ValueError("Reference remapping cannot cross species; use orthology")
        targets = tuple(self.targets)
        if any(not isinstance(target, ScopedTarget) for target in targets):
            raise TypeError("Mapping targets must be ScopedTargets")
        if any(target.reference.key != self.target_reference.key for target in targets):
            raise ValueError("Candidate targets must share the declared destination reference")
        if len({target.key for target in targets}) != len(targets):
            raise ValueError("Duplicate reference mapping targets")
        if self.status not in ("unique", "ambiguous", "unmapped", "unresolved"):
            raise ValueError("Unknown reference mapping status")
        if ((self.status == "unique" and len(targets) != 1)
                or (self.status == "ambiguous" and len(targets) < 2)
                or (self.status == "unmapped" and targets)):
            raise ValueError("Reference mapping status disagrees with candidate count")
        if self.reason is not None:
            _identifier(self.reason, "mapping reason")
        if self.status in ("unmapped", "unresolved") and not self.reason:
            raise ValueError("Unmapped/unresolved assignments require a reason")
        if not self.validation:
            raise ValueError("Mapping validation and its limits must be explicit")
        object.__setattr__(self, "targets", targets)
        _metadata(self, "validation", "provenance")
        for target in targets:
            kinds = self.source.kind, target.kind
            allowed = (kinds[0] == kinds[1] or kinds in (("protein", "orf"), ("orf", "protein"))
                       or (kinds[0] == "protein_group" and kinds[1] in ("protein", "orf")))
            if not allowed:
                raise ValueError("Reference mapping cannot change target resolution; use a context relationship")
            if self.source.kind == "protein_group" and target.kind != "protein_group" and self.status == "unique":
                raise ValueError("A protein-group observation cannot become a unique member measurement")
            if "orf" in kinds and not (self.sequence_validated and self.validation.get("frame")
                                        and self.validation.get("position")):
                raise ValueError("ORF remapping requires sequence, frame and position validation")
        for name in ("strain", "strain_taxon_id"):
            source = getattr(self.source.reference, name)
            destination = getattr(self.target_reference, name)
            if source is not None and destination is not None and source != destination:
                if not self.validation.get("strain_relationship"):
                    raise ValueError("Different reference strains require explicit relationship evidence")

    @property
    def sequence_validated(self):
        """Whether all candidate translations are declared exact sequence matches.

        Importers must actually check every candidate against the pinned source;
        this records that validation, rather than running a sequence search.
        """
        if not self.targets:
            return False
        evidence = self.validation.get("sequence", {})
        if not isinstance(evidence, Mapping):
            return False
        source, target = evidence.get("source_sha256"), evidence.get("target_sha256")
        return (isinstance(source, str) and re.fullmatch(r"[0-9a-f]{64}", source) is not None
                and source == target and evidence.get("all_candidates_checked") is True)

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["source"] = _target(values["source"])
        values["target_reference"] = ReferenceContext(**values["target_reference"])
        values["targets"] = tuple(_target(target) for target in values["targets"])
        return cls(**values)


@dataclass(frozen=True)
class MappedMeasurement(_Record):
    """An original observation plus a requested target and its explicit mapping.

    There is deliberately no replacement value/target: inspect ``measurement``
    for native identity, unit, dataset, source and sample/assay scope. Mapping a
    group to a candidate member does not quantify that member independently.
    """

    measurement: Measurement
    mapping: ReferenceMapping
    requested_target: ScopedTarget

    def __post_init__(self):
        if not isinstance(self.measurement, Measurement) or not isinstance(self.mapping, ReferenceMapping):
            raise TypeError("Mapped views require an original Measurement and ReferenceMapping")
        if not isinstance(self.requested_target, ScopedTarget):
            raise TypeError("Mapped views require an exact requested target")
        if self.measurement.target.key != self.mapping.source.key:
            raise ValueError("Observation and mapping source identities disagree")
        if (self.measurement.dataset_id, self.measurement.dataset_version) != (
                self.mapping.dataset_id, self.mapping.dataset_version):
            raise ValueError("Observation and mapping belong to different dataset releases")
        if self.requested_target.key not in {target.key for target in self.mapping.targets}:
            raise ValueError("Requested target is not among the declared candidates")

    @property
    def attribution(self):
        if self.measurement.target.kind == "protein_group" and self.requested_target.kind != "protein_group":
            return "group_membership_context"
        return self.mapping.status

    def to_dict(self):
        result = super().to_dict()
        result["attribution"] = self.attribution
        return result

    @classmethod
    def from_dict(cls, data):
        return cls(Measurement.from_dict(data["measurement"]), ReferenceMapping.from_dict(data["mapping"]),
                   _target(data["requested_target"]))
