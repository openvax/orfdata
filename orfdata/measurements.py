"""Immutable quantitative observations, separate from ORF existence evidence."""

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from types import MappingProxyType
from typing import Optional

from .reference import ReferenceContext, ScopedTarget

ASSAY_KINDS = ("rna_seq", "ribo_seq", "digest_ms", "mhc_ms", "targeted_ms")
QUANTITIES = (
    "rna_abundance", "ribosome_occupancy", "translation_efficiency",
    "protein_abundance", "peptide_abundance", "spectrum_count", "turnover",
    "protein_rna_ratio", "abundance_ratio",
)
MEASUREMENT_STATUSES = ("measured", "below_detection", "not_measured", "not_mapped", "missing")
MEASUREMENT_ORIGINS = ("observed", "derived", "predicted")


def _freeze(value):
    """Copy JSON metadata into immutable containers; reject non-JSON values."""
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("Metadata keys must be strings")
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float) and math.isfinite(value):
        return value
    raise TypeError("Metadata must contain finite JSON values")


def _plain(value):
    if is_dataclass(value):
        return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def _identifier(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")


def _metadata(instance, *names):
    for name in names:
        value = getattr(instance, name)
        if not isinstance(value, Mapping):
            raise TypeError(f"{name} must be a metadata mapping")
        object.__setattr__(instance, name, _freeze(value))


def _optional_text(instance, *names):
    for name in names:
        value = getattr(instance, name)
        if value is not None:
            _identifier(value, name)


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")


class _Record:
    def to_dict(self):
        """Return a mutable JSON copy without exposing stored metadata."""
        return _plain(self)

    @classmethod
    def from_dict(cls, data):
        """Reconstruct a record, validating nested reference and target identities."""
        values = dict(data)
        if "reference" in values:
            values["reference"] = ReferenceContext(**values["reference"])
        if "target" in values:
            target = dict(values["target"])
            target["reference"] = ReferenceContext(**target["reference"])
            values["target"] = ScopedTarget(**target)
        if "targets" in values:
            targets = []
            for row in values["targets"]:
                row = dict(row)
                row["reference"] = ReferenceContext(**row["reference"])
                targets.append(ScopedTarget(**row))
            values["targets"] = tuple(targets)
        if "input_references" in values:
            values["input_references"] = tuple(MeasurementReference.from_dict(row)
                                               for row in values["input_references"])
        return cls(**values)


@dataclass(frozen=True)
class SourceAsset(_Record):
    """A pinned original input; hashes and sizes describe its raw bytes."""

    asset_id: str
    url: str
    sha256: str
    size: int
    access_terms: str
    license_url: Optional[str] = None
    redistributable: bool = False
    metadata: Mapping = field(default_factory=dict)

    def __post_init__(self):
        for name in ("asset_id", "url", "access_terms"):
            _identifier(getattr(self, name), name)
        if not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
            raise ValueError("Source assets require a lowercase SHA-256")
        if type(self.size) is not int or self.size < 0:
            raise ValueError("Source asset size must be a nonnegative integer")
        if type(self.redistributable) is not bool:
            raise TypeError("redistributable must be an explicit boolean")
        _optional_text(self, "license_url")
        _metadata(self, "metadata")


@dataclass(frozen=True)
class DatasetManifest(_Record):
    """Provenance and measurable coverage for a concrete quantitative release."""

    dataset_id: str
    version: str
    title: str
    references: tuple
    sources: tuple
    citations: tuple
    assay_kinds: tuple
    coverage: Mapping
    processing: Mapping
    schema: str = "orfdata.measurements.v1"

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["references"] = tuple(ReferenceContext(**c) for c in values["references"])
        values["sources"] = tuple(SourceAsset.from_dict(s) for s in values["sources"])
        return cls(**values)

    def __post_init__(self):
        for name in ("dataset_id", "version", "title"):
            _identifier(getattr(self, name), name)
        for name in ("dataset_id", "version"):
            value = getattr(self, name)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
                raise ValueError(f"{name} must be a portable dataset/path identifier")
        if self.schema != "orfdata.measurements.v1":
            raise ValueError("Unsupported quantitative dataset schema")
        for name in ("references", "sources", "citations", "assay_kinds"):
            if isinstance(getattr(self, name), str):
                raise TypeError(f"Dataset {name} must be a sequence of entries, not a string")
            values = tuple(getattr(self, name))
            if not values:
                raise ValueError(f"Dataset {name} cannot be empty")
            object.__setattr__(self, name, values)
        if any(not isinstance(c, ReferenceContext) for c in self.references):
            raise TypeError("Dataset references must be ReferenceContexts")
        if len({c.key for c in self.references}) != len(self.references):
            raise ValueError("Duplicate dataset reference contexts")
        if any(not isinstance(s, SourceAsset) for s in self.sources):
            raise TypeError("Dataset sources must be SourceAssets")
        if len({s.asset_id for s in self.sources}) != len(self.sources):
            raise ValueError("Duplicate source assets")
        if any(kind not in ASSAY_KINDS for kind in self.assay_kinds):
            raise ValueError("Unknown dataset assay kind")
        for citation in self.citations:
            _identifier(citation, "citation")
        if not self.coverage or not self.processing:
            raise ValueError("Dataset measurable coverage and processing provenance are required")
        _metadata(self, "coverage", "processing")


@dataclass(frozen=True)
class Sample(_Record):
    """One source-defined biological sample or pool within a dataset release.

    ``reference`` identifies the coordinate/annotation system used for mapping;
    sample strain/breed belongs in ``attributes`` and need not equal the strain
    from which the reference genome was assembled. Pool members are source IDs,
    never inferred from the number of runs or channels.
    """

    sample_id: str
    reference: ReferenceContext
    biological_unit_id: Optional[str] = None
    tissue: Optional[str] = None
    cell_type: Optional[str] = None
    condition: Optional[str] = None
    collection_time: Optional[str] = None
    biological_replicate: Optional[str] = None
    pool_members: tuple = ()
    attributes: Mapping = field(default_factory=dict)
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        _identifier(self.sample_id, "sample_id")
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Sample reference must be a ReferenceContext")
        _optional_text(self, "biological_unit_id", "tissue", "cell_type", "condition",
                       "collection_time", "biological_replicate")
        if isinstance(self.pool_members, str):
            raise TypeError("Pool members must be a sequence of IDs")
        members = tuple(self.pool_members)
        for member in members:
            _identifier(member, "pool member")
        if len(set(members)) != len(members):
            raise ValueError("Pool members must be unique")
        object.__setattr__(self, "pool_members", members)
        _metadata(self, "attributes", "provenance")


@dataclass(frozen=True)
class Assay(_Record):
    """An acquisition/quantification assay, distinct from biological replication."""

    assay_id: str
    sample_id: str
    kind: str
    reference: ReferenceContext
    run_id: Optional[str] = None
    technical_replicate: Optional[str] = None
    fraction: Optional[str] = None
    channel: Optional[str] = None
    plex: Optional[str] = None
    bridge: Optional[str] = None
    search_reference: Optional[str] = None
    normalization: Mapping = field(default_factory=dict)
    calibration: Mapping = field(default_factory=dict)
    attributes: Mapping = field(default_factory=dict)
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        _identifier(self.assay_id, "assay_id")
        _identifier(self.sample_id, "sample_id")
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Assay reference must be a ReferenceContext")
        if self.kind not in ASSAY_KINDS:
            raise ValueError("Unknown assay kind")
        _optional_text(self, "run_id", "technical_replicate", "fraction", "channel", "plex",
                       "bridge", "search_reference")
        _metadata(self, "normalization", "calibration", "attributes", "provenance")


@dataclass(frozen=True)
class MeasurementReference(_Record):
    """An input observation identified across dataset releases and references."""

    dataset_id: str
    dataset_version: str
    reference: ReferenceContext
    measurement_id: str

    def __post_init__(self):
        for name in ("dataset_id", "dataset_version", "measurement_id"):
            _identifier(getattr(self, name), name)
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Measurement reference requires a ReferenceContext")

    @property
    def key(self):
        return self.dataset_id, self.dataset_version, self.reference.key, self.measurement_id

    @classmethod
    def from_measurement(cls, measurement):
        if not isinstance(measurement, Measurement):
            raise TypeError("Expected a Measurement")
        return cls(measurement.dataset_id, measurement.dataset_version,
                   measurement.target.reference, measurement.measurement_id)


@dataclass(frozen=True)
class Measurement(_Record):
    """A native value with units, missingness, target resolution and provenance.

    No unit conversion is implicit. A measured zero remains a zero; censored or
    missing measurements have value=None and an explicit status. Preserve any
    sentinel/native field in provenance, not as a fabricated numeric measurement.
    Spectrum counts are their own quantity and never protein abundance.
    """

    measurement_id: str
    target: ScopedTarget
    sample_id: str
    assay_id: str
    quantity: str
    value: Optional[float]
    unit: str
    dataset_id: str
    dataset_version: str
    status: str = "measured"
    origin: str = "observed"
    scale: str = "linear"
    log_base: Optional[float] = None
    uncertainty: Optional[float] = None
    uncertainty_type: Optional[str] = None
    detection_limit: Optional[float] = None
    input_measurements: tuple = ()
    derivation: Mapping = field(default_factory=dict)
    provenance: Mapping = field(default_factory=dict)
    input_references: tuple = ()

    def __post_init__(self):
        for name in ("measurement_id", "sample_id", "assay_id", "unit", "dataset_id", "dataset_version"):
            _identifier(getattr(self, name), name)
        if not isinstance(self.target, ScopedTarget):
            raise TypeError("Measurement target must be a ScopedTarget")
        if self.quantity not in QUANTITIES:
            raise ValueError("Unknown measurement quantity")
        if self.status not in MEASUREMENT_STATUSES or self.origin not in MEASUREMENT_ORIGINS:
            raise ValueError("Unknown measurement status or origin")
        if self.status == "measured":
            _number(self.value, "value")
        elif self.value is not None:
            raise ValueError("Unmeasured/censored values must be None; retain native sentinels in provenance")
        if self.scale not in ("linear", "log"):
            raise ValueError("Scale must be linear or log")
        if self.scale == "log":
            _number(self.log_base, "log_base")
            if self.log_base <= 0 or self.log_base == 1:
                raise ValueError("Logarithm base must be positive and not one")
        elif self.log_base is not None:
            raise ValueError("Linear values cannot have a log base")
        if self.uncertainty is not None:
            _number(self.uncertainty, "uncertainty")
            _identifier(self.uncertainty_type, "uncertainty_type")
            if self.uncertainty < 0 or not self.uncertainty_type:
                raise ValueError("Uncertainty needs a nonnegative value and a named definition")
        elif self.uncertainty_type is not None:
            raise ValueError("Uncertainty type requires an uncertainty value")
        if self.detection_limit is not None:
            _number(self.detection_limit, "detection_limit")
        if isinstance(self.input_measurements, str):
            raise TypeError("Input measurements must be a sequence of IDs")
        inputs = tuple(self.input_measurements)
        for identifier in inputs:
            _identifier(identifier, "input measurement")
        if len(set(inputs)) != len(inputs) or self.measurement_id in inputs:
            raise ValueError("Derived input IDs must be unique and exclude the output")
        object.__setattr__(self, "input_measurements", inputs)
        references = tuple(self.input_references)
        if inputs and references:
            raise ValueError("Choose legacy same-context input IDs or explicit input references, not both")
        if any(not isinstance(value, MeasurementReference) for value in references):
            raise TypeError("input_references must contain MeasurementReferences")
        own_key = (self.dataset_id, self.dataset_version, self.target.reference.key, self.measurement_id)
        keys = {value.key for value in references}
        if len(keys) != len(references) or own_key in keys:
            raise ValueError("Input references must be unique and exclude the output")
        for value in references:
            if (value.dataset_id, value.dataset_version) != (self.dataset_id, self.dataset_version):
                raise ValueError("Derived inputs must belong to the same dataset release")
            if value.reference.taxon_id != self.target.reference.taxon_id:
                raise ValueError("Derived abundance inputs cannot cross species")
        object.__setattr__(self, "input_references", references)
        if self.origin in ("derived", "predicted") and not self.derivation:
            raise ValueError("Derived/predicted values require a named method and assumptions")
        if self.origin in ("derived", "predicted") and isinstance(self.derivation, Mapping):
            _identifier(self.derivation.get("method"), "derivation method")
        if self.origin == "observed" and (inputs or references or self.derivation):
            raise ValueError("An observed value cannot contain a derivation")
        _metadata(self, "derivation", "provenance")


@dataclass(frozen=True)
class TargetMapping(_Record):
    """Source identity and all supported target assignments, including failures."""

    mapping_id: str
    source_namespace: str
    source_identifier: str
    reference: ReferenceContext
    targets: tuple = ()
    status: str = "unmapped"
    method: Optional[str] = None
    reason: Optional[str] = None
    evidence: Mapping = field(default_factory=dict)
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        for name in ("mapping_id", "source_namespace", "source_identifier"):
            _identifier(getattr(self, name), name)
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Mapping reference must be a ReferenceContext")
        _optional_text(self, "method", "reason")
        targets = tuple(self.targets)
        if any(not isinstance(t, ScopedTarget) for t in targets):
            raise TypeError("Mapping targets must be ScopedTargets")
        if any(t.reference.key != self.reference.key for t in targets):
            raise ValueError("Mappings cannot cross reference contexts; orthology/liftover is separate")
        if len({t.key for t in targets}) != len(targets):
            raise ValueError("Duplicate mapping targets")
        if self.status not in ("unique", "ambiguous", "unmapped", "unresolved"):
            raise ValueError("Unknown mapping status")
        if ((self.status == "unique" and len(targets) != 1)
                or (self.status == "ambiguous" and len(targets) < 2)
                or (self.status == "unmapped" and targets)):
            raise ValueError("Mapping status disagrees with candidate targets")
        if self.status in ("unmapped", "unresolved") and not self.reason:
            raise ValueError("An unmapped or unresolved identifier requires a reason")
        if self.status != "unmapped" and not self.method:
            raise ValueError("Mapped identifiers require an attribution method")
        object.__setattr__(self, "targets", targets)
        _metadata(self, "evidence", "provenance")
