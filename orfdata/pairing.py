"""Explicit assay relationships and comparisons of original observations."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Optional

from .measurements import Measurement, SourceAsset, _identifier, _metadata, _number, _Record
from .reference import ReferenceContext, ScopedTarget

PAIRING_BASES = (
    "same_material", "adjacent_donor_tissue", "biological_replicate",
    "independent_cultures", "tissue_aggregate", "condition_aggregate",
)
TEMPORAL_RELATIONS = ("concurrent", "longitudinal", "unknown")
TARGET_RELATIONS = ("same_target", "canonical_translation", "parent_rna_context", "shared_group_member")


def _target(data):
    values = dict(data)
    values["reference"] = ReferenceContext(**values["reference"])
    return ScopedTarget(**values)


def _release(record):
    for name in ("dataset_id", "dataset_version"):
        _identifier(getattr(record, name), name)


def _check_one_to_one(seen, source_key, target_key, one_to_one):
    """Check a claimed one-to-one edge against every supplied candidate edge."""
    targets, strict = seen.setdefault(source_key, (set(), False))
    targets.add(target_key)
    strict = strict or one_to_one
    if strict and len(targets) > 1:
        raise ValueError("Conflicting one-to-one target relationships")
    seen[source_key] = (targets, strict)


@dataclass(frozen=True)
class AssayPair(_Record):
    """A source-backed assay pairing, including an explicitly missing mate.

    The pairing basis and temporal relationship are independent. A shared tissue
    label supplies neither. Replicate/aliquot identity remains in the original
    samples and assays, rather than being counted from the number of pair rows.
    """

    pair_id: str
    reference: ReferenceContext
    dataset_id: str
    dataset_version: str
    left_assay_id: Optional[str]
    right_assay_id: Optional[str]
    basis: str
    evidence: Mapping
    temporal_relation: str = "unknown"
    interval: Optional[float] = None
    interval_unit: Optional[str] = None

    def __post_init__(self):
        _identifier(self.pair_id, "pair_id")
        _release(self)
        if not isinstance(self.reference, ReferenceContext):
            raise TypeError("Pair requires a ReferenceContext")
        for identifier in (self.left_assay_id, self.right_assay_id):
            if identifier is not None:
                _identifier(identifier, "assay ID")
        if self.left_assay_id is None and self.right_assay_id is None:
            raise ValueError("A pairing requires at least one declared assay")
        if self.left_assay_id == self.right_assay_id:
            raise ValueError("A pairing requires distinct assays")
        if self.basis not in PAIRING_BASES:
            raise ValueError("Unknown assay pairing basis")
        if self.temporal_relation not in TEMPORAL_RELATIONS:
            raise ValueError("Unknown temporal relationship")
        if not self.evidence:
            raise ValueError("Pairing requires source evidence")
        if self.interval is not None:
            _number(self.interval, "interval")
            _identifier(self.interval_unit, "interval_unit")
            if self.temporal_relation == "concurrent" and self.interval != 0:
                raise ValueError("Concurrent assays cannot have a nonzero interval")
        elif self.interval_unit is not None:
            raise ValueError("Interval units require a value")
        _metadata(self, "evidence")


@dataclass(frozen=True)
class TargetRelation(_Record):
    """Explain why two exact targets can be compared without changing either ID."""

    relation_id: str
    dataset_id: str
    dataset_version: str
    left: ScopedTarget
    right: ScopedTarget
    relation: str
    cardinality: str
    evidence: Mapping

    def __post_init__(self):
        _identifier(self.relation_id, "relation_id")
        _release(self)
        if not isinstance(self.left, ScopedTarget) or not isinstance(self.right, ScopedTarget):
            raise TypeError("Target relations require ScopedTargets")
        if self.left.reference.key != self.right.reference.key:
            raise ValueError("Within-dataset target relations cannot cross reference contexts")
        if self.relation not in TARGET_RELATIONS:
            raise ValueError("Unknown target relationship")
        if self.cardinality not in ("one_to_one", "one_to_many", "many_to_one", "many_to_many", "unknown"):
            raise ValueError("Unknown target relationship cardinality")
        if not self.evidence:
            raise ValueError("Target relations require mapping evidence")
        kinds = (self.left.kind, self.right.kind)
        if self.relation == "same_target" and self.left.key != self.right.key:
            raise ValueError("same_target requires identical scoped target identities")
        if self.relation == "same_target" and self.cardinality != "one_to_one":
            raise ValueError("An identical target has one-to-one cardinality")
        if self.relation == "canonical_translation" and kinds not in (("gene", "protein"), ("transcript", "protein")):
            raise ValueError("Canonical translation links a gene/transcript to a protein")
        if self.relation == "parent_rna_context" and kinds not in (("gene", "orf"), ("transcript", "orf")):
            raise ValueError("Parent RNA context links a gene/transcript to an ORF")
        if self.relation == "shared_group_member" and kinds not in (("protein_group", "protein"), ("protein_group", "orf")):
            raise ValueError("Group membership links a group to one candidate member")
        _metadata(self, "evidence")

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["left"], values["right"] = _target(values["left"]), _target(values["right"])
        return cls(**values)


@dataclass(frozen=True)
class MatchedMeasurements(_Record):
    """One declared relationship with all original values and explicit missingness.

    Tuples avoid multiplying technical/quantity records into pseudo-replicates.
    For parent_rna_context, left values remain measurements of the parent RNA.
    They are never measurements of the ORF named on the right.
    """

    pairing: AssayPair
    relation: TargetRelation
    left: tuple
    right: tuple

    @classmethod
    def from_dict(cls, data):
        return cls(AssayPair.from_dict(data["pairing"]), TargetRelation.from_dict(data["relation"]),
                   tuple(Measurement.from_dict(value) for value in data["left"]),
                   tuple(Measurement.from_dict(value) for value in data["right"]))

    def __post_init__(self):
        pair, relation = self.pairing, self.relation
        if not isinstance(pair, AssayPair) or not isinstance(relation, TargetRelation):
            raise TypeError("A match requires an AssayPair and TargetRelation")
        if (pair.dataset_id, pair.dataset_version) != (relation.dataset_id, relation.dataset_version):
            raise ValueError("Pairing and target relation belong to different dataset releases")
        if pair.reference.key != relation.left.reference.key:
            raise ValueError("Pairing and target relation have different references")
        for side in ("left", "right"):
            values = tuple(getattr(self, side))
            target = getattr(relation, side)
            assay_id = getattr(pair, side + "_assay_id")
            seen = set()
            for value in values:
                if not isinstance(value, Measurement):
                    raise TypeError("Matched values must be Measurements")
                if (value.dataset_id, value.dataset_version) != (pair.dataset_id, pair.dataset_version):
                    raise ValueError("Matched measurement belongs to a different dataset release")
                if value.target.key != target.key or value.assay_id != assay_id:
                    raise ValueError("Matched observation disagrees with its declared target/assay")
                if value.measurement_id in seen:
                    raise ValueError("Duplicate matched measurement")
                seen.add(value.measurement_id)
            object.__setattr__(self, side, values)

    def to_table_row(self):
        """Return named assay-quantity columns with original units and provenance.

        Columns contain lists, preserving multiple measurements without silently
        averaging them. Empty columns are missing observations, not measured zero.
        Pairing and target relationships remain present to qualify comparisons.
        """
        columns = {name: [] for name in (
            "rna_abundance", "ribosome_occupancy", "translation_efficiency",
            "protein_abundance", "peptide_abundance", "spectrum_count", "other")}
        for value in self.left + self.right:
            key = value.quantity if value.quantity in columns else "other"
            columns[key].append(value.to_dict())
        return dict(dataset_id=self.pairing.dataset_id, dataset_version=self.pairing.dataset_version,
                    pairing=self.pairing.to_dict(), target_relation=self.relation.to_dict(), **columns)

    def ratio(self, numerator, denominator, *, measurement_id, method, normalization):
        """Derive a named native-unit ratio from explicitly selected observations.

        No implicit log conversion, replicate aggregation, unit conversion or
        uncertainty propagation. Undefined ratios are not_measured with a reason.
        Parent RNA and ambiguous group membership cannot produce ORF ratios.
        """
        _identifier(method, "ratio method")
        if not isinstance(normalization, Mapping) or not normalization:
            raise ValueError("Ratio normalization and assumptions must be explicit")
        if self.relation.relation not in ("same_target", "canonical_translation"):
            raise ValueError("Parent RNA/group membership cannot establish target-specific ratios")
        if self.relation.cardinality != "one_to_one":
            raise ValueError("A target-specific ratio requires one-to-one attribution")
        if not ((numerator in self.left and denominator in self.right)
                or (numerator in self.right and denominator in self.left)):
            raise ValueError("Ratio inputs must be original measurements on opposite sides")
        if numerator.scale != "linear" or denominator.scale != "linear":
            raise ValueError("Ratio inputs must be linear; any log conversion must be explicit")
        if (numerator.quantity, denominator.quantity) == ("protein_abundance", "rna_abundance"):
            quantity = "protein_rna_ratio"
        elif (numerator.quantity, denominator.quantity) == ("ribosome_occupancy", "rna_abundance"):
            quantity = "translation_efficiency"
        elif numerator.quantity == denominator.quantity:
            quantity = "abundance_ratio"
        else:
            raise ValueError("Unsupported ratio quantity relationship")
        reason = None
        if numerator.status != "measured" or denominator.status != "measured":
            reason = "input_not_measured"
        elif denominator.value <= 0:
            reason = "nonpositive_denominator"
        derivation = dict(method=method, normalization=dict(normalization),
                          pairing=self.pairing.to_dict(), target_relation=self.relation.to_dict(),
                          uncertainty="not_propagated", unresolved_reason=reason)
        return Measurement(
            measurement_id=measurement_id, target=numerator.target,
            sample_id=numerator.sample_id, assay_id=numerator.assay_id, quantity=quantity,
            value=None if reason else numerator.value / denominator.value,
            unit=f"({numerator.unit})/({denominator.unit})",
            dataset_id=numerator.dataset_id, dataset_version=numerator.dataset_version,
            status="not_measured" if reason else "measured", origin="derived",
            input_measurements=(numerator.measurement_id, denominator.measurement_id),
            derivation=derivation,
        )


@dataclass(frozen=True)
class OrthologyMapping(_Record):
    """Pinned cross-species mapping; it never transfers measurements or evidence."""

    mapping_id: str
    left: ScopedTarget
    right: ScopedTarget
    source: SourceAsset
    source_release: str
    cardinality: str
    confidence: Mapping
    evidence: Mapping

    def __post_init__(self):
        _identifier(self.mapping_id, "mapping_id")
        _identifier(self.source_release, "source_release")
        if not isinstance(self.left, ScopedTarget) or not isinstance(self.right, ScopedTarget):
            raise TypeError("Orthology requires ScopedTargets")
        if self.left.reference.taxon_id == self.right.reference.taxon_id:
            raise ValueError("Orthology comparison requires different species")
        if not isinstance(self.source, SourceAsset):
            raise TypeError("Orthology source must have a pinned SourceAsset")
        if self.cardinality not in ("one_to_one", "one_to_many", "many_to_one", "many_to_many"):
            raise ValueError("Unknown orthology cardinality")
        if (self.left.kind, self.right.kind) not in (("gene", "gene"), ("orf", "orf")):
            raise ValueError("Orthology maps canonical genes or independently supported ORFs")
        if not self.evidence or not self.confidence:
            raise ValueError("Orthology evidence and source confidence definition are required")
        if self.left.kind == "orf" and not all(self.evidence.get(key) for key in
                                              ("sequence", "frame", "position")):
            raise ValueError("ORF orthology requires its own sequence, frame and position evidence")
        _metadata(self, "confidence", "evidence")

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["left"], values["right"] = _target(values["left"]), _target(values["right"])
        values["source"] = SourceAsset.from_dict(values["source"])
        return cls(**values)


@dataclass(frozen=True)
class OrthologyComparison(_Record):
    mapping: OrthologyMapping
    left: tuple
    right: tuple

    @classmethod
    def from_dict(cls, data):
        return cls(OrthologyMapping.from_dict(data["mapping"]),
                   tuple(Measurement.from_dict(value) for value in data["left"]),
                   tuple(Measurement.from_dict(value) for value in data["right"]))

    def __post_init__(self):
        if not isinstance(self.mapping, OrthologyMapping):
            raise TypeError("Comparison requires an OrthologyMapping")
        for side in ("left", "right"):
            values = tuple(getattr(self, side))
            target = getattr(self.mapping, side)
            if any(not isinstance(value, Measurement) or value.target.key != target.key for value in values):
                raise ValueError("Orthology observations must retain their original scoped targets")
            object.__setattr__(self, side, values)


def compare_orthologs(left_dataset, right_dataset, mappings, *,
                     allowed_cardinalities=("one_to_one",), allow_orfs=False,
                     left_filters=None, right_filters=None):
    """Return both species' original observations for explicitly opted-in mappings."""
    allowed = set(allowed_cardinalities)
    valid = {"one_to_one", "one_to_many", "many_to_one", "many_to_many"}
    if not allowed or allowed - valid:
        raise ValueError("Explicit valid orthology cardinalities are required")
    mappings = tuple(mappings)
    seen = ({}, {})
    for mapping in mappings:
        if not isinstance(mapping, OrthologyMapping):
            raise TypeError("Expected OrthologyMapping records")
        for index, (source, target) in enumerate(((mapping.left, mapping.right), (mapping.right, mapping.left))):
            _check_one_to_one(seen[index], (source.key, target.reference.key, target.namespace), target.key,
                              mapping.cardinality == "one_to_one")
    for mapping in mappings:
        if mapping.cardinality not in allowed or (mapping.left.kind == "orf" and not allow_orfs):
            continue
        yield OrthologyComparison(
            mapping,
            tuple(left_dataset.measurements(target=mapping.left, **(left_filters or {}))),
            tuple(right_dataset.measurements(target=mapping.right, **(right_filters or {}))),
        )
