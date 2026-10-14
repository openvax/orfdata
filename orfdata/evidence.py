"""Rules that select ORFs by Ribo-seq and MS evidence, separately or together.

``RiboSeq(min_studies=2)`` and ``MS(min_peptides=2, unique=True)`` are rules;
``a & b`` requires both and ``a | b`` accepts either. The evidence method
names used elsewhere (``"ms"``, ``"ribosome"``, ``"hla"`` ...) are shorthands
for the simplest rules. A rule tests what a source reports; passing one is not
proof that a stable protein exists.
"""

from dataclasses import dataclass
from typing import Optional

MS_METHODS = ("shotgun", "mhc", "hla", "ms_unclassified")
EVIDENCE_METHODS = ("ribosome",) + MS_METHODS + ("ms",)


class EvidenceRule:
    """A test on one ORF's evidence. Combine rules with ``&`` (both) and ``|`` (either)."""

    def __bool__(self):
        # ``rule_a and rule_b`` would silently keep only rule_b.
        raise TypeError("Combine evidence rules with & (both) or | (either), not 'and' or 'or'")

    def __and__(self, other):
        return AllOf((self, _checked(other)))

    def __or__(self, other):
        return AnyOf((self, _checked(other)))

    def matches(self, evidence):
        raise NotImplementedError

    def to_dict(self):
        raise NotImplementedError


def _checked(rule):
    if not isinstance(rule, EvidenceRule):
        raise TypeError("Combine evidence rules with other rules, such as RiboSeq() or MS()")
    return rule


def _count(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _flags(rule, *names):
    for name in names:
        if type(getattr(rule, name)) is not bool:
            raise TypeError(f"{name} must be True or False")


def _distinct(sequences):
    """Peptides counted once: I and L are indistinguishable by mass, and a peptide
    inside a longer one, such as an HLA length variant, is not separate support."""
    normalized = {s.replace("I", "L") for s in sequences}
    return {s for s in normalized if not any(s != t and s in t for t in normalized)}


@dataclass(frozen=True)
class RiboSeq(EvidenceRule):
    """Ribosome profiling in at least ``min_studies`` independent studies the source names.

    A source's aggregate count satisfies ``min_studies=1`` unless ``resolved``
    is set. It never stands for several studies, because its counting rules
    are the source's own.
    """

    min_studies: int = 1
    resolved: bool = False

    def __post_init__(self):
        _count(self.min_studies, "min_studies")
        _flags(self, "resolved")

    def matches(self, evidence):
        if len(set(evidence.ribosome_studies)) >= self.min_studies:
            return True
        return self.min_studies == 1 and not self.resolved and any(
            s.method == "ribosome" and s.is_supporting for s in evidence.summaries)

    def to_dict(self):
        return {"ribo_seq": dict(min_studies=self.min_studies, resolved=self.resolved)}


@dataclass(frozen=True)
class MS(EvidenceRule):
    """Mass spectrometry: at least ``min_peptides`` distinct supporting peptides.

    Peptides that differ only by I and L, or that lie inside a longer
    supporting peptide, count once. ``methods`` limits the peptide kinds (``"shotgun"``, ``"mhc"``, ``"hla"``,
    ``"ms_unclassified"``); ``"mhc"`` includes ``"hla"``. ``unique`` needs one
    protein mapping in the source, not uniqueness in every proteome.
    ``reviewed`` needs convincing reviewed spectra (benchmark scores 4-5).
    A source's aggregate count satisfies only ``min_peptides=1`` without
    ``unique``, ``reviewed`` or ``resolved``.
    """

    min_peptides: int = 1
    unique: bool = False
    reviewed: bool = False
    resolved: bool = False
    methods: Optional[tuple] = None

    def __post_init__(self):
        _count(self.min_peptides, "min_peptides")
        _flags(self, "unique", "reviewed", "resolved")
        if self.methods is not None:
            methods = {self.methods} if isinstance(self.methods, str) else set(self.methods)
            if not methods or methods - set(MS_METHODS):
                raise ValueError(f"MS methods must be drawn from {MS_METHODS}")
            if "mhc" in methods:  # mhc already includes human hla.
                methods.discard("hla")
            # One canonical form, so equal rules compare and serialize the same.
            object.__setattr__(self, "methods", tuple(m for m in MS_METHODS if m in methods))

    def _kind(self, method):
        return self.methods is None or method in self.methods or (method == "hla" and "mhc" in self.methods)

    def matches(self, evidence):
        peptides = _distinct(p.sequence for p in evidence.peptides
                             if p.is_supporting and self._kind(p.method)
                             and (not self.unique or p.unique_in_source)
                             and (not self.reviewed or p.quality == "reviewed_support"))
        if len(peptides) >= self.min_peptides:
            return True
        return (self.min_peptides == 1 and not (self.unique or self.reviewed or self.resolved) and any(
            s.method != "ribosome" and s.is_supporting and (self.methods is None or self._kind(s.method))
            for s in evidence.summaries))

    def to_dict(self):
        return {"ms": dict(min_peptides=self.min_peptides, unique=self.unique, reviewed=self.reviewed,
                           resolved=self.resolved, methods=None if self.methods is None else list(self.methods))}


@dataclass(frozen=True)
class AllOf(EvidenceRule):
    """Every rule must match."""

    rules: tuple

    def __post_init__(self):
        object.__setattr__(self, "rules", _flatten(AllOf, self.rules))

    def matches(self, evidence):
        return all(rule.matches(evidence) for rule in self.rules)

    def to_dict(self):
        return {"all": [rule.to_dict() for rule in self.rules]}


@dataclass(frozen=True)
class AnyOf(EvidenceRule):
    """At least one rule must match."""

    rules: tuple

    def __post_init__(self):
        object.__setattr__(self, "rules", _flatten(AnyOf, self.rules))

    def matches(self, evidence):
        return any(rule.matches(evidence) for rule in self.rules)

    def to_dict(self):
        return {"any": [rule.to_dict() for rule in self.rules]}


def _flatten(kind, rules):
    flat = []
    for rule in rules:
        flat.extend(rule.rules if isinstance(rule, kind) else (_checked(rule),))
    if len(flat) < 2:
        raise ValueError(f"{kind.__name__} combines two or more rules")
    return tuple(flat)


def evidence_rule(evidence, *, require_unique=False, require_reviewed=False, require_resolved=False):
    """The rule for ``evidence``: a rule, an evidence method name, or None for every ORF.

    A method name with ``require_*`` options becomes the equivalent rule:
    ``"ms"`` is ``MS()``, ``"hla"`` is ``MS(methods=("hla",))`` and
    ``"ribosome"`` is ``RiboSeq()``. Options that cannot apply raise.
    """
    if isinstance(evidence, EvidenceRule):
        if require_unique or require_reviewed or require_resolved:
            raise ValueError("Set unique, reviewed or resolved on the rule itself")
        return evidence
    if evidence is None:
        options = dict(require_unique=require_unique, require_reviewed=require_reviewed,
                       require_resolved=require_resolved)
        if any(options.values()):
            raise ValueError(f"{next(k for k, v in options.items() if v)} requires an evidence method")
        return None
    if evidence not in EVIDENCE_METHODS:
        raise ValueError(f"Unknown evidence method: {evidence}")
    if evidence == "ribosome":
        if require_unique or require_reviewed:
            raise ValueError("MS filters do not apply to ribosome evidence")
        return RiboSeq(resolved=require_resolved)
    return MS(unique=require_unique, reviewed=require_reviewed, resolved=require_resolved,
              methods=None if evidence == "ms" else (evidence,))
