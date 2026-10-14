"""Separate Ribo-seq and MS stringency, and the method names they generalize."""

import pytest
from test_evidence_summary import summary

from orfdata import MS, ORFEvidence, PeptideEvidence, RiboSeq, load_catalog
from orfdata.evidence import AllOf, AnyOf, evidence_rule


def legacy_has_evidence(evidence, method, require_unique=False, require_reviewed=False, require_resolved=False):
    """Single-method selection as released in 0.18.0, kept to check the shorthands against it."""
    if method == "ribosome":
        return bool(evidence.ribosome_studies) or (not require_resolved and any(
            s.method == "ribosome" and s.is_supporting for s in evidence.summaries))
    resolved = any(
        (method == "ms" or p.method == method or (method == "mhc" and p.method == "hla"))
        and p.is_supporting
        and (not require_unique or p.unique_in_source)
        and (not require_reviewed or p.quality == "reviewed_support")
        for p in evidence.peptides)
    return resolved or (not (require_resolved or require_unique or require_reviewed) and any(
        (method == "ms" and s.method != "ribosome" or s.method == method or method == "mhc" and s.method == "hla")
        and s.is_supporting for s in evidence.summaries))


def peptide(sequence, method="shotgun", **options):
    return PeptideEvidence(method, sequence, ((1, len(sequence)),), "synthetic", "fixture", **options)


@pytest.fixture(scope="module")
def catalog():
    return load_catalog()


def test_method_names_select_exactly_as_before(catalog):
    aggregates = [ORFEvidence(summaries=(summary(method=m, value=v),))
                  for m in ("ribosome", "ms", "shotgun", "mhc", "ms_unclassified") for v in (0, None, 3)]
    evidences = [r.evidence for r in catalog.records] + aggregates
    options = [dict(require_unique=u, require_reviewed=r, require_resolved=v)
               for u in (False, True) for r in (False, True) for v in (False, True)]
    cases = [("ribosome", dict(require_resolved=v)) for v in (False, True)] + [
        (m, o) for m in ("ms", "shotgun", "mhc", "hla", "ms_unclassified") for o in options]
    for method, option in cases:
        rule = evidence_rule(method, **option)
        assert [rule.matches(e) for e in evidences] == [legacy_has_evidence(e, method, **option) for e in evidences]


def test_separate_and_combined_stringency_on_the_bundled_catalogue(catalog):
    def count(rule):
        return len(catalog.query(evidence=rule))

    assert [count(RiboSeq(min_studies=n)) for n in (1, 2, 3)] == [7264, 3085, 1675]
    assert [count(MS(min_peptides=n)) for n in (1, 2, 3)] == [2186, 916, 417]
    assert (count(MS(unique=True)), count(MS(min_peptides=2, unique=True)), count(MS(reviewed=True))) == (1623, 612, 44)
    assert count(MS(methods="hla")) == count(MS(methods=("mhc",))) == 1881
    assert (count(MS(methods="shotgun")), count(MS(methods=("shotgun", "ms_unclassified")))) == (73, 610)
    both, either = count(RiboSeq(min_studies=2) & MS()), count(RiboSeq(min_studies=2) | MS())
    assert (both, either) == (1296, 3975) and both + either == 3085 + 2186
    assert count(RiboSeq(min_studies=2) & MS(min_peptides=2)) == 600
    assert count(RiboSeq(min_studies=3) | MS(min_peptides=2, unique=True)) == 2023


def test_peptides_count_once_and_aggregate_counts_stand_for_one_observation():
    repeated = ORFEvidence(peptides=(peptide("KPG"), peptide("KPG", "hla", mapping_count=1),
                                     peptide("KPGA", review_ratings=((2, 2),), quality="reviewed_low_quality")))
    assert MS().matches(repeated) and not MS(min_peptides=2).matches(repeated)
    assert MS(methods="mhc", unique=True).matches(repeated)
    assert not MS(methods="shotgun", unique=True).matches(repeated)
    studies = ORFEvidence(ribosome_studies=("a", "b", "a"))
    assert RiboSeq(min_studies=2).matches(studies) and not RiboSeq(min_studies=3).matches(studies)
    aggregate = ORFEvidence(summaries=(summary(method="ribosome", value=5), summary(method="shotgun", value=3)))
    assert RiboSeq().matches(aggregate) and MS().matches(aggregate) and MS(methods="shotgun").matches(aggregate)
    for rule in (RiboSeq(min_studies=2), RiboSeq(resolved=True), MS(min_peptides=2), MS(unique=True),
                 MS(resolved=True), MS(methods="hla"), RiboSeq() & MS(min_peptides=2)):
        assert not rule.matches(aggregate), rule


def test_rules_combine_flatten_and_serialize():
    strict = RiboSeq(min_studies=2) & MS(unique=True) & MS(methods="shotgun")
    assert isinstance(strict, AllOf) and len(strict.rules) == 3
    lenient = RiboSeq() | MS() | MS(reviewed=True)
    assert isinstance(lenient, AnyOf) and len(lenient.rules) == 3
    assert (RiboSeq(min_studies=2) | MS(methods="hla")).to_dict() == {"any": [
        {"ribo_seq": {"min_studies": 2, "resolved": False}},
        {"ms": {"min_peptides": 1, "unique": False, "reviewed": False, "resolved": False, "methods": ["hla"]}}]}
    assert MS(methods="hla") == MS(methods=("hla",))
    for bad in (lambda: RiboSeq(min_studies=0), lambda: MS(min_peptides=True), lambda: MS(methods=("ribosome",)),
                lambda: MS(methods=()), lambda: RiboSeq() & "ms", lambda: AllOf((RiboSeq(),))):
        with pytest.raises((ValueError, TypeError)):
            bad()
    with pytest.raises(ValueError, match="on the rule itself"):
        evidence_rule(MS(), require_unique=True)
    with pytest.raises(ValueError, match="do not apply to ribosome"):
        evidence_rule("ribosome", require_unique=True)
    assert evidence_rule(None) is None and evidence_rule("hla", require_unique=True) == MS(unique=True, methods="hla")
