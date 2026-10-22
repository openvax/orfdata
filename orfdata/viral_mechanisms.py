"""Curated virus--tumour mechanism selectors for the GDC viral collection.

These selectors describe evidence scoped to a virus and a relevant tumour
setting. They are not patient-level treatment recommendations: an integration
event, a viral transcript, and tumour dependence on a viral product must each
be measured in the sample under study.
"""

from dataclasses import dataclass


_IARC_GROUP_1_HPV = frozenset((
    "HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV39", "HPV45", "HPV51", "HPV52",
    "HPV56", "HPV58", "HPV59",
))


@dataclass(frozen=True)
class ViralMechanism:
    """One source-backed mechanism selector for relevant tumour settings."""

    name: str
    title: str
    description: str
    abbreviations: tuple[str, ...]
    source_urls: tuple[str, ...]
    limitations: str

    @property
    def references(self):
        """Number of GDC reference rows selected by this mechanism."""
        return len(self.abbreviations)


_MECHANISMS = {
    "host-genome-integrating": ViralMechanism(
        name="host-genome-integrating",
        title="Viruses that integrate into the host genome",
        description=(
            "Host-genome integration is part of viral persistence or is documented in a relevant "
            "tumour setting"),
        abbreviations=tuple(sorted(
            _IARC_GROUP_1_HPV | {"HBV", "HIV-1", "HIV-2", "HTLV-1", "MCV"})),
        source_urls=(
            "https://pubmed.ncbi.nlm.nih.gov/38624210/",
            "https://pubmed.ncbi.nlm.nih.gov/33563643/",
            "https://pubmed.ncbi.nlm.nih.gov/6328324/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2740911/",
        ),
        limitations=(
            "Integration is not universal across every tumour, and retroviral integration alone does "
            "not establish direct tumour-cell causality."
        ),
    ),
    "direct-tumor-cell-driver": ViralMechanism(
        name="direct-tumor-cell-driver",
        title="Direct tumour-cell viral drivers",
        description=(
            "IARC Group 1 viruses with a direct oncogenic role in at least one relevant tumour setting"),
        abbreviations=tuple(sorted(
            _IARC_GROUP_1_HPV | {"EBV (HHV-4)", "HBV", "HHV-8, KSHV", "HTLV-1", "MCV"})),
        source_urls=(
            "https://publications.iarc.who.int/119",
            "https://www.iarc.who.int/news-events/iarc-monographs-evaluation-of-the-carcinogenicity-of-hepatitis-d-virus-human-cytomegalovirus-and-merkel-cell-polyomavirus/",
        ),
        limitations=(
            "The relationship is specific to particular tumour types and does not imply that every "
            "infected cell or every tumour with a matching read is virus-driven."
        ),
    ),
    "viral-oncogene-dependence-evidence": ViralMechanism(
        name="viral-oncogene-dependence-evidence",
        title="Tumour-cell dependence on viral oncogene expression",
        description=(
            "Direct experimental evidence supports continued viral-oncogene expression in named "
            "virus-positive tumour settings"),
        abbreviations=("HPV16", "HPV18", "MCV"),
        source_urls=(
            "https://pubmed.ncbi.nlm.nih.gov/22056390/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2898224/",
        ),
        limitations=(
            "The evidence is subtype- and tumour-specific. It does not predict whether antiviral, "
            "immune, or cell-directed treatment will eradicate a patient's cancer."
        ),
    ),
}

_ALIASES = {
    "integrating": "host-genome-integrating",
    "integration": "host-genome-integrating",
    "direct-driver": "direct-tumor-cell-driver",
    "oncogene-dependence": "viral-oncogene-dependence-evidence",
}


def viral_mechanisms():
    """List source-backed, tumour-scoped viral mechanism selectors."""
    return tuple(_MECHANISMS[name] for name in sorted(_MECHANISMS))


def normalize_viral_mechanisms(mechanisms):
    """Normalize mechanism names and their concise aliases."""
    if mechanisms is None:
        return ()
    mechanisms = (mechanisms,) if isinstance(mechanisms, str) else tuple(mechanisms)
    if not mechanisms or any(not isinstance(value, str) or not value.strip() for value in mechanisms):
        raise ValueError("mechanisms must contain one or more nonempty names")
    normalized = tuple(dict.fromkeys(
        _ALIASES.get(value.strip().lower(), value.strip().lower()) for value in mechanisms))
    unknown = sorted(set(normalized) - set(_MECHANISMS))
    if unknown:
        raise ValueError(
            f"Unknown viral mechanism(s) {', '.join(unknown)}; available: "
            f"{', '.join(sorted(_MECHANISMS))}"
        )
    return normalized


def viral_mechanism_abbreviations(mechanisms, *, match="any"):
    """Return GDC abbreviations selected by mechanism names.

    Repeated mechanisms select a union by default and an intersection when
    ``match='all'``. The latter keeps only agents with every requested property.
    """
    names = normalize_viral_mechanisms(mechanisms)
    if match not in ("any", "all"):
        raise ValueError("match must be 'any' or 'all'")
    if not names:
        return ()
    groups = [set(_MECHANISMS[name].abbreviations) for name in names]
    selected = set.union(*groups) if match == "any" else set.intersection(*groups)
    return tuple(sorted(selected))
