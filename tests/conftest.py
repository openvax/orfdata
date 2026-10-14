"""Shared synthetic quantitative fixtures; these are not experimental evidence."""

import hashlib

import pytest

from orfdata import (
    Assay,
    DatasetManifest,
    Measurement,
    ReferenceContext,
    Sample,
    ScopedTarget,
    SourceAsset,
)


@pytest.fixture(autouse=True)
def private_cache(tmp_path_factory, monkeypatch):
    """Tests never see, or write to, the user's own downloads and builds."""
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path_factory.mktemp("orfdata-cache")))


@pytest.fixture
def experiment(tmp_path):
    source = tmp_path / "input.tsv"
    source.write_text("source fixture\n")
    contexts = tuple(ReferenceContext(taxon, species, assembly, "ensembl", "110")
                     for taxon, species, assembly in (
                         (9606, "human", "GRCh38"), (10090, "mouse", "GRCm39"),
                         (9615, "dog", "ROS_Cfam_1.0"), (4932, "yeast", "R64-1-1")))
    manifest = DatasetManifest(
        "fixture", "1", "Synthetic four-species measurement validation", contexts,
        (SourceAsset("input.tsv", source.as_uri(), hashlib.sha256(source.read_bytes()).hexdigest(),
                     source.stat().st_size, "Synthetic test data; Apache-2.0", redistributable=True),),
        ("Synthetic test fixture, not experimental evidence",), ("rna_seq", "digest_ms", "ribo_seq"),
        dict(searched="Deliberately colliding test identifiers", exclusions="No biological inference"),
        dict(software="orfdata tests", version="1", parameters={} ),
    )
    samples = tuple(Sample("sample", context, biological_unit_id="animal",
                           attributes={"source_taxon": context.taxon_id}) for context in contexts)
    assays = tuple(Assay("rna", "sample", "rna_seq", context, technical_replicate="1")
                   for context in contexts)
    observations = tuple(Measurement(
        "value", ScopedTarget(context, "gene", "same_id", "ensembl:110"),
        "sample", "rna", "rna_abundance", index, "TPM", "fixture", "1",
        provenance=dict(source_asset="input.tsv", row=index + 1),
    ) for index, context in enumerate(contexts))
    return manifest, samples, assays, observations


