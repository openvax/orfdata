"""Shared synthetic fixtures; these are not experimental evidence."""

import hashlib
import json

import pandas as pd
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
from orfdata import viral_panel as viral_panel_module
from orfdata.record_sources import GENBANK_SOURCE, reference_source
from orfdata.reference_ledger import PUBLISHED, empty_ledger, ledger_from_rows

VIRAL_CONTEXT = "cancer-transcriptome-screening"


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


class Unquoted(str):
    """A qualifier value GenBank writes without quotes, such as transl_except."""


def _genbank(accession, sequence, features, *, name, description="test genome", taxon=10000):
    """One GenBank flat file, laid out as NCBI EFetch returns it."""
    lines = [f"LOCUS       {name:<16}{len(sequence):>10} bp    DNA     linear   VRL 01-JAN-2020",
             f"DEFINITION  {description}.",
             f"ACCESSION   {accession.split('.')[0]}",
             f"VERSION     {accession}",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{len(sequence)}",
             f'                     /db_xref="taxon:{taxon}"']
    for location, qualifiers in features:
        lines.append(f"     CDS             {location}")
        for key, value in qualifiers.items():
            if value is None:  # A qualifier that is a flag, such as ribosomal_slippage.
                lines.append(f"                     /{key}")
            elif isinstance(value, Unquoted) or not isinstance(value, str):
                lines.append(f"                     /{key}={value}")
            else:
                lines.append(f'                     /{key}="{value}"')
    lines.append("ORIGIN")
    for offset in range(0, len(sequence), 60):
        block = sequence[offset:offset + 60].lower()
        chunks = " ".join(block[index:index + 10] for index in range(0, len(block), 10))
        lines.append(f"{offset + 1:>9} {chunks}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def _payload(records):
    """Serialize synthetic records the way an NCBI EFetch response carries them."""
    return "\n".join(records).encode("utf-8")


def viral_coverage_filter(inventory, *, contexts=None, exclude_contexts=None,
                          context_match="any", **ignored):
    """Select synthetic inventory rows the way viral_panel_coverage does."""
    keep = pd.Series(True, index=inventory.index)
    if contexts:
        match = any if context_match == "any" else all
        keep &= inventory["contexts"].map(
            lambda held: match(context in held for context in contexts))
    if exclude_contexts:
        keep &= ~inventory["contexts"].map(
            lambda held: any(context in held for context in exclude_contexts))
    return inventory[keep]


def ledger_entry(abbreviation, sequence, *, index=1, state=PUBLISHED, features=1, source="pave"):
    """One synthetic ledger row, measured off a synthetic molecule."""
    return {"index": index, "abbreviation": abbreviation, "source": source,
            "identifier": f"{abbreviation}REF", "state": state, "length_nt": len(sequence),
            "sequence_md5": hashlib.md5(sequence.encode()).hexdigest(),
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "coding_features": features, "taxon_id": 0, "title": ""}


def viral_documents():
    """Four synthetic source records standing in for the 200 published ones."""
    # An ORF with a real start and stop codon, so a translation round-trips.
    coding = "ATG" + "GCT" * 8 + "TAA"
    first = _genbank("VTEST001.1", coding + "ACGT" * 10, [
        ("1..30", {"gene": "E1", "product": "protein E1", "protein_id": "AAA00001.1",
                   "codon_start": 1, "translation": "MAAAAAAAA"}),
    ], name="VTEST001", description="Test virus 1", taxon=10001)
    # A spliced, minus-strand CDS that starts one base in: partial, so not complete.
    second = _genbank("VTEST002.1", "ACGT" * 25, [
        ("complement(join(11..40,61..90))",
         {"gene": "E6", "product": "protein E6", "codon_start": 2,
          "translation": "MCCCCCCCCCCCCCCCCCC"}),
        ("91..99", {"gene": "E6", "product": "protein E6 short", "protein_id": "AAA00001.1",
                    "codon_start": 1, "translation": "MD"}),
    ], name="VTEST002", description="Test virus 2", taxon=10002)
    # A non-standard genetic code the source declares.
    third = _genbank("VTEST003.1", "ATG" + "GCT" * 13 + "TAA" + "ACGT" * 5, [
        ("1..45", {"gene": "L1", "product": "protein L1", "protein_id": "AAA00003.1",
                   "codon_start": 1, "transl_table": 11, "translation": "MAAAAAAAAAAAAA"}),
    ], name="VTEST003", description="Test virus 3", taxon=10003)
    # Some published records annotate no CDS at all, as HPV2, HPV8 and HPV51 do,
    # even though the genome plainly carries a reading frame.
    fourth = _genbank("VTEST004.1", "ATG" + "GCT" * 5 + "TAA" + "ACGT" * 10, [],
                      name="VTEST004", description="Test virus 4", taxon=10004)
    records = {"VTEST001.1": first, "VTEST002.1": second, "VTEST003.1": third,
               "VTEST004.1": fourth}
    return records


@pytest.fixture
def viral_panel(monkeypatch):
    """The synthetic panel, served as a collection's primary source."""
    records = viral_documents()
    inventory = pd.DataFrame([
        {"index": index, "name": f"Test virus {index}", "abbreviation": abbreviation,
         "refseq": f"NC_00000{index}.1", "genbank": accession, "pave_id": "", "pave_notes": "",
         "bundled_reference_id": "", "coding_models": pd.NA, "status": "not_bundled",
         "available": False, "contexts": (VIRAL_CONTEXT,) if index < 3 else ()}
        for index, (abbreviation, accession) in enumerate(
            (("TV1", "VTEST001.1"), ("TV2", "VTEST002.1"), ("TV3", "VTEST003.1"),
             ("TV4", "VTEST004.1")), 1)])
    inventory["coding_models"] = inventory["coding_models"].astype("Int64")

    def coverage(**options):
        return viral_coverage_filter(inventory, **options)

    monkeypatch.setattr(viral_panel_module, "viral_panel_coverage", coverage)
    serve_records(monkeypatch, records)
    # The primary ledger is built from the served records the same way the
    # packaged one is, so a fixture pins what a real collection pins.
    ledger = reference_source(GENBANK_SOURCE).build_ledger(inventory, timeout=1)
    monkeypatch.setattr(viral_panel_module, "_primary_ledger",
                        lambda collection, rows: ledger.loc[rows["abbreviation"]])
    return inventory


def serve_alternate_ledger(monkeypatch, rows):
    """Answer for the alternate source's ledger, so the packaged one is never read."""
    ledger = ledger_from_rows(rows)
    collection = viral_panel_module.GDC_VIRUSES
    monkeypatch.setattr(type(collection), "ledger",
                        lambda self, source=None: ledger
                        if source in (None, self.alternate_source) else empty_ledger())
    return ledger


def serve_records(monkeypatch, genbank=(), alternate=()):
    """Answer a collection's source requests from synthetic records.

    Both sources read through one HTTP seam, so a test can serve the GenBank
    records a collection downloads and the alternate reference documents it
    asks for, and nothing else is faked.
    """
    from orfdata import record_sources
    from orfdata.net import SourceUnavailable

    genbank, alternate = dict(genbank), dict(alternate)
    asked = []

    def fetch(url, *, timeout=None, params=None, **options):
        asked.append(url)
        if params and "id" in params:
            return _payload([genbank[name] for name in params["id"].split(",")])
        identifier = url.rsplit("/", 1)[-1]
        if identifier in alternate:
            document = alternate[identifier]
            return document if isinstance(document, bytes) else json.dumps(document).encode()
        raise SourceUnavailable(f"{url} was refused: HTTP 404 Not Found")

    monkeypatch.setattr(record_sources, "fetch", fetch)
    return asked
