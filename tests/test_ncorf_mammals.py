"""Mouse Ribo-seq ORFs of Chang et al. from real source rows, bound to real Ensembl 107 transcripts."""

import gzip
import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import (
    MS,
    CatalogRegistry,
    RiboSeq,
    catalog_status,
    extra_orfs,
    install_species,
    load_catalog,
)

FIXTURES = Path(__file__).parent / "data/ncorf_mammals"
FILES = dict(bed="mouse.bed", proteins="mouse.pep.fa", table="S3.xlsx", gtf="ensembl107.gtf.gz")
# Six real rows: ATG, CTG, GTG and TTG starts, single- and multi-exon, and the one without a stop codon.
COVERAGE = dict(records=5, distinct_protein_sequences=5, evidence_summaries=10, source_rows=6, resolved_rows=5,
                unresolved_rows=1)


def definitions():
    packaged = CatalogRegistry().info("ncorf-mammals-mouse")
    sources = tuple(replace(s, url=f"https://example.invalid/{FILES[s.asset_id]}",
                            sha256=hashlib.sha256((FIXTURES / FILES[s.asset_id]).read_bytes()).hexdigest(),
                            size=(FIXTURES / FILES[s.asset_id]).stat().st_size) for s in packaged.sources)
    manifest = replace(packaged, sources=sources, coverage=COVERAGE)
    release = dict(manifest=manifest.to_dict(), distribution=dict(access_terms="fixture", redistributable=False),
                   assets={})
    return {manifest.catalog_id: dict(default_version=manifest.version, versions={manifest.version: release},
                                      priority=CatalogRegistry().priority(manifest.catalog_id))}


PATHS = {asset: str(FIXTURES / name) for asset, name in FILES.items()}


@pytest.fixture
def registry(tmp_path):
    registry = CatalogRegistry(definitions(), cache_root=tmp_path / "cache")
    registry.build("ncorf-mammals-mouse", source_paths=PATHS)
    return registry


@pytest.fixture
def ensembl107(tmp_path):
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    gtf, fasta = tmp_path / "e107.gtf", tmp_path / "e107.fa"
    gtf.write_bytes(gzip.decompress((FIXTURES / "ensembl107.gtf.gz").read_bytes()))
    fasta.write_bytes(gzip.decompress((FIXTURES / "ensembl107.cdna.fa.gz").read_bytes()))
    genome = Genome("GRCm39", "ensembl", "107", gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    cache_directory_path=str(tmp_path / "pyensembl"))
    genome.species = find_species_by_name("mouse")
    genome.index()
    yield genome
    genome.close()


def test_coding_paths_come_from_the_transcript_and_near_cognate_starts_read_as_methionine(registry):
    catalog = registry.open("ncorf-mammals-mouse")
    starts = {r.orf_id: r.start_codon for r in catalog.records}
    # GTG is declared; CTG and TTG are already starts of the standard code.
    assert starts == {"10_100323107_100323051_-": None, "10_10111062_10129911_+": None,
                      "10_100323141_100313880_-": "GTG", "10_128649914_128650436_+": None,
                      "10_22693539_22693429_-": None}
    gtg = catalog.get("10_100323141_100313880_-")
    assert gtg.protein_sequence.startswith("M") and len(gtg.coding_blocks) > 1  # The source writes V.
    assert catalog.get("10_10111062_10129911_+").gene_name  # From Ensembl, so --gene works by name.
    assert catalog.metadata["unresolved"] == [dict(
        orf_id="13_74704301_74702220_-", transcript=catalog.metadata["unresolved"][0]["transcript"],
        reason="Coding length is not the protein plus a stop codon")]


def test_library_counts_are_aggregate_ribo_seq_support(registry):
    catalog = registry.open("ncorf-mammals-mouse")
    libraries, callers = catalog.records[0].evidence.summaries
    assert (libraries.metric, callers.metric, libraries.scope, libraries.method) == (
        "Ribo-seq libraries", "ORF callers", "orf", "ribosome")
    assert len(catalog.query(evidence=RiboSeq())) == 5 and catalog.query(evidence=MS()) == ()
    assert catalog.query(evidence=RiboSeq(min_studies=2)) == ()  # Library counts name no studies.


def test_every_record_binds_to_its_ensembl_107_transcript(registry, ensembl107):
    catalog = registry.open("ncorf-mammals-mouse")
    for record in catalog.records:
        bound = record.bind(ensembl107.transcript_by_id(record.transcript_id.split(".")[0]))
        assert bound.protein_sequence == record.protein_sequence and bound.on_source_transcript
    found = extra_orfs(ensembl107, None, catalog=catalog)
    assert sorted({o.orf_id for o in found}) == sorted(r.orf_id for r in catalog.records) and found.excluded == ()


def test_every_registered_catalogue_for_a_species_is_merged_and_installed_together(tmp_path, monkeypatch):
    from test_openprot import fixture_registry

    openprot, openprot_paths = fixture_registry(tmp_path)
    both = CatalogRegistry({**openprot._definitions, **definitions()}, cache_root=tmp_path / "cache")
    builds = []

    def build(name, version=None, *, force=False, source_paths=None, download_options=None):
        builds.append(name)
        paths = openprot_paths if name == "openprot-mouse" else PATHS
        return CatalogRegistry.build(both, name, version, force=force, source_paths=paths)

    monkeypatch.setattr(both, "build", build)
    from orfdata import liftover

    monkeypatch.setattr(liftover, "download_chain", lambda *a, **k: "offline")  # No network in tests.
    result = install_species("mouse", registry=both)
    assert result["catalogs"] == {"ncorf-mammals-mouse 1": "built", "openprot-mouse 2.2-1": "built"}
    assert sorted(result["chains"]["skipped"]) == ["GRCm38_to_GRCm39", "GRCm39_to_GRCm38"]  # UCSC terms not accepted.
    assert install_species("mouse", registry=both)["catalogs"]["openprot-mouse 2.2-1"] == "already installed"
    merged = load_catalog(species="mouse", registry=both)
    # The registry's priority, not the names, orders the merge: OpenProt's IDs name shared ORFs.
    assert merged.source == "openprot-mouse:2.2-1+ncorf-mammals-mouse:1" and len(merged.records) == 5 + 8
    assert all(s.selected for s in catalog_status(registry=both) if s.species == "mus_musculus")
    # A damaged build is repaired, not only reported.
    built = next(p for p in (tmp_path / "cache").rglob("catalog.json.gz") if "ncorf-mammals-mouse" in str(p))
    built.chmod(0o644)
    built.write_bytes(built.read_bytes() + b"x")
    assert install_species("mouse", registry=both)["catalogs"]["ncorf-mammals-mouse 1"] == "repaired"
    assert len(load_catalog(species="mouse", registry=both).records) == 13
    with pytest.raises(ValueError, match="no ORF catalogue or chain file to install for canis_lupus_familiaris"):
        install_species("dog", registry=both)
