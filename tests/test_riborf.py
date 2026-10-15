"""Yang et al.'s mouse RibORF 2.0 ORFs from real source rows, placed through real Ensembl 107 transcripts."""

import gzip
import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import MS, CatalogRegistry, RiboSeq, extra_orfs

FIXTURES = Path(__file__).parent / "data/riborf"
FILES = dict(table="SD3.xlsx", gtf="ensembl107.gtf.gz", cdna="ensembl107.cdna.fa.gz", ncrna="ensembl107.ncrna.fa.gz")
# Eight real rows: canonical and readthrough are not imported; one transcript changed since GENCODE vM25.
COVERAGE = dict(records=5, distinct_protein_sequences=5, evidence_summaries=5, source_rows=6, resolved_rows=5,
                unresolved_rows=1)
PATHS = {asset: str(FIXTURES / name) for asset, name in FILES.items()}


@pytest.fixture
def catalog(tmp_path):
    packaged = CatalogRegistry().info("riborf-mouse")
    sources = tuple(replace(s, url=f"https://example.invalid/{FILES[s.asset_id]}",
                            sha256=hashlib.sha256((FIXTURES / FILES[s.asset_id]).read_bytes()).hexdigest(),
                            size=(FIXTURES / FILES[s.asset_id]).stat().st_size) for s in packaged.sources)
    manifest = replace(packaged, sources=sources, coverage=COVERAGE)
    release = dict(manifest=manifest.to_dict(), distribution=dict(access_terms="fixture", redistributable=False),
                   assets={})
    registry = CatalogRegistry({manifest.catalog_id: dict(default_version=manifest.version,
                                                          versions={manifest.version: release})},
                               cache_root=tmp_path / "cache")
    registry.build("riborf-mouse", source_paths=PATHS)
    return registry.open("riborf-mouse")


def test_rows_map_through_their_transcript_and_translate_there(catalog):
    types = sorted((r.biotype, r.start_codon or "") for r in catalog.records)
    assert types == [("extension", ""), ("noncoding", "GTG"), ("truncation", ""), ("uORF", ""), ("uORF", "ACG")]
    assert all(r.protein_sequence.startswith("M") and r.transcript_id.count(".") == 1 for r in catalog.records)
    unresolved, = catalog.metadata["unresolved"]
    assert unresolved["transcript"] == "ENSMUST00000020909.3"
    assert unresolved["reason"].startswith("Transcript changed length since GENCODE vM25")
    assert all(r.biotype != "canonical" for r in catalog.records)


def test_presence_is_aggregate_ribo_seq_support(catalog):
    summary, = catalog.records[0].evidence.summaries
    assert (summary.metric, summary.value, summary.scope) == ("RibORF 2.0 translated ORF", 1, "orf")
    assert len(catalog.query(evidence=RiboSeq())) == 5 and catalog.query(evidence=MS()) == ()


def test_every_record_binds_to_its_ensembl_107_transcript(catalog, tmp_path):
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    gtf, fasta = tmp_path / "e107.gtf", tmp_path / "e107.fa"
    gtf.write_bytes(gzip.decompress((FIXTURES / "ensembl107.gtf.gz").read_bytes()))
    fasta.write_bytes(gzip.decompress((FIXTURES / "ensembl107.cdna.fa.gz").read_bytes()))
    genome = Genome("GRCm39", "ensembl", "107", gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(fasta)],
                    cache_directory_path=str(tmp_path / "pyensembl"))
    genome.species = find_species_by_name("mouse")
    genome.index()
    try:
        for record in catalog.records:
            bound = record.bind(genome.transcript_by_id(record.transcript_id.split(".")[0]))
            assert bound.protein_sequence == record.protein_sequence and bound.on_source_transcript
        assert extra_orfs(genome, None, catalog=catalog).excluded == ()
    finally:
        genome.close()
