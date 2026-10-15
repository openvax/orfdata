"""Yang et al.'s mouse RibORF 2.0 ORFs from real source rows, placed through real Ensembl 107 transcripts."""

import gzip
import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import MS, CatalogRegistry, RiboSeq, extra_orfs

FIXTURES = Path(__file__).parent / "data/riborf"
FILES = dict(table="SD3.xlsx", gtf="ensembl107.gtf.gz", cdna="ensembl107.cdna.fa.gz", ncrna="ensembl107.ncrna.fa.gz")
# Nineteen real rows; see data/riborf/README.md. Canonical and readthrough rows are not imported, two are too short,
# and three are unresolved: a transcript that changed length, a TEC transcript without sequence and one whose
# sequence changed at the same length.
COVERAGE = dict(records=12, distinct_protein_sequences=12, evidence_summaries=12, source_rows=17, short_rows=2,
                resolved_rows=12, unresolved_rows=3)
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
    # ouORF, odORF and noncoding take Chang et al.'s names; ACG and GTG starts are declared.
    assert types == [("doORF", ""), ("extension", ""), ("iORF", ""), ("lncORF", ""), ("lncORF", "ACG"),
                     ("lncORF", "GTG"), ("truncation", ""), ("uORF", ""), ("uORF", ""), ("uORF", ""),
                     ("uORF", "ACG"), ("uoORF", "")]
    assert all(r.protein_sequence.startswith("M") and r.transcript_id.count(".") == 1 for r in catalog.records)
    assert min(len(r.protein_sequence) for r in catalog.records) == 5  # The 1- and 4-residue rows are left out.
    noncoding = next(r for r in catalog.records if r.transcript_id.startswith("ENSMUST00000235655"))
    assert len(noncoding.protein_sequence) == 26  # Translated from the ncRNA FASTA.
    reasons = {u["transcript"]: u["reason"] for u in catalog.metadata["unresolved"]}
    assert reasons == {
        "ENSMUST00000020909.3": "Transcript changed length since GENCODE vM25 (1773 nt, now 1839)",
        "ENSMUST00000192542.1": "Transcript has no sequence in Ensembl 107's cDNA or ncRNA FASTA",
        "ENSMUST00000170111.1": "The transcript's sequence doesn't start with the stated codon and end at its first stop"}
    # The same transcript's unchanged part still places.
    assert any(r.transcript_id.startswith("ENSMUST00000170111") for r in catalog.records)


def test_ids_never_take_a_name_chang_et_al_use(catalog):
    assert all(r.orf_id.startswith("riborf_") for r in catalog.records)
    shared = sorted(r.orf_id for r in catalog.records if "@" in r.orf_id)
    assert [i.split("@")[1] for i in shared] == ["ENSMUST00000226495.2", "ENSMUST00000228521.2"]
    # Two paths with one start and stop each answer to their own ID; the shared name names both.
    assert all(catalog.get(i).orf_id == i for i in shared)
    with pytest.raises(ValueError, match="names several ORFs"):
        catalog.get(shared[0].split("@")[0])


def test_presence_is_aggregate_ribo_seq_support(catalog):
    summary, = catalog.records[0].evidence.summaries
    assert (summary.metric, summary.value, summary.scope) == ("RibORF 2.0 translated ORF", 1, "orf")
    assert len(catalog.query(evidence=RiboSeq())) == 12 and catalog.query(evidence=MS()) == ()
    assert catalog.query(evidence=RiboSeq(min_studies=2)) == ()


def test_every_record_binds_to_its_ensembl_107_transcript(catalog, tmp_path):
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    gtf = tmp_path / "e107.gtf"
    gtf.write_bytes(gzip.decompress((FIXTURES / "ensembl107.gtf.gz").read_bytes()))
    fastas = []
    for kind in ("cdna", "ncrna"):
        fastas.append(tmp_path / f"e107.{kind}.fa")
        fastas[-1].write_bytes(gzip.decompress((FIXTURES / f"ensembl107.{kind}.fa.gz").read_bytes()))
    genome = Genome("GRCm39", "ensembl", "107", gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(f) for f in fastas],
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


def test_setup_downloads_the_ensembl_gtf_both_mouse_ribo_seq_catalogues_pin_once(monkeypatch, tmp_path):
    from orfdata import install

    fetched = []
    monkeypatch.setattr(install, "_fetch", lambda url, sha256, size, staging, options: fetched.append(url) or "local")
    registry = CatalogRegistry()
    manifests = [registry.info("ncorf-mammals-mouse"), registry.info("riborf-mouse"), registry.info("openprot-mouse")]
    shared = install._shared_sources(manifests, tmp_path, {})
    gtf, = (s for s in manifests[1].sources if s.asset_id == "gtf")
    assert shared == {gtf.sha256: "local"} and fetched == [gtf.url]
