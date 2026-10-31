"""OpenProt mouse ORFs from real source rows, bound to real Ensembl 106 transcripts (#18, #55)."""

import gzip
import hashlib
import io
import zipfile
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import MS, CatalogRegistry, RiboSeq, catalog_status, extra_orfs, load_catalog

FIXTURES = Path(__file__).parent / "data/openprot"
STEM = "mouse-openprot-2_2-altprots+isoforms-min_2_pep-GRCm39.106"
# Sixteen rows for eight accessions: 10 on Ensembl transcripts, 6 on RefSeq ones.
COVERAGE = dict(records=8, distinct_protein_sequences=7, evidence_summaries=14, source_rows=16, resolved_rows=10,
                unresolved_rows=6, accessions=8, resolved_accessions=7)


def fixture_registry(tmp_path, paths=None, **coverage):
    """The packaged release, with the fixture archives in place of OpenProt's."""
    paths = paths or {kind: FIXTURES / f"{STEM}.{kind}.zip" for kind in ("tsv", "bed", "fasta")}
    packaged = CatalogRegistry().info("openprot-mouse")
    sources = tuple(replace(s, url=f"https://example.invalid/{STEM}.{s.asset_id}.zip",
                            sha256=hashlib.sha256(paths[s.asset_id].read_bytes()).hexdigest(),
                            size=paths[s.asset_id].stat().st_size) for s in packaged.sources)
    manifest = replace(packaged, sources=sources, coverage=dict(COVERAGE, **coverage))
    release = dict(manifest=manifest.to_dict(), distribution=dict(access_terms="fixture", redistributable=False),
                   assets={})
    registry = CatalogRegistry({manifest.catalog_id: dict(default_version=manifest.version,
                                                          versions={manifest.version: release},
                                                          priority=CatalogRegistry().priority(manifest.catalog_id))},
                               cache_root=tmp_path / "cache")
    return registry, {k: str(v) for k, v in paths.items()}


@pytest.fixture
def catalog(tmp_path):
    registry, paths = fixture_registry(tmp_path)
    registry.build("openprot-mouse", source_paths=paths)
    return registry.open_build("openprot-mouse")


@pytest.fixture
def ensembl106(tmp_path):
    """Real Ensembl 106 GTF lines and cDNA/ncRNA sequences for the fixture's transcripts."""
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    gtf, fasta = tmp_path / "e106.gtf", tmp_path / "e106.fa"
    gtf.write_bytes(gzip.decompress((FIXTURES / "ensembl106.gtf.gz").read_bytes()))
    fasta.write_bytes(gzip.decompress((FIXTURES / "ensembl106.cdna.fa.gz").read_bytes()))
    genome = Genome("GRCm39", "ensembl", "106", gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path / "pyensembl"))
    genome.species = find_species_by_name("mouse")
    genome.index()
    yield genome
    genome.close()


def test_rows_become_records_from_their_bed_blocks(catalog):
    assert sorted(r.orf_id for r in catalog.records) == [
        "II_085990", "II_1000149", "II_912148", "IP_1064920", "IP_824214",
        "IP_824771@ENSMUST00000187264.2", "IP_824771@ENSMUST00000188314.2", "IP_825804"]
    rab1a = catalog.get("II_085990")
    # The TSV says 20151629-20174844; the BED blocks, and Ensembl, put the ORF one base later (#55).
    assert rab1a.coding_blocks == ((20151630, 20151652), (20165630, 20165702), (20168813, 20168908),
                                   (20174648, 20174845))
    assert (rab1a.transcript_id, rab1a.gene_id, rab1a.biotype, rab1a.strand) == (
        "ENSMUST00000109601.8", "ENSMUSG00000020149", "Isoform", "+")
    assert rab1a.reference_context.key == catalog.records[0].reference_context.key
    assert (rab1a.reference_context.annotation_name, rab1a.reference_context.annotation_version) == ("ensembl", "106")
    # One protein on transcripts of three overlapping genes is one record, on the first transcript.
    assert catalog.get("IP_1064920").transcript_id == "ENSMUST00000167713.2"
    # RefSeq rows are listed, not dropped silently; IP_957586 has only RefSeq rows.
    unresolved = catalog.metadata["unresolved"]
    assert len(unresolved) == 6 and {u["protein"] for u in unresolved} == {"IP_957586", "II_1000149", "IP_1064920"}
    assert catalog.metadata["unresolved_reasons"] == {
        "RefSeq transcript assignment; no RefSeq annotation is pinned": 6}
    with pytest.raises(KeyError):
        catalog.get("IP_957586")


def test_scores_are_aggregate_counts_for_the_protein(catalog):
    copies = [catalog.get(f"IP_824771@ENSMUST0000018{n}.2") for n in ("7264", "8314")]
    assert copies[0].protein_sequence == copies[1].protein_sequence and copies[0].contig == "Y"
    assert copies[0].evidence.summaries == copies[1].evidence.summaries  # One count, not two.
    ms, te = catalog.get("II_912148").evidence.summaries
    assert (ms.metric, ms.value, te.metric, te.value) == ("MS score", 3, "TE Translation Event score", 6)
    assert ms.scope == "protein_sequence" and ms.target.identifier == "II_912148"
    assert len(catalog.evidence_summaries) == 14

    def selected(rule):
        return sorted(r.orf_id.split("@")[0] for r in catalog.query(evidence=rule))

    assert len(selected(MS())) == 8
    # A sum of per-study unique peptides is not two distinct peptides, nor six studies.
    assert selected(MS(min_peptides=2)) == [] and selected(RiboSeq(min_studies=2)) == []
    assert selected(RiboSeq()) == selected(RiboSeq() & MS()) == ["II_912148", "IP_825804"]


def test_records_bind_to_their_native_ensembl_106_transcripts(catalog, ensembl106):
    for record in catalog.records:
        transcript = ensembl106.transcript_by_id(record.transcript_id.split(".")[0])
        if record.orf_id == "IP_825804":  # A TEC transcript: in Ensembl's GTF, absent from its cDNA/ncRNA files.
            with pytest.raises(ValueError, match="sequence unavailable"):
                record.bind(transcript)
            continue
        bound = record.bind(transcript)
        assert bound.protein_sequence == record.protein_sequence and bound.on_source_transcript
    found = extra_orfs(ensembl106, catalog=catalog, evidence=None)
    assert sorted({o.orf_id for o in found}) == sorted(r.orf_id for r in catalog.records if r.orf_id != "IP_825804")
    assert [(e.orf_id, e.reason) for e in found.excluded] == [("IP_825804", "ORF transcript sequence unavailable")]
    # IP_1064920's path is on three transcripts of three genes.
    assert sorted(o.transcript.id for o in found.for_orf("IP_1064920")) == [
        "ENSMUST00000167713", "ENSMUST00000232739", "ENSMUST00000233672"]
    assert len(found.for_orf("II_912148")) == 1 and len(catalog.get("II_912148").coding_blocks) == 15


def test_load_catalog_opens_the_build_and_status_reports_it(catalog, tmp_path):
    registry, paths = fixture_registry(tmp_path)  # Same cache root as the fixture's build.
    assert registry.published("openprot-mouse") is False
    mouse = load_catalog(species="mouse", registry=registry)
    assert mouse.records == catalog.records and mouse.source == "openprot-mouse:2.2-1"
    built = next((tmp_path / "cache").rglob("catalog.json.gz"))
    assert load_catalog(species="mouse", registry=registry,
                        expected_sha256=hashlib.sha256(built.read_bytes()).hexdigest()).records == catalog.records
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_catalog(species="mouse", registry=registry, expected_sha256="0" * 64)
    status = {s.name: s for s in catalog_status(registry=registry)}
    assert (status["openprot-mouse"].installed, status["openprot-mouse"].acquire) == ("built", "build")
    assert status["openprot-mouse"].selected and status["bundled"].installed == "bundled"


def rewrite(tmp_path, kind, change):
    path = tmp_path / f"{STEM}.{kind}.zip"
    with zipfile.ZipFile(FIXTURES / path.name) as source:
        text = change(source.read(f"{STEM}.{kind}").decode())
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as target:
        target.writestr(f"{STEM}.{kind}", text)
    path.write_bytes(buffer.getvalue())
    return path


def test_sources_that_disagree_fail_the_build(tmp_path):
    paths = {kind: FIXTURES / f"{STEM}.{kind}.zip" for kind in ("tsv", "bed", "fasta")}
    # A BED row moved one base, whole, no longer matches its TSV row's convention.
    moved = rewrite(tmp_path, "bed", lambda text: text.replace("chr11\t20151629\t20174845\t",
                                                               "chr11\t20151630\t20174846\t", 1))
    registry, sources = fixture_registry(tmp_path, dict(paths, bed=moved))
    with pytest.raises(ValueError, match=r"line 2 \(II_085990 on ENSMUST00000109601.8\) disagrees with its BED row"):
        registry.build("openprot-mouse", source_paths=sources)
    assert registry.inspect_build("openprot-mouse").status == "missing"
    mito = rewrite(tmp_path, "tsv", lambda text: text.replace("\t11\t20151629\t", "\tMT\t20151629\t", 1))
    registry, sources = fixture_registry(tmp_path, dict(paths, tsv=mito))
    with pytest.raises(ValueError, match="line 2 .II_085990. is mitochondrial"):
        registry.build("openprot-mouse", source_paths=sources)
    registry, sources = fixture_registry(tmp_path, source_rows=17)
    with pytest.raises(ValueError, match="source_rows is 16, but the manifest declares 17"):
        registry.build("openprot-mouse", source_paths=sources)


def test_the_packaged_mouse_release_is_built_from_openprot_never_downloaded():
    registry = CatalogRegistry()
    manifest = registry.info("openprot-mouse")
    tsv = next(s for s in manifest.sources if s.asset_id == "tsv")
    # The pin recorded when the source was first qualified (#18).
    assert tsv.sha256 == "4b4f466084224d4796a68d77913997eace60a31b43e1c9497d85f5af5e62c1b0"
    assert tsv.url.startswith("https://api.openprot.org/api/2.0/MM/downloads/") and not tsv.redistributable
    assert manifest.coverage["records"] == 2976 and manifest.coverage["unresolved_rows"] == 4571
    with pytest.raises(ValueError, match="no published artifact; build it"):
        registry.download("openprot-mouse")
    row = next(s for s in catalog_status() if s.name == "openprot-mouse")
    assert (row.species, row.genome_build, row.installed, row.acquire, row.selected) == (
        "mus_musculus", "GRCm39", None, "build", True)


def test_a_built_release_opens_inspects_and_queries_like_a_download(catalog, tmp_path, capsys):
    import json

    from orfdata.cli import run

    registry, _ = fixture_registry(tmp_path)
    assert registry.open("openprot-mouse").records == catalog.records
    assert registry.inspect("openprot-mouse").status == "available"
    definitions = tmp_path / "registry.json"
    definitions.write_text(json.dumps(registry._definitions))
    run(["catalogs", "--registry", str(definitions), "--cache-root", str(tmp_path / "cache"), "query",
         "openprot-mouse", "--orf-id", "II_085990", "--json"])
    assert [r["orf_id"] for r in json.loads(capsys.readouterr().out)["records"]] == ["II_085990"]
    # A split accession's bare ID names its records.
    with pytest.raises(ValueError, match=r"IP_824771 names several ORFs: IP_824771@ENSMUST00000187264.2, "
                                         r"IP_824771@ENSMUST00000188314.2"):
        catalog.get("IP_824771")
    with pytest.raises(ValueError, match="build it from its pinned sources .*'openprot-mouse', '2.2-1'"):
        registry.download("openprot-mouse")


def test_a_damaged_build_says_to_rebuild(catalog, tmp_path):
    registry, paths = fixture_registry(tmp_path)
    built = next((tmp_path / "cache").rglob("catalog.json.gz"))
    built.chmod(0o644)
    built.write_bytes(built.read_bytes() + b"x")
    assert next(s for s in catalog_status(registry=registry) if s.name == "openprot-mouse").installed == "damaged"
    for expected in (None, "0" * 64):
        with pytest.raises(ValueError, match=r"catalogue, which is damaged: build it again .*force=True.* --force\."):
            load_catalog(species="mouse", registry=registry, expected_sha256=expected)
    registry.build("openprot-mouse", source_paths=paths, force=True)
    assert load_catalog(species="mouse", registry=registry).records == catalog.records


def test_missing_build_requests_data_without_an_installation_extra(monkeypatch, tmp_path):
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    with pytest.raises(ValueError, match="need the openprot-mouse 2.2-1 catalogue, which is not "
                         "installed: build it.*the ncorf-mammals-mouse 1 catalogue.*orfdata setup --species mus_musculus"):
        load_catalog(species="mouse")
    row = next(s for s in catalog_status() if s.name == "openprot-mouse")
    assert (row.installed, row.acquire) == (None, "build")


def test_a_build_only_default_leaves_published_releases_working(tmp_path):
    from orfdata import CatalogRegistry

    registry, _ = fixture_registry(tmp_path)
    definitions = registry._definitions
    release = definitions["openprot-mouse"]["versions"]["2.2-1"]
    older = dict(release, manifest=dict(release["manifest"], version="2.2-0"),
                 assets={"catalog.json.gz": dict(url="https://example.invalid/c.json.gz", sha256="a" * 64, size=1)})
    definitions["openprot-mouse"]["versions"]["2.2-0"] = older
    registry = CatalogRegistry(definitions, cache_root=tmp_path / "cache")
    assert registry.inspect("openprot-mouse", "2.2-0").status == "missing"
    assert registry.inspect("openprot-mouse").status == "missing"  # The build-only default, not yet built.
    assert {s.version for s in catalog_status(registry=registry) if s.name == "openprot-mouse"} == {"2.2-0", "2.2-1"}
