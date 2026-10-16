"""Yang et al.'s RibORF ORFs for human, zebrafish, worm and yeast, and Chang et al.'s human ncORFs, from real rows.

See data/riborf/README.md and data/ncorf_mammals/README.md for the fixture rows.
"""

import gzip
import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import CatalogRegistry, RiboSeq, catalog_status, extra_orfs, load_catalog
from orfdata.reference import canonical_species

RIBORF, CHANG = Path(__file__).parent / "data/riborf", Path(__file__).parent / "data/ncorf_mammals"
FILES = {
    "riborf-zebrafish": dict(table=RIBORF / "SD3.xlsx", gtf=RIBORF / "zebrafish-ensembl115.gtf.gz",
                             cdna=RIBORF / "zebrafish-ensembl115.cdna.fa.gz",
                             ncrna=RIBORF / "zebrafish-ensembl115.ncrna.fa.gz"),
    "riborf-worm": dict(table=RIBORF / "SD3.xlsx", gtf=RIBORF / "worm-ensembl105.gtf.gz",
                        cdna=RIBORF / "worm-ensembl105.cdna.fa.gz", ncrna=RIBORF / "worm-ensembl105.ncrna.fa.gz"),
    "riborf-human": dict(table=RIBORF / "SD3.xlsx", gtf=RIBORF / "human-ensembl105.gtf.gz",
                         cdna=RIBORF / "human-ensembl105.cdna.fa.gz", ncrna=RIBORF / "human-ensembl105.ncrna.fa.gz"),
    "riborf-yeast": dict(table=RIBORF / "SD3.xlsx", gtf=RIBORF / "yeast-ensembl115.gtf.gz",
                         genome=RIBORF / "yeast-ensembl115.dna.fa.gz"),
    "ncorf-mammals-human": dict(bed=CHANG / "human.bed", proteins=CHANG / "human.pep.fa", table=CHANG / "S3.xlsx",
                                gtf=CHANG / "human-ensembl107.gtf.gz"),
}
# Each sheet's rows: an annotated CDS (not imported), a short ORF and the rest; see the READMEs.
COVERAGE = {
    "riborf-zebrafish": dict(records=5, distinct_protein_sequences=5, evidence_summaries=5, source_rows=7,
                             short_rows=1, resolved_rows=5, unresolved_rows=1),
    "riborf-worm": dict(records=5, distinct_protein_sequences=5, evidence_summaries=5, source_rows=6, short_rows=1,
                        resolved_rows=5, unresolved_rows=0),
    "riborf-human": dict(records=5, distinct_protein_sequences=5, evidence_summaries=5, source_rows=7, short_rows=1,
                         resolved_rows=5, unresolved_rows=1),
    "riborf-yeast": dict(records=4, distinct_protein_sequences=4, evidence_summaries=4, source_rows=6, short_rows=1,
                         resolved_rows=4, unresolved_rows=1),
    "ncorf-mammals-human": dict(records=3, distinct_protein_sequences=3, evidence_summaries=6, source_rows=4,
                                resolved_rows=3, unresolved_rows=1),
}


def definitions(*names):
    """Fixture registry entries: the packaged manifests with the fixtures pinned in place of the sources."""
    packaged, found = CatalogRegistry(), {}
    for name in names:
        manifest = packaged.info(name)
        files = FILES[name]
        sources = tuple(replace(s, url=f"https://example.invalid/{files[s.asset_id].name}",
                                sha256=hashlib.sha256(files[s.asset_id].read_bytes()).hexdigest(),
                                size=files[s.asset_id].stat().st_size) for s in manifest.sources)
        manifest = replace(manifest, sources=sources, coverage=COVERAGE[name])
        release = dict(manifest=manifest.to_dict(), distribution=dict(access_terms="fixture", redistributable=False),
                       assets={})
        found[name] = dict(default_version=manifest.version, versions={manifest.version: release},
                           priority=packaged.priority(name), optional=packaged.optional(name))
    return found


def built(tmp_path, *names):
    registry = CatalogRegistry(definitions(*names), cache_root=tmp_path / "cache")
    for name in names:
        registry.build(name, source_paths={k: str(v) for k, v in FILES[name].items()})
    return registry


def genome(tmp_path, build, release, species, gtf, fastas):
    from pyensembl import Genome
    from pyensembl.species import find_species_by_name

    paths = []
    for path in (gtf, *fastas):
        paths.append(tmp_path / path.name.removesuffix(".gz"))
        paths[-1].write_bytes(gzip.decompress(path.read_bytes()))
    found = Genome(build, "ensembl", str(release), gtf_path_or_url=str(paths[0]),
                   transcript_fasta_paths_or_urls=[str(p) for p in paths[1:]],
                   cache_directory_path=str(tmp_path / "pyensembl"))
    found.species = find_species_by_name(species)
    found.index()
    return found


def test_zebrafish_and_worm_are_known_species():
    assert canonical_species("zebrafish") == "danio_rerio"
    assert canonical_species("worm") == canonical_species("c_elegans") == "caenorhabditis_elegans"


@pytest.mark.parametrize("name, species, release", [("riborf-zebrafish", "zebrafish", 115),
                                                   ("riborf-worm", "caenorhabditis_elegans", 105)])
def test_new_species_map_through_their_transcripts_and_place_on_their_genome(tmp_path, name, species, release):
    catalog = built(tmp_path, name).open(name)
    assert all(r.orf_id.startswith("riborf_") and r.protein_sequence.startswith("M") for r in catalog.records)
    assert catalog.query(evidence=RiboSeq()) == catalog.records
    summary, = catalog.records[0].evidence.summaries
    samples = {"riborf-zebrafish": "zebrafish Ribo-seq data of Supplementary Data 1 (20 samples, 2 GEO series)",
               "riborf-worm": "worm Ribo-seq data of Supplementary Data 1 (39 samples, 1 GEO series)"}[name]
    assert samples in summary.definition
    files = FILES[name]
    loaded = genome(tmp_path, catalog.genome_build, release, species, files["gtf"], (files["cdna"], files["ncrna"]))
    try:
        assert extra_orfs(loaded, None, catalog=catalog).excluded == ()
    finally:
        loaded.close()


def test_zebrafish_leaves_out_an_orf_with_ambiguous_bases(tmp_path):
    catalog = built(tmp_path, "riborf-zebrafish").open("riborf-zebrafish")
    unresolved, = catalog.metadata["unresolved"]
    assert unresolved["reason"] == "The transcript's sequence has ambiguous bases in the ORF"
    # CTG is a start of the standard code, read as methionine without a declaration.
    truncation = next(r for r in catalog.records if r.transcript_id.startswith("ENSDART00000185589"))
    assert truncation.start_codon is None and truncation.protein_sequence.startswith("MEAPET")


def test_worm_names_with_dots_of_their_own_are_kept_whole(tmp_path):
    catalog = built(tmp_path, "riborf-worm").open("riborf-worm")
    assert {r.transcript_id for r in catalog.records} == {"F45H11.11.1", "F08G2.13", "T07C12.13.1", "C45B2.8.1",
                                                          "M163.15"}
    gtg = next(r for r in catalog.records if r.transcript_id == "M163.15")
    assert (gtg.start_codon, gtg.biotype) == ("GTG", "lncORF")


def test_yeast_orfs_come_from_the_genome_and_place_inside_ensembl_transcripts(tmp_path):
    catalog = built(tmp_path, "riborf-yeast").open("riborf-yeast")
    found = {r.transcript_id: r for r in catalog.records}
    # Only an ORF on Ensembl's (CDS-only) transcript takes its name; the others keep the source's transcript name,
    # with the Ensembl gene where there is one, and an unstable transcript Ensembl lacks keeps its own.
    assert set(found) == {"YNL046W_mRNA", "YNR004W", "YNL130C", "CUT808"}
    assert (found["YNR004W"].biotype, found["YNR004W"].gene_name) == ("uORF", "SWM2")
    assert (found["CUT808"].biotype, found["CUT808"].gene_id) == ("lncORF", "CUT808")
    # An extension across CPT1's intron is rebuilt from its two block sizes.
    extension = found["YNL130C"]
    assert (extension.biotype, extension.gene_name, len(extension.coding_blocks)) == ("extension", "CPT1", 2)
    assert len(extension.protein_sequence) == 407 and extension.protein_sequence.startswith("MDSKTV")
    # Ty1's frameshift joins its blocks across 1 nt, not an intron.
    unresolved, = catalog.metadata["unresolved"]
    assert (unresolved["transcript"], unresolved["reason"]) == (
        "YNL284C-B", "The source's two blocks are joined across 1 nt that aren't an intron (a programmed frameshift?)")
    loaded = genome(tmp_path, "R64-1-1", 115, "saccharomyces_cerevisiae", FILES["riborf-yeast"]["gtf"],
                    (RIBORF / "yeast-ensembl115.cdna.fa.gz",))
    try:
        placed = extra_orfs(loaded, None, catalog=catalog)
        assert {o.orf_id for o in placed} == {found["YNL046W_mRNA"].orf_id}
        assert {e.orf_id: e.reason for e in placed.excluded} == {
            found[t].orf_id: "Transcript absent from this annotation" for t in ("YNR004W", "YNL130C", "CUT808")}
        # No record names an Ensembl transcript that doesn't carry it.
        for transcript in loaded.transcripts():
            assert [o.orf_id for o in catalog.for_transcript(transcript)] == (
                [found["YNL046W_mRNA"].orf_id] if transcript.id == "YNL046W_mRNA" else [])
    finally:
        loaded.close()


def test_human_ribo_seq_catalogues_build_from_their_rows(tmp_path):
    registry = built(tmp_path, "riborf-human", "ncorf-mammals-human")
    riborf, chang = registry.open("riborf-human"), registry.open("ncorf-mammals-human")
    assert {r.biotype for r in riborf.records} == {"lncORF", "uORF", "iORF", "truncation"}
    unresolved, = riborf.metadata["unresolved"]
    assert unresolved["reason"] == "Transcript has no sequence in Ensembl 105's cDNA or ncRNA FASTA"
    assert sorted(r.start_codon or "" for r in chang.records) == ["", "", "GTG"]
    unresolved, = chang.metadata["unresolved"]
    assert unresolved["reason"] == "Coding length is not the protein plus a stop codon"


def test_human_ribo_seq_catalogues_are_used_only_when_asked_for(tmp_path):
    registry = built(tmp_path, "riborf-human", "ncorf-mammals-human")
    bundled = load_catalog()
    assert load_catalog(species="human", registry=registry).records == bundled.records
    everything = load_catalog(species="human", registry=registry, catalogs="all")
    assert len(everything.records) == len(bundled.records) + 5 + 3
    only = load_catalog(species="human", registry=registry, catalogs=["riborf-human"])
    assert {r.source for r in only.records} == {"riborf-human"}
    with pytest.raises(ValueError, match="No homo_sapiens catalogue named 'riborf-mouse'"):
        load_catalog(species="human", registry=registry, catalogs=["riborf-mouse"])
    with pytest.raises(ValueError, match="empty selection"):
        load_catalog(species="human", registry=registry, catalogs=[])
    assert load_catalog(species="human", registry=registry, catalogs=["all"]).records == everything.records
    with pytest.raises(ValueError, match="orfdata setup --species homo_sapiens --catalogs all"):
        load_catalog(species="human", catalogs="all", registry=CatalogRegistry(
            definitions("riborf-human", "ncorf-mammals-human"), cache_root=tmp_path / "empty"))
    with pytest.raises(ValueError, match="No ORF catalogue covers canis_lupus_familiaris"):
        load_catalog(species="dog", registry=registry, catalogs=["riborf-human"])
    # Chang et al.'s IDs name ORFs both human Ribo-seq catalogues have, as for mouse.
    assert registry.priority("ncorf-mammals-human") < registry.priority("riborf-human")
    status = {s.name: s for s in catalog_status(registry=registry)}
    assert status["riborf-human"].optional and not status["riborf-human"].selected
    assert status["bundled"].selected and not status["bundled"].optional
    # Everywhere else, every catalogue of the species is in the default set.
    packaged = CatalogRegistry()
    assert [n for n in ("riborf-zebrafish", "riborf-worm", "riborf-yeast", "riborf-mouse") if packaged.optional(n)] == []
