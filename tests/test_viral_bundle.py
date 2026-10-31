"""Viral reference files keep genome, transcriptome and proteome layers separate."""

import gzip
import json
import shutil
import threading
import time

import pytest

from orfdata import viral_panel as viral_panel_module
from orfdata.data_resources import data_path
from orfdata.datasets import get_data_dir
from orfdata.downloads import collection_state
from orfdata.gtf import transcript_blocks
from orfdata.viral_annotation import split_genbank
from orfdata.viral_panel import (
    VIRAL_LAYERS,
    ViralReferenceFiles,
    bundle_version,
    download_viral_references,
    open_viral_references,
    prune_viral_references,
    select_viral_references,
    viral_reference_summary,
)
from orfdata.viral_registry import GDC_VIRUSES

CONTEXT = "cancer-transcriptome-screening"


def _install(tmp_path):
    return download_viral_references(species="human", cache_root=str(tmp_path))


def test_bundle_separates_the_three_layers(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    assert isinstance(bundle, ViralReferenceFiles)
    assert (bundle.genome_count, bundle.transcript_count, bundle.protein_count) == (4, 4, 4)
    assert [layer.name for layer in bundle.layers] == ["genome", "transcriptome", "proteome"]
    assert bundle.layer("genome").molecule == "dna"
    assert bundle.layer("transcriptome").annotation == bundle.annotation_gtf
    assert bundle.layer("proteome").records == 4  # One record annotates no CDS at all.
    with pytest.raises(ValueError, match="Unknown viral layer"):
        bundle.layer("exome")

    genomes = bundle.genome_fasta.read_text().splitlines()
    assert genomes[0].startswith(">VTEST001.1 abbreviation=TV1 taxon=10001 molecule=genome")
    assert bundle.to_dict()["layers"][1]["fasta"] == str(bundle.transcript_fasta)


def test_transcripts_and_proteins_share_identifiers(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    transcripts = {line.split()[0][1:]: line for line in
                   bundle.transcript_fasta.read_text().splitlines() if line.startswith(">")}
    proteins = {line.split()[0][1:] for line in
                bundle.protein_fasta.read_text().splitlines() if line.startswith(">")}
    assert set(transcripts) == {"VTEST001.1_cds1", "VTEST002.1_cds1", "VTEST002.1_cds2",
                                "VTEST003.1_cds1"}
    # A GenBank protein accession names a protein only when it names exactly one.
    assert proteins == {"VTEST001.1_cds1.p", "VTEST002.1_cds1.p", "VTEST002.1_cds2.p",
                        "AAA00003.1"}

    rows = {row["transcript_id"]: row for row in
            [dict(zip(*(line.split("\t") for line in
                        [bundle.transcript_table.read_text().splitlines()[0], body])))
             for body in bundle.transcript_table.read_text().splitlines()[1:]]}
    assert rows["VTEST001.1_cds1"]["annotation_complete"] == "yes"
    assert rows["VTEST002.1_cds1"]["codon_start"] == "2"
    assert rows["VTEST002.1_cds1"]["annotation_complete"] == "no"  # 59 coding bases, not whole codons
    assert rows["VTEST002.1_cds1"]["exon_count"] == "2"
    assert rows["VTEST003.1_cds1"]["translation_table"] == "11"
    assert rows["VTEST001.1_cds1"]["source_protein_id"] == "AAA00001.1"


def test_bundle_keeps_the_source_genbank_text_byte_for_byte(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    stored = sorted(bundle.record_directory.glob("*.gz"))
    # Every document is named for the entry, the source that served it and its identifier.
    assert [path.name for path in stored] == [
        "TV1__genbank__VTEST001.1.gb.gz", "TV2__genbank__VTEST002.1.gb.gz",
        "TV3__genbank__VTEST003.1.gb.gz", "TV4__genbank__VTEST004.1.gb.gz"]
    text = gzip.decompress(stored[0].read_bytes()).decode("utf-8")
    assert text.startswith("LOCUS") and text.rstrip().endswith("//")
    # Each stored record is what EFetch returns for that accession alone, blank
    # separator line and all, so splitting one back out is the identity.
    assert split_genbank(text.encode("utf-8")) == (text,)
    batch = b"".join(gzip.decompress(path.read_bytes()) for path in stored)
    assert split_genbank(batch) == tuple(
        gzip.decompress(path.read_bytes()).decode("utf-8") for path in stored)


def test_annotation_places_transcripts_on_the_viral_genomes(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    rows = [line.split("\t") for line in bundle.annotation_gtf.read_text().splitlines()
            if not line.startswith("#")]
    assert {row[0] for row in rows} == {"VTEST001.1", "VTEST002.1", "VTEST003.1"}
    assert {row[1] for row in rows} == {"orfdata"}

    first = [row for row in rows if 'transcript_id "VTEST001.1_cds1"' in row[8]]
    kinds = [row[2] for row in first]
    assert kinds == ["transcript", "exon", "CDS", "start_codon", "stop_codon"]
    transcript, exon, cds, start, stop = first
    assert (transcript[3], transcript[4], transcript[6]) == ("1", "30", "+")
    assert (cds[3], cds[4], cds[7]) == ("1", "27", "0")  # CDS stops before the stop codon.
    assert (start[3], start[4]) == ("1", "3") and (stop[3], stop[4]) == ("28", "30")
    assert 'transl_table "1";' in transcript[8]

    # The spliced, minus-strand, 5'-partial CDS: exons in transcript order, no codon rows.
    second = [row for row in rows if 'transcript_id "VTEST002.1_cds1"' in row[8]]
    exons = [(row[3], row[4]) for row in second if row[2] == "exon"]
    assert exons == [("61", "90"), ("11", "40")]
    assert all(row[6] == "-" for row in second)
    assert [row[2] for row in second].count("start_codon") == 0
    coding = [(row[3], row[4], row[7]) for row in second if row[2] == "CDS"]
    # codon_start=2 drops the first transcript base, which is the highest coordinate;
    # 29 coding bases then precede the second block, so one of its bases finishes a codon.
    assert coding == [("61", "89", "0"), ("11", "40", "1")]
    assert 'transl_table "11";' in [row[8] for row in rows
                                    if 'transcript_id "VTEST003.1_cds1"' in row[8]][0]

    genes = [row for row in rows if row[2] == "gene"]
    assert len(genes) == 4  # E6 on two strands of VTEST002 stays two genes.


def test_empty_trailing_table_values_stay_quoted(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    for table in (bundle.genome_table, bundle.transcript_table):
        assert not any(line.endswith("\t") for line in table.read_text().splitlines())
    assert bundle.genome_table.read_text().splitlines()[1].endswith('\t""')


def test_reuse_is_offline_and_tampering_is_caught(viral_panel, tmp_path):
    bundle = _install(tmp_path)
    reopened = open_viral_references(species="human", cache_root=str(tmp_path))
    assert reopened == bundle

    bundle.protein_fasta.write_text(bundle.protein_fasta.read_text() + ">EXTRA\nMM\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        open_viral_references(species="human", cache_root=str(tmp_path))
    with pytest.raises(FileNotFoundError, match="rebuild them with force=True"):
        download_viral_references(species="human", cache_root=str(tmp_path))


def test_missing_bundle_and_unknown_collection_explain_themselves(viral_panel, tmp_path):
    with pytest.raises(FileNotFoundError, match="are not installed"):
        open_viral_references(species="human", cache_root=str(tmp_path))
    with pytest.raises(ValueError, match="Unknown viral collection"):
        open_viral_references(species="human", collection="gdc", cache_root=str(tmp_path))
    with pytest.raises(ValueError, match="is a different target"):
        open_viral_references(species="mouse", cache_root=str(tmp_path))


def test_index_selects_whole_genomes_with_their_transcripts(viral_panel, tmp_path):
    _install(tmp_path)
    index = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert isinstance(index, ViralReferenceFiles)
    assert (index.genome_count, index.transcript_count, index.protein_count) == (2, 3, 3)
    assert index.subsets == () and index.contexts == (CONTEXT,)
    genomes = {line.split()[0][1:] for line in index.genome_fasta.read_text().splitlines()
               if line.startswith(">")}
    assert genomes == {"VTEST001.1", "VTEST002.1"}
    assert {line.split("\t")[0] for line in index.annotation_gtf.read_text().splitlines()
            if not line.startswith("#")} == genomes
    assert not any(line.endswith("\t") for line in index.genome_table.read_text().splitlines())

    assert select_viral_references(contexts=[CONTEXT], species="human",
                                  cache_root=str(tmp_path)).path == index.path
    with pytest.raises(ValueError, match="No viral references match"):
        select_viral_references(contexts=[CONTEXT], exclude_contexts=[CONTEXT],
                               species="human", cache_root=str(tmp_path))


def test_a_damaged_index_is_rebuilt_from_the_local_bundle(viral_panel, tmp_path):
    _install(tmp_path)
    index = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    index.protein_fasta.write_text("")
    rebuilt = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    # The rebuilt selection is published beside the damaged one, never over it.
    assert rebuilt.path != index.path and rebuilt.protein_count == 3
    assert (rebuilt.genome_count, rebuilt.contexts) == (index.genome_count, index.contexts)
    assert select_viral_references(contexts=[CONTEXT], species="human",
                                  cache_root=str(tmp_path)) == rebuilt
    assert json.loads((rebuilt.path / "metadata.json").read_text())["genome_ids"] == [
        "VTEST001.1", "VTEST002.1"]


def test_transcript_blocks_map_intervals_onto_both_strands():
    plus = ((11, 20, "+"), (31, 40, "+"))
    assert transcript_blocks(plus, 0, 20) == ((11, 20, "+", 0), (31, 40, "+", 10))
    assert transcript_blocks(plus, 5, 10) == ((16, 20, "+", 0), (31, 35, "+", 5))
    minus = ((31, 40, "-"), (11, 20, "-"))
    assert transcript_blocks(minus, 0, 20) == ((31, 40, "-", 0), (11, 20, "-", 10))
    assert transcript_blocks(minus, 5, 10) == ((31, 35, "-", 0), (16, 20, "-", 5))
    assert transcript_blocks(plus, 0, 0) == ()


def test_viral_genome_opens_for_downstream_tools(viral_panel, tmp_path):
    pytest.importorskip("pyensembl")
    bundle = _install(tmp_path)
    genome = bundle.pyensembl_genome(cache_root=str(tmp_path))
    # The fourth record annotates no CDS, so the annotation places nothing on it.
    assert sorted(genome.contigs()) == ["VTEST001.1", "VTEST002.1", "VTEST003.1"]

    transcript = genome.transcript_by_id("VTEST001.1_cds1")
    assert transcript.contig == "VTEST001.1" and transcript.strand == "+"
    assert transcript.complete
    assert transcript.sequence == "ATG" + "GCT" * 8 + "TAA"
    assert str(transcript.coding_sequence) == "ATG" + "GCT" * 8 + "TAA"
    assert transcript.protein_sequence == "MAAAAAAAA"
    assert genome.transcript_by_id("VTEST002.1_cds1").strand == "-"
    assert [exon.start for exon in genome.transcript_by_id("VTEST002.1_cds1").exons] == [61, 11]


def test_viral_genome_annotates_variants_for_downstream_tools(viral_panel, tmp_path):
    pytest.importorskip("varcode")
    from varcode import Variant

    genome = _install(tmp_path).pyensembl_genome(cache_root=str(tmp_path))
    # VTEST001.1 codes ATG GCT... from base 1 on the plus strand.
    effect = Variant(contig="VTEST001.1", start=4, ref="G", alt="A",
                     genome=genome).effects().top_priority_effect()
    assert effect.short_description == "p.A2T"
    assert effect.transcript_id == "VTEST001.1_cds1"


def test_summary_reads_the_installed_layers_without_writing(viral_panel, tmp_path):
    _install(tmp_path)
    whole = viral_reference_summary(species="human", cache_root=str(tmp_path))
    assert whole["collection"] == "gdc-viruses" and whole["selection"]["selected"] is False
    assert (whole["genome_count"], whole["transcript_count"]) == (4, 4)
    assert whole["complete_transcript_count"] == 3  # One CDS is 5'-partial.
    # A genome whose source record annotates no CDS is reported, not hidden.
    assert whole["genomes_without_transcripts"] == ["TV4"]
    assert [layer["name"] for layer in whole["layers"]] == list(VIRAL_LAYERS)

    selected = viral_reference_summary(subsets=None, contexts=[CONTEXT], species="human",
                                       cache_root=str(tmp_path))
    assert selected["genome_count"] == 2 and selected["layers"][0]["records"] == 2
    assert selected["selection"]["contexts"] == [CONTEXT]
    # Summarizing never materializes a derived index.
    assert not (get_data_dir(str(tmp_path)) / "viral" / "selections").exists()


def _slowed(monkeypatch, seconds=0.2):
    """Make every collection build take a while, and record how many run at once."""
    write, running, counts = viral_panel_module._write_bundle, [0], []
    guard = threading.Lock()

    def slow_write(*args, **kwargs):
        with guard:
            running[0] += 1
            counts.append(running[0])
        try:
            time.sleep(seconds)
            return write(*args, **kwargs)
        finally:
            with guard:
                running[0] -= 1

    monkeypatch.setattr(viral_panel_module, "_write_bundle", slow_write)
    return counts


def _store_path(tmp_path):
    return (get_data_dir(str(tmp_path)) / "viral" / "collections" / "homo_sapiens" /
            "gdc-viruses" / bundle_version(GDC_VIRUSES))


def test_a_refresh_never_removes_the_files_a_reader_holds(viral_panel, tmp_path, monkeypatch):
    old = _install(tmp_path)
    gtf = old.annotation_gtf.read_text()
    _slowed(monkeypatch)
    seen, done = [], threading.Event()

    def read_while_refreshing():
        while not done.is_set():
            try:
                seen.append(open_viral_references(species="human", cache_root=str(tmp_path)).path)
            except Exception as error:  # Any failure is a gap a reader would have hit.
                seen.append(error)

    reader = threading.Thread(target=read_while_refreshing)
    with old.annotation_gtf.open() as held:
        reader.start()
        try:
            new = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
        finally:
            done.set()
            reader.join()
        assert held.read() == gtf
    assert seen and all(path in (old.path, new.path) for path in seen), seen
    assert new.path != old.path and open_viral_references(species="human", cache_root=str(tmp_path)) == new
    # Nothing is pruned automatically: a program that resolved the old bundle
    # can still open it by path after one refresh, and after another.
    newer = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    for bundle in (old, new):
        assert bundle.annotation_gtf.read_text() == gtf
        assert viral_panel_module._verified_existing_bundle(
            bundle.path, species="homo_sapiens", collection=GDC_VIRUSES) == bundle
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == newer


def test_concurrent_refreshes_take_turns_and_pruning_is_explicit(viral_panel, tmp_path, monkeypatch):
    first = _install(tmp_path)
    counts = _slowed(monkeypatch)
    results, errors = [], []

    def refresh():
        try:
            results.append(download_viral_references(
                species="human", cache_root=str(tmp_path), force=True))
        except Exception as error:
            errors.append(error)

    threads = [threading.Thread(target=refresh) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors and len(results) == 2 and max(counts) == 1
    newest = max(results, key=lambda result: result.path.name)
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == newest
    bundles = newest.path.parent
    assert sorted(p.name for p in bundles.iterdir()) == sorted(
        bundle.path.name for bundle in (first, *results))

    selection = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    rebuilt = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path),
                                      force=True)
    index = viral_panel_module._annotation_index(str(tmp_path), GDC_VIRUSES, first.path)
    index.mkdir(parents=True)
    with pytest.raises(ValueError, match="keep must be a positive integer"):
        prune_viral_references(species="human", cache_root=str(tmp_path), keep=0)
    pruned = prune_viral_references(species="human", cache_root=str(tmp_path))
    assert set(pruned["deleted"]) == {first.path, *(r.path for r in results if r != newest), selection.path}
    assert pruned["skipped"] == {}
    assert [p.name for p in bundles.iterdir()] == [newest.path.name] and not index.exists()
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == newest
    assert select_viral_references(contexts=[CONTEXT], species="human",
                                  cache_root=str(tmp_path)) == rebuilt
    assert prune_viral_references(species="human", cache_root=str(tmp_path)) == dict(deleted=(), skipped={})


def test_a_collection_from_orfdata_044_opens_and_is_replaced_safely(
        viral_panel, tmp_path, monkeypatch):
    # ORFdata 0.44 kept a collection as one plain directory.
    installed = _install(tmp_path)
    legacy = (get_data_dir(str(tmp_path)) / "viral" / "references" / "homo_sapiens" /
              "gdc-viruses" / bundle_version(GDC_VIRUSES))
    legacy.parent.mkdir(parents=True)
    store = installed.path.parent.parent
    shutil.move(installed.path, legacy)
    shutil.rmtree(store)

    opened = open_viral_references(species="human", cache_root=str(tmp_path))
    assert opened.path == legacy
    assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "available"
    assert data_path("gdc-viruses", cache_root=str(tmp_path)) == legacy
    with monkeypatch.context() as offline:
        offline.setattr(viral_panel_module, "_write_bundle",
                        lambda *a, **k: pytest.fail("a verified 0.44 collection is reused"))
        assert download_viral_references(species="human", cache_root=str(tmp_path)) == opened

    with opened.annotation_gtf.open() as held:
        refreshed = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
        assert held.read() and legacy.is_dir()  # Kept as the previous collection.
    assert refreshed.path.parent == store / "bundles"
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == refreshed
    assert data_path("gdc-viruses", cache_root=str(tmp_path)) == refreshed.path

    again = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert legacy.is_dir() and refreshed.path.is_dir() and again.path.is_dir()
    assert set(prune_viral_references(species="human", cache_root=str(tmp_path))["deleted"]) == {
        legacy, refreshed.path}
    assert not legacy.exists() and again.path.is_dir()


def test_publishing_hashes_the_new_bundle_once(viral_panel, tmp_path, monkeypatch):
    _install(tmp_path)
    calls = []
    verify = viral_panel_module._verify_checksums
    monkeypatch.setattr(viral_panel_module, "_verify_checksums",
                        lambda path: calls.append(path) or verify(path))
    refreshed = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert len(calls) == 1 and calls[0].name == "bundle"  # The staged copy, before publication.
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == refreshed


def test_a_store_that_lost_its_marker_is_repaired_by_force(viral_panel, tmp_path):
    installed = _install(tmp_path)
    marker = _store_path(tmp_path) / ".datacache-store.json"
    marker.unlink()
    assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "corrupt"
    for attempt in (open_viral_references, download_viral_references):
        with pytest.raises(ValueError, match=r"lost its .datacache-store.json marker; "
                                             r"download_viral_references\(force=True\), or "
                                             r"'orfdata data fetch gdc-viruses --force', repairs it"):
            attempt(species="human", cache_root=str(tmp_path))
    repaired = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert marker.is_file() and repaired.path != installed.path
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == repaired
    assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "available"


def test_a_broken_newest_bundle_is_replaced_by_force(viral_panel, tmp_path):
    installed = _install(tmp_path)
    broken = installed.path.parent / "2999-01-01T00-00-00Z"
    broken.symlink_to(tmp_path)
    with pytest.raises(ValueError, match=r"isn't a directory; download_viral_references\(force=True\)"):
        open_viral_references(species="human", cache_root=str(tmp_path))
    repaired = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert repaired.path.name > broken.name
    assert open_viral_references(species="human", cache_root=str(tmp_path)) == repaired


def test_store_damage_force_cannot_fix_names_the_path_to_remove(viral_panel, tmp_path):
    installed = _install(tmp_path)
    bundles = installed.path.parent
    shutil.rmtree(bundles)
    bundles.write_text("not a directory\n")
    for force in (False, True):
        with pytest.raises(ValueError, match=f"remove {bundles} and install again"):
            download_viral_references(species="human", cache_root=str(tmp_path), force=force)
    assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "corrupt"


def test_a_damaged_selection_store_is_rebuilt(viral_panel, tmp_path):
    _install(tmp_path)
    selection = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    (selection.path.parent.parent / ".datacache-store.json").unlink()
    rebuilt = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert rebuilt.path != selection.path and rebuilt.genome_count == selection.genome_count


def test_an_interrupted_install_leaves_nothing_behind_for_long(viral_panel, tmp_path, monkeypatch):
    import os

    _install(tmp_path)
    store = _store_path(tmp_path)
    stale = store / f".staging-{os.getuid()}-dead"
    (stale / "bundle").mkdir(parents=True)
    stale.chmod(0o700)
    monkeypatch.setattr(viral_panel_module, "_write_bundle",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("interrupted")))
    with pytest.raises(RuntimeError, match="interrupted"):
        download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert not stale.exists()
    assert not [entry for entry in store.iterdir() if entry.name.startswith(".staging-")]


def test_an_unreadable_044_directory_does_not_block_an_install(viral_panel, tmp_path):
    import os

    if os.geteuid() == 0:
        pytest.skip("root reads directories without read permission")
    legacy_root = get_data_dir(str(tmp_path)) / "viral" / "references"
    legacy_root.mkdir(parents=True)
    legacy_root.chmod(0)
    try:
        installed = _install(tmp_path)
        assert installed.path.parent == _store_path(tmp_path) / "bundles"
        assert collection_state("gdc-viruses", cache_root=str(tmp_path)) == "available"
    finally:
        legacy_root.chmod(0o755)


def test_a_moved_datacache_store_module_is_reported(monkeypatch):
    from datacache import bundle_store

    monkeypatch.delattr(bundle_store, "user_key")
    with pytest.raises(ImportError, match="need datacache.bundle_store"):
        viral_panel_module._bundle_store()


def test_prune_skips_a_damaged_store_untouched_and_goes_on(viral_panel, tmp_path):
    first = _install(tmp_path)
    selection = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path), force=True)
    newest = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    selection_store = selection.path.parent.parent
    (selection_store / ".datacache-store.json").unlink()
    pruned = prune_viral_references(species="human", cache_root=str(tmp_path))
    assert pruned["deleted"] == (first.path,) and newest.path.is_dir()
    assert selection.path.is_dir()  # Nothing is deleted from a store that doesn't validate.
    assert list(pruned["skipped"]) == [str(selection_store)]
    reason = pruned["skipped"][str(selection_store)]
    assert "lost its .datacache-store.json marker" in reason and "force" not in reason


def test_a_selection_follows_a_refreshed_collection(viral_panel, tmp_path):
    import json

    installed = _install(tmp_path)
    selection = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert selection.source_collection == installed.path
    assert select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path)) == selection
    refreshed = download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    again = select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path))
    assert again.path != selection.path and again.source_collection == refreshed.path
    metadata = json.loads((again.path / "metadata.json").read_text())
    assert metadata["source_bundle"] == str(refreshed.path)
    # The replaced selection keeps its own provenance, not the collection's newest.
    assert json.loads((selection.path / "metadata.json").read_text())["source_bundle"] == str(installed.path)


def test_a_damaged_collection_seen_through_a_selection_names_what_repairs_it(viral_panel, tmp_path):
    _install(tmp_path)
    (_store_path(tmp_path) / ".datacache-store.json").unlink()
    for force in (False, True):
        with pytest.raises(ValueError, match=r"download_viral_references\(force=True\), or "
                                             r"'orfdata data fetch gdc-viruses --force', repairs it"):
            select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path), force=force)
    download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    assert select_viral_references(contexts=[CONTEXT], species="human", cache_root=str(tmp_path)).genome_count


def test_a_linked_044_directory_still_opens(viral_panel, tmp_path):
    installed = _install(tmp_path)
    moved = tmp_path / "elsewhere" / "collection"
    moved.parent.mkdir()
    store = installed.path.parent.parent
    shutil.move(installed.path, moved)
    shutil.rmtree(store)
    legacy = (get_data_dir(str(tmp_path)) / "viral" / "references" / "homo_sapiens" /
              "gdc-viruses" / bundle_version(GDC_VIRUSES))
    legacy.parent.mkdir(parents=True)
    legacy.symlink_to(moved, target_is_directory=True)
    opened = open_viral_references(species="human", cache_root=str(tmp_path))
    assert opened.path == legacy and collection_state("gdc-viruses", cache_root=str(tmp_path)) == "available"
    download_viral_references(species="human", cache_root=str(tmp_path), force=True)
    pruned = prune_viral_references(species="human", cache_root=str(tmp_path))
    assert pruned["deleted"] == () and "link to a 0.44 directory" in pruned["skipped"][str(legacy)]
    assert moved.is_dir()
