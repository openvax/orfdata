"""Lifting catalogue ORFs between assemblies with chain files, and placing across annotation releases (#72)."""

import gzip
import json
from dataclasses import asdict, replace
from types import SimpleNamespace

import pytest
from test_proteome import genome, make_catalog, orf  # noqa: F401

from orfdata import ORFCatalog, ReferenceContext, extra_orfs, load_catalog
from orfdata.cli import run
from orfdata.liftover import (
    ChainFile,
    PinnedChain,
    download_chain,
    lift_catalog,
    load_chain,
    pinned_chains,
)

# GRCm38-like chr1 [1000, 1400) aligns ungapped to GRCm39-like 1 [0, 400); 2 [0, 90) aligns reversed onto 5.
CHAINS = """##matrix=axtChain 16 91,-114,-31,-123,-114,100,-125,-31,-31,-125,100,-114,-123,-31,-114,91
chain 100 chr1 10000 + 1000 1400 1 10000 + 0 400 1
400

chain 100 1 10000 + 1000 1400 1 10000 + 0 400 1
400

chain 50 chr2 5000 + 0 100 5 1000 - 900 1000 2
40 10 10
50

chain 10 chr3 5000 + 0 100 3 1000 + 0 100 3
100

chain 10 chr3 5000 + 50 100 HG1_PATCH 1000 + 0 50 4
50

chain 10 chr6 5000 + 0 100 7 1000 + 0 100 5
100

chain 10 chr6 5000 + 0 100 8 1000 + 0 100 6
100

chain 10 chr9 5000 + 0 100 NC_000009.12 1000 + 0 100 7
100
"""


@pytest.fixture
def chain(tmp_path):
    path = tmp_path / "GRCm38ToGRCm39.over.chain.gz"
    path.write_bytes(gzip.compress(CHAINS.encode()))
    chain = ChainFile(path, source_build="GRCm38", target_build="GRCm39")
    chain.path = path
    return chain


def chain_path(chain):
    return chain.path


def test_chain_blocks_map_whole_intervals_only(chain):
    assert chain.map("1", 1104, 1217) == ("1", "+", 104, 217)  # Repeated under chr1 and 1, one alignment.
    assert chain.map("chr1", 1001, 1400) == ("1", "+", 1, 400)
    assert chain.map("1", 1390, 1410) is None  # Runs past the aligned block.
    # Reversed alignment: source [0, 40) is reverse-strand [900, 940), forward [60, 100) of the 1,000 bases.
    assert chain.map("2", 1, 40) == ("5", "-", 61, 100)
    assert chain.map("2", 51, 60) == ("5", "-", 41, 50)
    assert chain.map("2", 35, 55) is None  # Spans the alignment gap.
    # chr3 also aligns to a fix patch: the same-named contig wins. chr6 aligns to two others: no answer.
    assert chain.map("3", 60, 70) == ("3", "+", 60, 70) and chain.map("3", 10, 20) == ("3", "+", 10, 20)
    assert chain.locate("6", 10, 20) == ([("7", "+", 10, 20), ("8", "+", 10, 20)], False)
    assert chain.map("6", 10, 20) is None and chain.locate("5", 10, 20) == ([], False)
    renamed = ChainFile(chain_path(chain), source_build="GRCm38", target_build="GRCm39",
                        contig_names={"NC_000009.12": "9"})
    assert renamed.map("9", 1, 10) == ("9", "+", 1, 10) and chain.map("9", 1, 10)[0] == "NC_000009.12"
    assert chain.liftover.source_build == "GRCm38" and len(chain.liftover.sha256) == 64


def shifted(record, offset=1000):
    context = record.initiation_context
    if context is not None:
        context = replace(context, genomic_positions=tuple(p + offset for p in context.genomic_positions))
    return replace(record, contig="chr1", genome_build="GRCm38", initiation_context=context,
                   coding_blocks=tuple((a + offset, b + offset) for a, b in record.coding_blocks))


def test_lifted_catalogue_places_like_the_native_one(chain, genome, tmp_path):  # noqa: F811
    native = make_catalog()
    beyond = replace(orf("beyond", blocks=((390, 410),), protein="MKPGKP"), contig="1")
    old = ORFCatalog([shifted(r) for r in (*native.records, beyond)], dict(native.metadata, genome_build="GRCm38"))
    lifted = lift_catalog(old, chain)
    assert lifted.genome_build == "GRCm39" and lifted.source == "orfdata-catalog:fixture-v1@GRCm39"
    assert [(u["orf_id"], u["reason"]) for u in lifted.metadata["liftover"]["unmapped"]] == [
        ("beyond", "A coding block crosses an alignment gap or edge")]
    assert lifted.metadata["liftover"]["input"]["metadata"]["genome_build"] == "GRCm38"
    for record in lifted.records:
        original = native.get(record.orf_id)
        assert (record.contig, record.coding_blocks) == ("1", original.coding_blocks)
        assert record.initiation_context == original.initiation_context
        assert record.lifted_from == old.get(record.orf_id) and record.liftover == chain.liftover

    def placements(catalog):
        return sorted((o.orf_id, o.transcript.id, o.start, o.end) for o in extra_orfs(genome, catalog=catalog))

    assert placements(lifted) == placements(native)
    path = tmp_path / "lifted.json"
    path.write_text(json.dumps(lifted.to_dict()))
    assert load_catalog(path).records == lifted.records
    with pytest.raises(ValueError, match="maps from GRCm38, not GRCm39"):
        lift_catalog(lifted, chain)


def test_a_reversed_alignment_flips_the_orf_strand(chain):
    record = replace(orf("rev", blocks=((1, 9), (51, 62)), protein="MKPGKP"), contig="2", genome_build="GRCm38")
    catalog = ORFCatalog([record], dict(make_catalog().metadata, genome_build="GRCm38"))
    lifted, = lift_catalog(catalog, chain).records
    assert (lifted.contig, lifted.strand, lifted.coding_blocks) == ("5", "-", ((39, 50), (92, 100)))


def test_placement_crosses_annotation_releases_but_binding_does_not(genome):  # noqa: F811
    older = ReferenceContext(10090, "mus_musculus", "GRCm39", "synthetic", "0")
    record = replace(make_catalog().get("ms-orf"), reference_context=older)
    catalog = ORFCatalog([record], make_catalog().metadata)
    with pytest.raises(ValueError, match="reference context mismatch"):
        record.bind(genome.transcript_by_id("test_rna"))
    found = extra_orfs(genome, catalog=catalog)
    assert sorted((o.transcript.id, o.on_source_transcript) for o in found) == [
        ("test_rna", True), ("test_rna_long", False)]
    other_strain = replace(older, strain="C3H")
    stranded = ORFCatalog([replace(record, reference_context=other_strain)], make_catalog().metadata)
    assert len(extra_orfs(genome, catalog=stranded)) == 0


def test_merged_catalogues_lift_with_their_members(chain):
    from orfdata import merge_catalogs

    native = make_catalog()
    first = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    ms = first.get("ms-orf")
    second = ORFCatalog([replace(ms, orf_id="IP_1", source="other")], dict(
        native.metadata, genome_build="GRCm38", dataset_version="other-v1", sources=[dict(id="other", scope="t")]))
    lifted = lift_catalog(merge_catalogs([first, second]), chain)
    record = lifted.get("IP_1")
    assert record.orf_id == "ms-orf" and [m.coding_blocks for m in record.members] == [
        native.get("ms-orf").coding_blocks] * 2
    assert all(m.lifted_from is not None for m in record.members)


def test_a_missing_chain_names_the_command_that_downloads_it():
    with pytest.raises(ValueError, match=r"homo_sapiens GRCh37 ORFs need the chain files GRCh38_to_GRCh37 to reach "
                                         r"GRCh37 by liftover: download them once with 'orfdata chains download "
                                         r"--species homo_sapiens' or download_chains\(species='homo_sapiens'\)"):
        load_catalog(species="human", genome_build="GRCh37")


def test_flanks_out_of_order_strains_and_unaligned_blocks(tmp_path):
    # The coding path lifts in place; the upstream flank lands elsewhere, out of order.
    path = tmp_path / "split.chain"
    path.write_text("chain 1 chr1 10000 + 102 300 1 10000 + 1102 1300 1\n198\n\n"
                    "chain 1 chr1 10000 + 90 102 1 10000 + 5000 5012 2\n12\n")
    chain = ChainFile(path, source_build="GRCm38", target_build="GRCm39")
    ms = make_catalog().get("ms-orf")
    strain = ReferenceContext(10090, "mus_musculus", "GRCm38", "synthetic", "1", strain="C3H")
    old = ORFCatalog([replace(ms, genome_build="GRCm38", reference_context=strain),
                      replace(orf("far", blocks=((4001, 4009),), protein="MK"), genome_build="GRCm38")],
                     dict(make_catalog().metadata, genome_build="GRCm38"))
    lifted = lift_catalog(old, chain)
    record, = lifted.records
    assert record.coding_blocks == ((1103, 1104), (1205, 1217))
    assert record.initiation_context.genomic_positions == () and record.initiation_context.upstream == "CCC"
    pin = record.reference_context
    assert (pin.strain, pin.genome_build, pin.annotation_name, pin.source_sha256) == (
        "C3H", "GRCm39", None, chain.liftover.sha256)
    assert [(u["orf_id"], u["reason"]) for u in lifted.metadata["liftover"]["unmapped"]] == [
        ("far", "A coding block is not aligned")]
    with pytest.raises(ValueError, match="keeps its original"):
        replace(record, protein_sequence="MKPA")


def test_orfs_that_do_not_lift_are_reported_where_they_would_be_placed(chain, genome):  # noqa: F811
    native = make_catalog()
    beyond = replace(orf("beyond", blocks=((390, 410),), protein="MKPGKP", ribosome_studies=("s",)), contig="1")
    old = ORFCatalog([shifted(r) for r in (*native.records, beyond)], dict(native.metadata, genome_build="GRCm38"))
    found = extra_orfs(genome, None, catalog=lift_catalog(old, chain))
    assert ("beyond", "Did not lift to GRCm39: A coding block crosses an alignment gap or edge") in [
        (e.orf_id, e.reason) for e in found.excluded]
    # Under a rule it doesn't meet, an unlifted ORF isn't reported.
    assert "beyond" not in [e.orf_id for e in extra_orfs(genome, "ms", catalog=lift_catalog(old, chain)).excluded]


def test_a_lifted_record_can_come_first_in_a_merge(chain):
    from orfdata import merge_catalogs

    native = make_catalog()
    old = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    lifted = lift_catalog(old, chain)
    other = ORFCatalog([replace(r, source="other", ribosome_studies=("extra",)) for r in native.records],
                       dict(native.metadata, dataset_version="other-v1", sources=[dict(id="other", scope="t")]))
    for catalogs in ((lifted, other), (other, lifted)):
        record = merge_catalogs(catalogs).get("ms-orf")
        assert record.liftover is None and "extra" in record.ribosome_studies
        assert [m.liftover is not None for m in record.members] == [c is lifted for c in catalogs]


def test_ucsc_scaffold_names_and_strands_in_chain_files(tmp_path):
    path = tmp_path / "scaffold.chain"
    path.write_text("chain 1 chr1_GL456210v1_random 1000 + 0 100 chrUn_JH584304v1 1000 + 0 100 1\n100\n")
    chain = ChainFile(path, source_build="GRCm38", target_build="GRCm39")
    assert chain.map("GL456210.1", 1, 10) == ("JH584304.1", "+", 1, 10)
    path.write_text("chain 1 chr1 1000 - 0 100 chr1 1000 + 0 100 1\n100\n")
    with pytest.raises(ValueError, match="source strand is -"):
        ChainFile(path, source_build="GRCm38", target_build="GRCm39")


def test_several_places_are_listed_apart_from_the_reason(chain):
    record = replace(orf("twice", blocks=((10, 18),), protein="MK"), contig="6", genome_build="GRCm38")
    lifted = lift_catalog(ORFCatalog([record], dict(make_catalog().metadata, genome_build="GRCm38")), chain)
    entry, = lifted.metadata["liftover"]["unmapped"]
    assert (entry["reason"], entry["contigs"]) == ("A coding block aligns to several places", ["7", "8"])


def test_pinned_chains_name_their_terms():
    pins = {c.name: c for c in pinned_chains()}
    assert set(pins) == {"GRCh37_to_GRCh38", "GRCh38_to_GRCh37", "GRCh38_to_NCBI36", "NCBI36_to_GRCh38",
                         "GRCh38_to_T2T-CHM13v2.0", "T2T-CHM13v2.0_to_GRCh38",
                         "GRCm38_to_GRCm39", "GRCm39_to_GRCm38", "GRCm38_to_NCBIM37", "NCBIM37_to_GRCm38"}
    assert pins["GRCh38_to_T2T-CHM13v2.0"].provider == "T2T Consortium" and not pins["GRCh38_to_T2T-CHM13v2.0"].restricted
    assert [c.name for c in pinned_chains(species="mouse")] == [
        "GRCm38_to_GRCm39", "GRCm38_to_NCBIM37", "GRCm39_to_GRCm38", "NCBIM37_to_GRCm38"]
    assert [c.provider for c in pinned_chains() if c.restricted] == ["UCSC", "UCSC"]
    assert not pins["GRCh38_to_GRCh37"].restricted and pins["GRCm39_to_GRCm38"].restricted
    assert pins["GRCm39_to_GRCm38"].filename == "mm39ToMm10.over.chain.gz"
    assert all(len(c.sha256) == 64 and c.size > 0 and c.url.startswith("https://") for c in pins.values())


@pytest.fixture
def pinned(chain, monkeypatch):
    from orfdata import liftover

    pin = PinnedChain("GRCm38_to_GRCm39", "mus_musculus", "GRCm38", "GRCm39",
                      "https://example.invalid/" + chain.path.name, chain.liftover.sha256,
                      chain.path.stat().st_size, "UCSC", "Non-commercial only", "https://example.invalid/README", True)
    entry = asdict(pin)
    pins = {entry.pop("name"): entry}
    monkeypatch.setattr(liftover, "_pins", lambda: pins)
    return pin


def test_a_pinned_chain_installs_verified_and_opens_offline(chain, pinned, tmp_path):
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match="not downloaded; download it with download_chain"):
        load_chain("GRCm38", "GRCm39", cache_root=cache)
    with pytest.raises(ValueError, match="UCSC with these terms: Non-commercial only.*accept_terms=True"):
        download_chain("GRCm38", "GRCm39", source_path=chain.path, cache_root=cache)
    with pytest.raises(ValueError, match="No pinned chain file maps GRCm39 to GRCm38 .pinned: GRCm38_to_GRCm39."):
        load_chain("GRCm39", "GRCm38", cache_root=cache)
    wrong = tmp_path / "wrong.chain.gz"
    wrong.write_bytes(b"x" * pinned.size)
    with pytest.raises(ValueError, match="mismatch"):
        download_chain("GRCm38", "GRCm39", source_path=wrong, accept_terms=True, cache_root=cache)
    path = download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True, cache_root=cache)
    assert path.read_bytes() == chain.path.read_bytes()
    # The same file gives the same provenance, downloaded or not.
    assert load_chain("GRCm38", "GRCm39", cache_root=cache).liftover == chain.liftover


def test_cli_lifts_with_a_chain_file_or_a_downloaded_pin(chain, pinned, tmp_path, capsys):
    native = make_catalog()
    old = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    source = tmp_path / "old.json"
    old.save(source)
    by_file, by_pin = tmp_path / "file.json.gz", tmp_path / "pin.json.gz"
    run(["lift", "--catalog", str(source), "--to", "GRCm39", "--chain", str(chain.path), "--output", str(by_file)])
    summary = json.loads(capsys.readouterr().err)
    assert (summary["records"], summary["lifted_records"], summary["unmapped"]) == (4, 4, 0)
    assert summary["lifted"] == "orfdata-catalog:fixture-v1@GRCm39"
    assert load_catalog(by_file).records == lift_catalog(old, chain).records
    cache = str(tmp_path / "cache")
    run(["chains", "--cache-root", cache, "download", "GRCm38", "GRCm39", "--source-path", str(chain.path),
         "--accept-terms"])
    assert json.loads(capsys.readouterr().out)["path"].endswith(chain.path.name)
    run(["lift", "--catalog", str(source), "--to", "GRCm39", "--cache-root", cache, "--output", str(by_pin)])
    assert by_pin.read_bytes() == by_file.read_bytes()
    for options, message in ((["--species", "mouse"], "not both"),
                             (["--to", "GRCm40"], "No pinned chain files lift mus_musculus GRCm38 to GRCm40")):
        argv = ["lift", "--catalog", str(source), "--to", "GRCm39", "--output", str(by_pin)]
        with pytest.raises(SystemExit):
            run(argv[:3] + argv[5:] + options if "--to" in options else argv + options)
        assert message in capsys.readouterr().err


def test_cli_lists_pins_without_downloading(capsys):
    run(["chains", "list", "--species", "human"])
    listed = json.loads(capsys.readouterr().out)
    assert [c["name"] for c in listed] == ["GRCh37_to_GRCh38", "GRCh38_to_GRCh37", "GRCh38_to_NCBI36",
                                           "GRCh38_to_T2T-CHM13v2.0", "NCBI36_to_GRCh38", "T2T-CHM13v2.0_to_GRCh38"]
    assert listed[0]["provider"] == "Ensembl" and listed[0]["restricted"] is False


def test_a_pinned_file_must_be_used_in_its_own_direction_and_species(chain, pinned):
    with pytest.raises(ValueError, match="is the pinned mus_musculus GRCm38 to GRCm39 chain, not GRCm39 to GRCm38"):
        ChainFile(chain.path, source_build="GRCm39", target_build="GRCm38")
    with pytest.raises(ValueError, match="not from GRCm38 to itself"):
        ChainFile(chain.path, source_build="GRCm38", target_build="GRCm38")
    pinned_file = ChainFile(chain.path, source_build="GRCm38", target_build="GRCm39")
    assert pinned_file.species == "mus_musculus" and chain.species is None  # The fixture was read before the pin.
    human = [replace(shifted(r), taxon_id=9606, species="homo_sapiens") for r in make_catalog().records]
    catalog = ORFCatalog(human, dict(make_catalog().metadata, genome_build="GRCm38", taxon_id=9606,
                                     species="homo_sapiens"))
    with pytest.raises(ValueError, match="is a mus_musculus chain; the catalogue is homo_sapiens"):
        lift_catalog(catalog, pinned_file)


def test_cli_rejects_a_lift_to_the_same_assembly_or_against_the_pin(chain, pinned, tmp_path, capsys):
    native = make_catalog()
    source = tmp_path / "old.json"
    ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38")).save(source)
    for to, message in (("GRCm38", "already on GRCm38"), ("GRCm40", "is the pinned mus_musculus GRCm38 to GRCm39")):
        with pytest.raises(SystemExit):
            run(["lift", "--catalog", str(source), "--to", to, "--chain", str(chain.path),
                 "--output", str(tmp_path / "out.json")])
        assert message in capsys.readouterr().err


def test_damaged_and_missing_chains_say_how_to_download_again(chain, pinned, tmp_path):
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match=rf"download_chain\('GRCm38', 'GRCm39', accept_terms=True, "
                                         rf"cache_root='{cache}'\) or 'orfdata chains --cache-root {cache} "
                                         r"download GRCm38 GRCm39 --accept-terms'"):
        load_chain("GRCm38", "GRCm39", cache_root=cache)
    path = download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True, cache_root=cache)
    path.chmod(0o644)
    path.write_bytes(b"x" * path.stat().st_size)
    for call in (lambda: load_chain("GRCm38", "GRCm39", cache_root=cache),
                 lambda: download_chain("GRCm38", "GRCm39", accept_terms=True, cache_root=cache)):
        with pytest.raises(ValueError, match="damaged.*force=True.*--force'"):
            call()
    download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True, force=True, cache_root=cache)
    assert load_chain("GRCm38", "GRCm39", cache_root=cache).liftover == chain.liftover


def test_an_assembly_without_a_pinned_chain_suggests_a_chain_file():
    with pytest.raises(ValueError, match=r"--to HG002v1.1 --chain FILE --output FILE', or lift_catalog\(.*"
                                         r"ChainFile\(path, source_build='GRCh38', target_build='HG002v1.1'\)"):
        load_catalog(species="human", genome_build="HG002v1.1")
    # T2T-CHM13 is reached through the T2T Consortium's public-domain chain.
    with pytest.raises(ValueError, match="need the chain files GRCh38_to_T2T-CHM13v2.0 to reach T2T-CHM13v2.0"):
        load_catalog(species="human", genome_build="T2T-CHM13v2.0")


def test_lifted_and_merged_catalogues_export_as_tsv(chain, tmp_path, capsys):
    import csv
    import io

    from orfdata import merge_catalogs

    native = make_catalog()
    old = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    other = ORFCatalog([replace(old.get("ms-orf"), orf_id="IP_1", source="other")], dict(
        old.metadata, dataset_version="other-v1", sources=[dict(id="other", scope="t")]))
    path = tmp_path / "lifted.json"
    lift_catalog(merge_catalogs([old, other]), chain).save(path)
    run(["--catalog", str(path), "--format", "tsv"])
    rows = list(csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t"))
    merged = next(r for r in rows if r["orf_id"] == "ms-orf")
    assert json.loads(merged["liftover_json"])["target_build"] == "GRCm39"
    assert [m["orf_id"] for m in json.loads(merged["members_json"])] == ["ms-orf", "IP_1"]
    assert json.loads(merged["lifted_from_json"])["genome_build"] == "GRCm38"


def test_a_genome_build_alone_picks_the_catalogue(capsys):
    with pytest.raises(SystemExit):
        run(["--genome-build", "GRCh37", "--gene", "TPST1"])
    assert "homo_sapiens GRCh37 ORFs need the chain files GRCh38_to_GRCh37" in capsys.readouterr().err


def test_placement_on_another_version_of_the_source_transcript_says_so(genome):  # noqa: F811
    from orfdata import Proteome

    older = ReferenceContext(10090, "mus_musculus", "GRCm39", "synthetic", "0")
    base = Proteome.from_genome(genome)
    for transcript_id, binding in (("test_rna.1", None), ("test_rna.2", "ORF transcript version mismatch")):
        record = replace(make_catalog().get("ms-orf"), transcript_id=transcript_id, reference_context=older)
        proteome = base.with_extra_orfs(None, catalog=ORFCatalog([record], make_catalog().metadata))
        entry, = [e for e in proteome if e.kind == "extra_orf"]
        assert entry.source_binding == binding
        assert any(o.on_source_transcript for o in entry.occurrences) is (binding is None)
        assert proteome.layers[-1]["not_on_source_transcript"] == (binding is not None)
        transcripts = [o.transcript_id for o in entry.occurrences]
        assert len(transcripts) == len(set(transcripts))  # Each transcript once.


def test_a_lifted_orf_is_on_its_source_transcript_under_the_new_assemblys_version(chain, genome):  # noqa: F811
    native = make_catalog()
    old = ORFCatalog([replace(shifted(native.get("ms-orf")), transcript_id="test_rna.7")],
                     dict(native.metadata, genome_build="GRCm38"))
    found = extra_orfs(genome, None, catalog=lift_catalog(old, chain))
    assert [o.on_source_transcript for o in found.for_orf("ms-orf") if o.transcript.id == "test_rna"] == [True]
    assert found.source_failures == ()


def test_a_bare_accession_that_names_several_records_is_never_guessed():
    from orfdata import merge_catalogs

    split = ORFCatalog([replace(orf("IP_9@test_rna"), source="b"),
                        replace(orf("IP_9@test_rna_long", transcript_id="test_rna_long", blocks=((218, 226),),
                                    protein="MK"), source="b")],
                       dict(make_catalog().metadata, dataset_version="b-v1", sources=[dict(id="b", scope="t")]))
    bare = ORFCatalog([replace(orf("IP_9", blocks=((206, 214),), protein="KP"), source="c")],
                      dict(make_catalog().metadata, dataset_version="c-v1", sources=[dict(id="c", scope="t")]))
    with pytest.raises(ValueError, match="IP_9 names several ORFs"):
        merge_catalogs([split, bare]).get("IP_9")


class Registry:
    """Installed catalogues by name, as load_catalog sees a CatalogRegistry."""

    def __init__(self, **catalogs):
        self.catalogs = catalogs

    def list(self):
        return tuple(SimpleNamespace(catalog_id=name, version="1", genome_build=c.genome_build,
                                     references=(SimpleNamespace(species=c.species),))
                     for name, c in self.catalogs.items())

    def info(self, name, version=None):
        return next(m for m in self.list() if m.catalog_id == name)

    def published(self, name, version=None):
        return True

    def inspect(self, name, version=None, verify_files=True):
        return SimpleNamespace(status="available")

    def open(self, name, version=None):
        return self.catalogs[name]


def placements(found):
    return sorted((o.orf_id, o.transcript.id, o.start, o.end) for o in found)


def test_load_catalog_lifts_catalogues_from_other_assemblies(chain, pinned, genome):  # noqa: F811
    native = make_catalog()
    old = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    with pytest.raises(ValueError, match=r"need the chain files GRCm38_to_GRCm39 to reach GRCm39 by liftover: "
                                         r"download them once with 'orfdata chains download --species mus_musculus "
                                         r"--accept-terms'"):
        load_catalog(species="mouse", genome_build="GRCm39", registry=Registry(old=old))
    download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True)
    lifted = load_catalog(species="mouse", genome_build="GRCm39", registry=Registry(old=old))
    assert lifted.source == "orfdata-catalog:fixture-v1@GRCm39" and lifted.records == lift_catalog(old, chain).records
    assert placements(extra_orfs(genome, catalog=lifted)) == placements(extra_orfs(genome, catalog=native))
    with pytest.raises(ValueError, match="expected_sha256 checks one catalogue file"):
        load_catalog(species="mouse", genome_build="GRCm39", registry=Registry(old=old), expected_sha256="0" * 64)
    # Another source's catalogue on this assembly too: both are used, merged, this assembly's first.
    other = ORFCatalog([replace(r, source="other") for r in native.records],
                       dict(native.metadata, dataset_version="other-v1", sources=[dict(id="other", scope="t")]))
    merged = load_catalog(species="mouse", genome_build="GRCm39", registry=Registry(old=old, new=other))
    assert merged.source == "orfdata-catalog:other-v1+orfdata-catalog:fixture-v1@GRCm39"
    record = merged.get("ms-orf")
    assert [m.liftover is not None for m in record.members] == [False, True]
    assert placements(extra_orfs(genome, catalog=merged)) == placements(extra_orfs(genome, catalog=native))


def test_routes_take_the_fewest_hops_and_lift_through_each(chain, tmp_path, monkeypatch, capsys):
    from orfdata import liftover

    identity = tmp_path / "GRCm39ToGRCm40.chain"
    identity.write_text("chain 1 1 10000 + 0 10000 1 10000 + 0 10000 9\n10000\n")
    pins = {}
    for name, path, source, target in (("GRCm38_to_GRCm39", chain.path, "GRCm38", "GRCm39"),
                                       ("GRCm39_to_GRCm40", identity, "GRCm39", "GRCm40")):
        pins[name] = dict(species="mus_musculus", source_build=source, target_build=target,
                          url=f"https://example.invalid/{path.name}", sha256=ChainFile(
                              path, source_build=source, target_build=target).liftover.sha256,
                          size=path.stat().st_size, provider="test", terms="none", terms_url="https://example.invalid",
                          restricted=False)
    monkeypatch.setattr(liftover, "_pins", lambda: pins)
    routes = liftover.chain_routes("mouse", "GRCm38")
    assert {b: [c.name for c in r] for b, r in routes.items()} == {
        "GRCm39": ["GRCm38_to_GRCm39"], "GRCm40": ["GRCm38_to_GRCm39", "GRCm39_to_GRCm40"]}
    assert liftover.chain_routes("mouse", "GRCm40") == {}
    for name, path, source, target in (("a", chain.path, "GRCm38", "GRCm39"), ("b", identity, "GRCm39", "GRCm40")):
        download_chain(source, target, source_path=path)
    native = make_catalog()
    old = ORFCatalog([shifted(r) for r in native.records], dict(native.metadata, genome_build="GRCm38"))
    far = liftover.lift_to(old, "GRCm40")
    assert far.source == "orfdata-catalog:fixture-v1@GRCm39@GRCm40"
    assert [r.coding_blocks for r in far.records] == [r.coding_blocks for r in native.records]
    assert far.get("ms-orf").lifted_from.lifted_from == old.get("ms-orf")
    source = tmp_path / "old.json"
    old.save(source)
    run(["lift", "--catalog", str(source), "--to", "GRCm40", "--output", str(tmp_path / "far.json")])
    summary = json.loads(capsys.readouterr().err)
    assert [c["chain"] for c in summary["chains"]] == [chain.path.name, identity.name]
    assert summary["lifted_records"] == 4 and load_catalog(tmp_path / "far.json").records == far.records
    # Lifting a lifted file reports only this command's hops.
    lift_catalog(old, chain).save(tmp_path / "mid.json")
    run(["lift", "--catalog", str(tmp_path / "mid.json"), "--to", "GRCm40", "--output", str(tmp_path / "far2.json")])
    assert [c["chain"] for c in json.loads(capsys.readouterr().err)["chains"]] == [identity.name]


def test_download_all_skips_restricted_chains_without_consent(pinned, capsys):
    from orfdata import download_chains

    result = download_chains()
    assert result["downloaded"] == {} and result["skipped"]["GRCm38_to_GRCm39"].startswith(
        "UCSC's terms need accept_terms=True")
    with pytest.raises(ValueError, match="No chain files are pinned for dog"):
        download_chains(species="dog")
    for argv, message in ((["chains", "download", "GRCm38"], "Name both assemblies"),
                          (["chains", "download", "--source-path", "x"], "--source-path installs one chain")):
        with pytest.raises(SystemExit):
            run(argv)
        assert message in capsys.readouterr().err


def test_without_datacache_a_lift_says_what_to_install(monkeypatch):
    monkeypatch.setitem(__import__("sys").modules, "datacache", None)
    with pytest.raises(ValueError, match="need orfdata.data. installed, then the chain files GRCh38_to_GRCh37"):
        load_catalog(species="human", genome_build="GRCh37")
    assert load_catalog(species="human").genome_build == "GRCh38"  # The bundled catalogue needs nothing.


def test_status_lists_the_assemblies_each_catalogue_reaches():
    from orfdata import catalog_status

    rows = {row.name: row for row in catalog_status()}
    assert rows["bundled"].lifts_to == (
        dict(genome_build="GRCh37", chains=["GRCh38_to_GRCh37"], downloaded=False),
        dict(genome_build="NCBI36", chains=["GRCh38_to_NCBI36"], downloaded=False),
        dict(genome_build="T2T-CHM13v2.0", chains=["GRCh38_to_T2T-CHM13v2.0"], downloaded=False))
    assert rows["openprot-mouse"].lifts_to[1] == dict(
        genome_build="NCBIM37", chains=["GRCm39_to_GRCm38", "GRCm38_to_NCBIM37"], downloaded=False)


def mouse(records, build, version, source="synthetic"):
    metadata = dict(make_catalog().metadata, genome_build=build, dataset_version=version,
                    sources=[dict(id=source, scope="t")])
    return ORFCatalog([replace(r, genome_build=build, source=source) for r in records], metadata)


def test_load_catalog_ignores_what_it_cannot_reach_and_counts_a_source_once(chain, pinned, tmp_path):
    native = make_catalog()
    old = mouse([shifted(r) for r in native.records], "GRCm38", "old-v1")
    far = [mouse([orf("F1")], "GRCm40", "f1-v1", "f1"), mouse([orf("F2")], "GRCm40", "f2-v1", "f2")]
    # Two catalogues on an assembly no chain reaches are never candidates, so not ambiguous.
    assert load_catalog(species="mouse", genome_build="GRCm38",
                        registry=Registry(old=old, f1=far[0], f2=far[1])).records == old.records
    # Two catalogues that would be lifted are ambiguous, and the error says so.
    shifted_records = [shifted(r) for r in native.records]
    twins = Registry(a=mouse(shifted_records, "GRCm38", "a-v1", "a"), b=mouse(shifted_records, "GRCm38", "b-v1", "b"))
    with pytest.raises(ValueError, match="Several ORF catalogues cover mus_musculus GRCm38 to lift to GRCm39 "
                                         r"\(a 1, b 1\), and none is chosen for you"):
        load_catalog(species="mouse", genome_build="GRCm39", registry=twins)
    # One source registered on two assemblies under two names counts once: from the nearer assembly.
    registry = Registry(op39=native, op38=old)
    registry._cache_root = tmp_path / "root"  # Chains are found where the registry's catalogues are.
    download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True, cache_root=tmp_path / "root")
    assert load_catalog(species="mouse", genome_build="GRCm39", registry=registry).records == native.records
    assert load_catalog(species="mouse", genome_build="GRCm38", registry=registry).records == old.records


def test_other_names_and_other_species_assemblies_are_explained():
    with pytest.raises(ValueError, match="orfdata names that assembly GRCh38, as Ensembl does; give a custom genome "
                                         "the reference name 'GRCh38'"):
        load_catalog(species="human", genome_build="hg38")
    with pytest.raises(ValueError, match="orfdata names that assembly GRCh38"):
        load_catalog(species="human", genome_build="GRCh38.p14")
    with pytest.raises(ValueError, match="GRCm39 is a mus_musculus assembly, not homo_sapiens"):
        load_catalog(species="human", genome_build="GRCm39")


def test_an_explicit_catalogue_on_another_assembly_is_lifted(chain, pinned, genome):  # noqa: F811
    native = make_catalog()
    old = mouse([shifted(r) for r in native.records], "GRCm38", "fixture-v1")
    with pytest.raises(ValueError, match="is on GRCm38; lifting it to GRCm39 needs the chain files GRCm38_to_GRCm39"):
        extra_orfs(genome, None, catalog=old)
    download_chain("GRCm38", "GRCm39", source_path=chain.path, accept_terms=True)
    assert placements(extra_orfs(genome, None, catalog=old)) == placements(extra_orfs(genome, None, catalog=native))
    other = mouse(native.records, "GRCm39", "other-v1", "other")
    merged = extra_orfs(genome, None, catalog=[old, other])
    assert merged.catalog.source == "orfdata-catalog:fixture-v1@GRCm39+orfdata-catalog:other-v1"
