"""Lifting catalogue ORFs between assemblies with chain files, and placing across annotation releases (#72)."""

import gzip
import json
from dataclasses import asdict, replace

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


def test_a_missing_assembly_names_the_catalogue_to_lift():
    with pytest.raises(ValueError, match=r"homo_sapiens GRCh37.*'orfdata lift --species homo_sapiens --genome-build "
                                         r"GRCh38 --to GRCh37 --output FILE', or lift_catalog\(load_catalog\("
                                         r"species='homo_sapiens', genome_build='GRCh38'\), load_chain\('GRCh38', 'GRCh37'\)\)"):
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
    assert set(pins) == {"GRCh37_to_GRCh38", "GRCh38_to_GRCh37", "GRCm38_to_GRCm39", "GRCm39_to_GRCm38"}
    assert [c.name for c in pinned_chains(species="mouse")] == ["GRCm38_to_GRCm39", "GRCm39_to_GRCm38"]
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
                             (["--chain", str(chain.path), "--cache-root", cache], "does not apply with --chain"),
                             (["--to", "GRCm40"], "No pinned chain file maps GRCm38 to GRCm40")):
        argv = ["lift", "--catalog", str(source), "--to", "GRCm39", "--output", str(by_pin)]
        with pytest.raises(SystemExit):
            run(argv[:3] + argv[5:] + options if "--to" in options else argv + options)
        assert message in capsys.readouterr().err


def test_cli_lists_pins_without_downloading(capsys):
    run(["chains", "list", "--species", "human"])
    listed = json.loads(capsys.readouterr().out)
    assert [c["name"] for c in listed] == ["GRCh37_to_GRCh38", "GRCh38_to_GRCh37"]
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
    with pytest.raises(ValueError, match=r"--to T2T-CHM13v2.0 --chain FILE --output FILE', or lift_catalog\(.*"
                                         r"ChainFile\(path, source_build='GRCh38', target_build='T2T-CHM13v2.0'\)"):
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
