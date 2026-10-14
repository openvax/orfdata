"""Lifting catalogue ORFs between assemblies with chain files, and placing across annotation releases (#72)."""

import gzip
import json
from dataclasses import replace

import pytest
from test_proteome import genome, make_catalog, orf  # noqa: F401

from orfdata import ORFCatalog, ReferenceContext, extra_orfs, load_catalog
from orfdata.liftover import ChainFile, lift_catalog

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
    assert lifted.metadata["liftover"]["unmapped"] == [
        dict(orf_id="beyond", reason="A coding block crosses an alignment gap or edge")]
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
    with pytest.raises(ValueError, match=r"homo_sapiens GRCh37.*lift_catalog\(load_catalog\(species='homo_sapiens', "
                                         r"genome_build='GRCh38'\), ChainFile"):
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
    assert lifted.metadata["liftover"]["unmapped"] == [dict(orf_id="far", reason="A coding block is not aligned")]
    with pytest.raises(ValueError, match="keeps its original"):
        replace(record, protein_sequence="MKPA")
