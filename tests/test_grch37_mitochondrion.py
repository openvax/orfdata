"""Under GRCh37, hg19's chrM is not the rCRS that Ensembl names MT (#200).

hg19's chrM is NC_001807 (16,571 bp); the rCRS (NC_012920, 16,569 bp) is
Ensembl's MT, and GENCODE's GRCh37 and hg38 call it chrM. Their positions are
offset: hg19 chrM:12995 is rCRS MT:12994. Other builds' chrM is the rCRS.
"""

from dataclasses import replace

from orfdata import ORFCatalog, merge_catalogs
from orfdata.catalog import ORFEvidence, ORFRecord
from orfdata.liftover import ChainFile
from orfdata.model import ORFModel

# Ensembl's GRCh37_to_GRCh38 chain repeats its MT chain as chrM, with the
# rCRS length. The third, hg19 chrM chain is synthetic: hg19 bases 12991-12992
# have no rCRS counterpart, so later positions shift by two.
CHAINS = """chain 1 MT 16569 + 0 16569 MT 16569 + 0 16569 1
16569

chain 1 chrM 16569 + 0 16569 MT 16569 + 0 16569 2
16569

chain 1 chrM 16571 + 0 16571 chrM 16569 + 0 16569 3
12990 2 0
3579
"""

METADATA = dict(genome_build="GRCh37", taxon_id=9606, species="homo_sapiens",
                coordinates="1-based-inclusive", dataset_version="fixture-v1",
                sources=[dict(id="synthetic", scope="test fixture")])


def record(orf_id, contig):
    model = ORFModel(orf_id, "tx", "gene", "GENE", contig, "+", ((100, 111),), "MKP", "lncRNA",
                     "GRCh37", 9606, "homo_sapiens", source="synthetic")
    return ORFRecord.from_model(model, ORFEvidence())


def test_grch37_chains_keep_hg19_chrm_apart_from_the_rcrs(tmp_path):
    path = tmp_path / "GRCh37_to_GRCh38.chain"
    path.write_text(CHAINS)
    chain = ChainFile(path, source_build="GRCh37", target_build="GRCh38")
    # Ensembl's two names for the rCRS are one chain.
    assert chain.map("MT", 12995, 12995) == ("MT", "+", 12995, 12995)
    # hg19's chrM lifts only through its own chain.
    assert chain.map("chrM", 12995, 12995) == ("MT", "+", 12993, 12993)
    assert chain.map("M", 12995, 12995) == ("MT", "+", 12993, 12993)
    # hg38's chrM is the rCRS.
    path.write_text("chain 1 chrM 16569 + 0 16569 MT 16569 + 0 16569 1\n16569\n")
    back = ChainFile(path, source_build="GRCh38", target_build="GRCh37")
    assert back.map("chrM", 12994, 12994) == back.map("MT", 12994, 12994) == ("MT", "+", 12994, 12994)


def test_grch37_queries_and_merges_keep_hg19_chrm_apart_from_mt():
    catalog = ORFCatalog([record("rcrs", "MT"), record("hg19", "chrM")], METADATA)
    for contig, expected in (("MT", "rcrs"), ("chrMT", "rcrs"), ("chrM", "hg19"), ("M", "hg19")):
        found = catalog.query(contig=contig, start=100, end=100, genome_build="GRCh37")
        assert [r.orf_id for r in found] == [expected], contig
    # The same path and protein on the two sequences are two ORFs.
    other = ORFCatalog([replace(record("hg19-other", "chrM"), source="other")],
                       dict(METADATA, dataset_version="other-v1", sources=[dict(id="other", scope="t")]))
    merged = merge_catalogs([ORFCatalog([record("rcrs", "MT")], METADATA), other])
    assert sorted(r.orf_id for r in merged.records) == ["hg19-other", "rcrs"]
