"""The whole gdc-viruses collection, built offline from its pinned real source records.

Opt in with ``pytest --real-data``: the one download is the pinned archive from
GitHub (``orfdata data fetch gdc-viruses-records``), never NCBI or PaVE. This is
the check the release process runs, because defects in real records (a protein
ID shared after sanitizing, a lineage whose last name is not the genus, a
slippage CDS that reads a base twice) have not shown in the synthetic suite.
"""

import gzip

import pytest
from conftest import serve_records

from orfdata import download_viral_references, viral_references
from orfdata.source_records import install_records, read_records

pytestmark = pytest.mark.real_data

NAME = "gdc-viruses-records"
REFERENCES = 201
ORFS = 1758


@pytest.fixture(scope="module")
def collection(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("full-collection"))
    install_records(NAME, cache_root=root)
    genbank, alternate = {}, {}
    for filename, raw in read_records(NAME, cache_root=root).items():
        source, identifier = filename.removesuffix(".gz").split("__", 1)
        if source == "genbank":
            accession = identifier.removesuffix(".gb")
            genbank[accession] = genbank[accession.split(".")[0]] = gzip.decompress(raw).decode()
        else:
            alternate[identifier.removesuffix(".json")] = gzip.decompress(raw)
    patch = pytest.MonkeyPatch()
    serve_records(patch, genbank, alternate)
    try:
        download_viral_references(species="human", cache_root=root)
    finally:
        patch.undo()
    return root


def test_the_collection_offers_the_published_counts(collection):
    references = viral_references(cache_root=collection)
    assert len(references) == REFERENCES
    assert sum(len(item.orfs) for item in references) == ORFS


def test_eight_source_cds_are_left_out_and_each_for_a_known_reason(collection):
    """Three ambiguous bases, two starts off the source protein, two HTLV-1 slippages, one EBV splice."""
    import re

    left_out = sum(int(match.group(1)) for item in viral_references(cache_root=collection)
                   if (match := re.search(r"(\d+) of \d+ source CDS not offered", item.scope)))
    assert left_out == 8
