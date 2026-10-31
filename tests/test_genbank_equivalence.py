"""ORFdata's GenBank reader is checked against an independent implementation.

The flat-file format has more edges than the records ORFdata pins, so the
reader is compared with `gb-io`, a maintained parser built on a Rust
implementation of the format. It is a development dependency: the comparison
runs wherever it is installed and is skipped elsewhere, so a published install
needs no parser of its own.

The two do not agree on everything, and where they differ ORFdata keeps the
more faithful reading. gb-io reads an ``order(...)`` location as a join, losing
a distinction the format makes, and returns a location it cannot parse — such
as the ``complement(order(>104408,...))`` on an EBV gene — as an unparsed
string. Those are recorded here rather than treated as failures of either side.
"""

import gzip
import io
from pathlib import Path

import pytest
from test_genbank import RECORD, WRAPPED

from orfdata import genbank

gb_io = pytest.importorskip("gb_io", reason="the differential check needs the dev extra")

FIXTURE = Path(__file__).parent / "data/viral/NC_001526.4.gb.gz"
# Published records chosen for what they exercise: a circular genome whose
# polymerase crosses the origin, 5'-partial CDSs, a four-base partial feature,
# a 3'-partial one, and a record that annotates no CDS at all.
CORPUS = sorted((Path(__file__).parent / "data/viral/genbank").glob("*.gb.gz"))


# ``order(...)`` says the segments are in this order without asserting that they
# are joined. ORFdata keeps that distinction; gb-io reads both as a join, so the
# comparison allows it rather than calling the stricter model a difference.
_JOINED = {"join": "join", "order": "join"}


def _ours(location):
    """One comparable form of a parsed location."""
    kind = type(location).__name__
    if kind == "Range":
        return (location.start, location.end, location.before, location.after)
    if kind == "Complement":
        return ("complement", _ours(location.location))
    return (_JOINED[kind.lower()], tuple(_ours(part) for part in location.locations))


def _theirs(location):
    kind = type(location).__name__
    if kind == "Range":
        return (location.start + 1, location.end, bool(location.before), bool(location.after))
    if kind == "Complement":
        return ("complement", _theirs(location.location))
    return (_JOINED[kind.lower()], tuple(_theirs(part) for part in location.locations))


declined = []


def _compare(text):
    mine = genbank.parse(text)
    theirs = next(iter(gb_io.iter(io.BytesIO(text.encode("utf-8")))))
    assert mine.accession_version == theirs.version
    assert mine.sequence == theirs.sequence.decode("ascii").upper()
    assert [f.kind for f in mine.features] == [f.kind for f in theirs.features]
    for feature, other in zip(mine.features, theirs.features, strict=True):
        if type(other.location).__name__ == "External":
            # The oracle declined to parse this one; ORFdata's reading stands.
            declined.append((mine.accession_version, feature.kind, str(other.location)))
            continue
        assert _ours(feature.location) == _theirs(other.location), feature.kind
        assert list(feature.qualifiers) == [(q.key, q.value) for q in other.qualifiers], feature.kind
    return mine


def test_a_real_published_record_reads_the_same_way():
    record = _compare(gzip.decompress(FIXTURE.read_bytes()).decode("utf-8"))
    # The HPV16 RefSeq record: nine CDSs, one of them spliced.
    assert record.accession_version == "NC_001526.4"
    assert sum(feature.kind == "CDS" for feature in record.features) == 9


@pytest.mark.parametrize("text", [RECORD, WRAPPED])
def test_the_awkward_cases_read_the_same_way(text):
    _compare(text)


@pytest.mark.parametrize("path", CORPUS, ids=lambda path: path.name.split("__")[0])
def test_the_awkward_published_records_read_the_same_way(path):
    _compare(gzip.decompress(path.read_bytes()).decode("utf-8"))


def test_a_whole_collection_reads_the_same_way():
    """Every record under ``ORFDATA_GENBANK_CORPUS``, for checking a real collection.

    The committed records above run always; point this at a downloaded
    collection's ``records`` directory to compare the readers over all 200.
    """
    import os

    corpus = os.environ.get("ORFDATA_GENBANK_CORPUS")
    records = sorted(Path(corpus).rglob("*.gb.gz")) if corpus else CORPUS
    assert records, f"no GenBank records under {corpus}"
    declined.clear()
    for path in records:
        _compare(gzip.decompress(path.read_bytes()).decode("utf-8"))
    # Every location the oracle could not parse is one ORFdata reads; none is a CDS.
    assert all(kind != "CDS" for _, kind, _ in declined), declined
