"""Transcript exon and CDS structure read in bulk from PyEnsembl's SQLite index.

Building a proteome or placing every catalogue ORF needs the structure of
hundreds of thousands of transcripts. Reading it one transcript at a time
through PyEnsembl objects costs several queries each; these helpers read the
same index tables contig by contig instead. Nothing is downloaded or indexed.
"""

import sqlite3
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Structure:
    """One transcript's exons, in transcript order, and its annotated CDS path."""

    transcript_id: str
    versioned_id: str
    gene_id: str
    gene_name: Optional[str]
    biotype: Optional[str]
    contig: str
    strand: str
    exons: tuple
    offsets: tuple
    length: int
    cds_blocks: tuple = ()
    protein_id: Optional[str] = None
    protein_version: Optional[int] = None
    start_complete: bool = False
    stop_complete: bool = False

    @property
    def layout(self):
        return self.exons, self.offsets, self.length


def _rows(db, sql, params=()):
    try:
        return db.connection.execute(sql, params).fetchall()
    except sqlite3.OperationalError:  # For example no CDS table in a noncoding annotation.
        return []


def _column(db, table, column, default="NULL"):
    try:
        return column if db.column_exists(table, column) else default
    except sqlite3.OperationalError:
        return default


def _merge(ranges):
    merged = []
    for a, b in sorted(ranges):
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return tuple(tuple(x) for x in merged)


def _build(transcripts, exons, cds, starts, stops):
    for transcript_id, (version, gene_id, gene_name, biotype, contig, strand) in transcripts.items():
        pieces = exons.get(transcript_id)
        if not pieces:
            continue
        ordered = tuple(sorted(pieces, reverse=strand == "-"))
        offsets, total = [], 0
        for a, b in ordered:
            offsets.append(total)
            total += b - a + 1
        coding = cds.get(transcript_id)
        protein_id = protein_version = None
        blocks = ()
        if coding:
            protein_id, protein_version = coding[0][2], coding[0][3]
            blocks = _merge([(a, b) for a, b, _, _ in coding] + list(stops.get(transcript_id, ())))
        yield Structure(
            transcript_id, f"{transcript_id}.{version}" if version not in (None, "") else transcript_id,
            gene_id, gene_name or None, biotype, contig, strand, ordered, tuple(offsets), total, blocks,
            protein_id, int(protein_version) if protein_version not in (None, "") else None,
            transcript_id in starts, transcript_id in stops)


def structures(genome, *, transcript_ids=None):
    """Yield a :class:`Structure` per transcript.

    Without ``transcript_ids``, every transcript is read contig by contig;
    with them, only those transcripts are read. Raises ``AttributeError`` when
    the genome has no SQLite index to read, so callers can fall back to objects.
    """
    db = genome.db
    version = _column(db, "transcript", "transcript_version")
    gene_name = _column(db, "transcript", "gene_name")
    biotype = _column(db, "transcript", "transcript_biotype")
    protein_version = _column(db, "CDS", "protein_version")
    queries = dict(
        transcript=f"SELECT transcript_id, {version}, gene_id, {gene_name}, {biotype}, seqname, strand FROM transcript",
        exon="SELECT transcript_id, start, end FROM exon",
        cds=f"SELECT transcript_id, start, end, protein_id, {protein_version} FROM CDS",
        start_codon="SELECT transcript_id FROM start_codon",
        stop_codon="SELECT transcript_id, start, end FROM stop_codon",
    )
    if transcript_ids is None:
        groups = [("seqname = ?", (row[0],)) for row in _rows(db, "SELECT DISTINCT seqname FROM transcript")]
    else:
        ids = sorted(set(transcript_ids))
        groups = [(f"transcript_id IN ({','.join('?' * len(chunk))})", tuple(chunk))
                  for chunk in (ids[i : i + 500] for i in range(0, len(ids), 500))]
    for condition, params in groups:
        rows = {name: _rows(db, f"{sql} WHERE {condition}", params) for name, sql in queries.items()}
        transcripts = {row[0]: row[1:] for row in rows["transcript"]}
        exons, cds, stops = {}, {}, {}
        for transcript_id, a, b in rows["exon"]:
            exons.setdefault(transcript_id, []).append((a, b))
        for transcript_id, a, b, protein, protein_v in rows["cds"]:
            cds.setdefault(transcript_id, []).append((a, b, protein, protein_v))
        for transcript_id, a, b in rows["stop_codon"]:
            stops.setdefault(transcript_id, []).append((a, b))
        yield from _build(transcripts, exons, cds, {row[0] for row in rows["start_codon"]}, stops)


def transcripts_containing(genome, spans):
    """Transcript IDs whose genomic span contains each requested span.

    ``spans`` maps a key to ``(contig, strand, start, end)``, 1-based inclusive.
    Contig names are compared as PyEnsembl normalizes them. Raises
    ``AttributeError`` when the genome has no SQLite index to read.
    """
    from .model import _contig

    width = 1 << 16
    bins = {}
    for transcript_id, contig, strand, start, end in _rows(
            genome.db, "SELECT transcript_id, seqname, strand, start, end FROM transcript"):
        key = (_contig(contig), strand)
        for b in range(start // width, end // width + 1):
            bins.setdefault((key, b), []).append((start, end, transcript_id))
    found = {}
    for name, (contig, strand, start, end) in spans.items():
        candidates = bins.get(((_contig(contig), strand), start // width), ())
        found[name] = sorted({t for a, b, t in candidates if a <= start and b >= end})
    return found
