"""Lift catalogue ORFs between assemblies of one species with a chain file.

Liftover is an explicit catalogue build, never a side effect of loading. A
record lifts only when every coding block lies inside one ungapped alignment
block, and all blocks land on one contig and strand in order; otherwise it is
left out with the reason. A lifted record keeps the original record in
``lifted_from`` and the chain file's identity in ``liftover``. Placement on the
target genome then re-checks the exact coding path and protein, so a lift that
changes the sequence is excluded rather than trusted.
"""

import gzip
import hashlib
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class Liftover:
    """The chain file that moved an ORF from one assembly to another."""

    source_build: str
    target_build: str
    chain: str
    sha256: str


class ChainFile:
    """Ungapped alignment blocks from a UCSC or Ensembl chain file.

    ``source_build`` is the assembly the chain maps from (UCSC's "target" or
    reference side) and ``target_build`` the one it maps to. Ensembl's
    ``GRCh37_to_GRCh38.chain.gz`` and UCSC's ``mm10ToMm39.over.chain.gz`` both
    work; contig names are compared as PyEnsembl normalizes them. Check each
    file's terms of use before redistributing what it produces.
    """

    def __init__(self, path, *, source_build, target_build, name=None):
        from .model import _contig

        raw = Path(path).read_bytes()
        self.liftover = Liftover(source_build, target_build, name or Path(path).name,
                                 hashlib.sha256(raw).hexdigest())
        text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode()
        blocks = {}
        header = None
        for line in text.splitlines():
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "chain":
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                header = (_contig(fields[2]), fields[7], fields[9], int(fields[8]))
                t, q = int(fields[5]), int(fields[10])
                continue
            size = int(fields[0])
            source, contig, strand, length = header
            blocks.setdefault(source, set()).add((t, t + size, _contig(contig), strand, q, length))
            if len(fields) == 3:
                t, q = t + size + int(fields[1]), q + size + int(fields[2])
        self._index = {}
        # Ensembl's files repeat each chain under UCSC and Ensembl contig names; the copies coincide here.
        for source, items in blocks.items():
            items = sorted(items)
            self._index[source] = ([s for s, *_ in items], max(e - s for s, e, *_ in items), items)

    def map(self, contig, start, end):
        """Where a 1-based inclusive interval lands: ``(contig, orientation, start, end)``.

        ``orientation`` is ``"+"`` when the target keeps the source's strand and
        ``"-"`` when it reverses it. None unless exactly one ungapped block
        contains the whole interval.
        """
        from .model import _contig

        starts, widest, items = self._index.get(_contig(contig), ((), 0, ()))
        low, high = start - 1, end
        found = [item for item in items[bisect_left(starts, low - widest) : bisect_right(starts, low)]
                 if item[0] <= low and high <= item[1]]
        if len(found) != 1:
            return None
        t, _, target, strand, q, length = found[0]
        a, b = q + low - t, q + high - t
        if strand == "-":
            a, b = length - b, length - a
        return target, strand, a + 1, b


def lift_catalog(catalog, chain):
    """A copy of ``catalog`` on ``chain``'s target assembly, with every record that lifts.

    The result's metadata lists the records left out and why under
    ``liftover["unmapped"]``. Its ``source`` is the original catalogue's
    followed by ``@`` and the target assembly.
    """
    from .catalog import ORFCatalog

    if catalog.genome_build != chain.liftover.source_build:
        raise ValueError(f"The chain maps from {chain.liftover.source_build}, not {catalog.genome_build}")
    lifted, unmapped = [], []
    for record in catalog.records:
        try:
            lifted.append(_lift(record, chain))
        except ValueError as error:
            unmapped.append(dict(orf_id=record.orf_id, reason=str(error)))
    metadata = catalog.metadata
    for key in ("catalog_manifest", "merged_from"):
        metadata.pop(key, None)
    metadata.update(genome_build=chain.liftover.target_build, liftover=dict(
        source=catalog.source, **asdict(chain.liftover), lifted=len(lifted), unmapped=unmapped))
    return ORFCatalog(lifted, metadata)


def _lift(record, chain):
    mapped = [chain.map(record.contig, a, b) for a, b in record.coding_blocks]
    if any(m is None for m in mapped):
        raise ValueError("A coding block is not inside one ungapped alignment block")
    if len({(contig, strand) for contig, strand, _, _ in mapped}) != 1:
        raise ValueError("Coding blocks land on different contigs or strands")
    contig, orientation = mapped[0][:2]
    blocks = tuple((a, b) for _, _, a, b in mapped)
    if orientation == "-":
        blocks = blocks[::-1]
    if any(b >= c for (_, b), (c, _) in zip(blocks, blocks[1:])):
        raise ValueError("Coding blocks land out of order")
    strand = record.strand if orientation == "+" else {"+": "-", "-": "+"}[record.strand]
    members = tuple(_lift(m, chain) for m in record.members)
    lifted = replace(record, contig=contig, strand=strand, coding_blocks=blocks,
                     genome_build=chain.liftover.target_build, reference_context=None,
                     initiation_context=None, members=members, lifted_from=record, liftover=chain.liftover)
    context = record.initiation_context
    if context is not None:
        positions = [chain.map(record.contig, p, p) for p in context.genomic_positions]
        moved = tuple(m[2] for m in positions if m is not None and m[:2] == (contig, orientation))
        if len(moved) == len(positions) and all(
                (lifted._coding_offset(p) is None) if relative < 0 else lifted._coding_offset(p) == relative - 1
                for relative, p in zip(context.relative_positions, moved)):
            context = replace(context, genomic_positions=moved)
        else:  # The flanks did not move with the coding path; keep their sequence only.
            context = replace(context, genomic_positions=())
        lifted = replace(lifted, initiation_context=context)
    return lifted
