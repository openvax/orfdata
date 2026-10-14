"""Lift catalogue ORFs between assemblies of one species with a chain file.

Liftover is an explicit catalogue build, never a side effect of loading. A
record lifts only when every coding block lies inside one ungapped alignment
block, and all blocks land on one contig and strand in order; otherwise it is
left out with the reason. A lifted record keeps the original record in
``lifted_from`` and the chain file's identity in ``liftover``. Placement on the
target genome then re-checks the exact coding path and protein, so a lift that
changes the sequence is excluded rather than trusted.

Chain files are never downloaded implicitly. :func:`download_chain` fetches
one of the pinned files listed by :func:`pinned_chains` and verifies its checksum;
:func:`load_chain` then opens it offline.
"""

import gzip
import hashlib
import json
import re
from bisect import bisect_left
from dataclasses import asdict, dataclass, replace
from importlib.resources import files
from pathlib import Path


@dataclass(frozen=True)
class Liftover:
    """The chain file that moved an ORF from one assembly to another."""

    source_build: str
    target_build: str
    chain: str
    sha256: str


# UCSC names unlocalized and unplaced scaffolds chr1_GL456210v1_random or chrUn_JH584304v1; Ensembl GL456210.1.
_UCSC_SCAFFOLD = re.compile(r"chr(?:[0-9]+|[XYM]|Un)_([A-Za-z]+[0-9]+)v([0-9]+)(?:_random)?")


def _chain_contig(name):
    found = _UCSC_SCAFFOLD.fullmatch(name)
    return name if found is None else f"{found.group(1)}.{found.group(2)}"


class _NotLifted(ValueError):
    """Why a record did not lift; ``contigs`` names the places a block aligns to, when several."""

    def __init__(self, reason, contigs=()):
        super().__init__(reason)
        self.contigs = tuple(contigs)


class ChainFile:
    """Ungapped alignment blocks from a UCSC or Ensembl chain file.

    ``source_build`` is the assembly the chain maps from (UCSC's "target" or
    reference side) and ``target_build`` the one it maps to. Ensembl's
    ``GRCh37_to_GRCh38.chain.gz`` and UCSC's ``mm10ToMm39.over.chain.gz`` both
    work; contig names are compared as PyEnsembl normalizes them.
    ``contig_names`` renames the chain's target contigs to the target genome's,
    for example from the RefSeq accessions UCSC uses for GCF assemblies, using
    the NCBI assembly report. Check each file's terms of use before
    redistributing what it produces.
    """

    def __init__(self, path, *, source_build, target_build, name=None, contig_names=None):
        from .model import _contig

        raw = Path(path).read_bytes()
        self.liftover = Liftover(source_build, target_build, name or Path(path).name,
                                 hashlib.sha256(raw).hexdigest())
        text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode()
        renamed = dict(contig_names or {})
        blocks = {}
        header = None
        for line in text.splitlines():
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if fields[0] == "chain":
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                if fields[4] != "+":
                    raise ValueError(f"{self.liftover.chain} has a chain whose source strand is {fields[4]}; "
                                     "chain files align the source's forward strand")
                header = (_contig(_chain_contig(fields[2])),
                          _contig(_chain_contig(renamed.get(fields[7], fields[7]))), fields[9], int(fields[8]))
                t, q = int(fields[5]), int(fields[10])
                continue
            size = int(fields[0])
            source, contig, strand, length = header
            blocks.setdefault(source, set()).add((t, t + size, contig, strand, q, length))
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
        ``"-"`` when it reverses it. None unless the whole interval lies in one
        ungapped block. Where a region also aligns to fix patches or alternate
        haplotypes, the block on the same-named contig is used.
        """
        hits, _ = self.locate(contig, start, end)
        return hits[0] if len(hits) == 1 else None

    def locate(self, contig, start, end):
        """Every placement of a whole 1-based inclusive interval, and whether a block covers only part of it."""
        from .model import _contig

        source = _contig(contig)
        starts, widest, items = self._index.get(source, ((), 0, ()))
        low, high = start - 1, end
        hits, partial = [], False
        for t, t_end, target, strand, q, length in items[bisect_left(starts, low - widest) : bisect_left(starts, high)]:
            if t <= low and high <= t_end:
                a, b = q + low - t, q + high - t
                hits.append((target, strand) + ((length - b, length - a) if strand == "-" else (a, b)))
            elif t < high and low < t_end:
                partial = True
        hits = [(target, strand, a + 1, b) for target, strand, a, b in hits]
        same = [h for h in hits if h[0] == source]
        return (same if len(hits) > 1 and len(same) == 1 else hits), partial


@dataclass(frozen=True)
class PinnedChain:
    """A published chain file orfdata can download, pinned by URL, size and SHA-256."""

    name: str
    species: str
    source_build: str
    target_build: str
    url: str
    sha256: str
    size: int
    provider: str
    terms: str
    terms_url: str
    restricted: bool

    @property
    def filename(self):
        return self.url.rsplit("/", 1)[-1]


def _pins():
    return json.loads(files("orfdata").joinpath("data/chains.json").read_text())


def pinned_chains(*, species=None):
    """The pinned chain files, without downloading any."""
    from .catalog import _species

    found = tuple(PinnedChain(name, **entry) for name, entry in sorted(_pins().items()))
    return found if species is None else tuple(c for c in found if c.species == _species(species))


def _pinned(source_build, target_build, species):
    found = [c for c in pinned_chains(species=species) if (c.source_build, c.target_build) == (source_build, target_build)]
    if not found:
        known = ", ".join(c.name for c in pinned_chains(species=species)) or "none"
        raise ValueError(f"No pinned chain file maps {source_build} to {target_build} (pinned: {known}); "
                         "pass your own file as ChainFile(path, source_build=..., target_build=...)")
    return found[0]


def _chain_registry(chain, cache_root):
    from .datasets import _datacache, _digest, get_data_dir

    assets = {chain.filename: dict(url=chain.url, sha256=chain.sha256, size=chain.size)}
    return _datacache().VersionedDatasetRegistry(
        {chain.name: dict(default_version="1", versions={"1": assets})},
        store_path=lambda name, version: get_data_dir(cache_root) / "chains" / name / _digest(assets))


def download_chain(source_build, target_build, *, species=None, accept_terms=False, force=False,
                   source_path=None, cache_root=None):
    """Download the pinned chain file from ``source_build`` to ``target_build``; return its path.

    The file's SHA-256 and size must match the pin. ``source_path`` installs a
    copy you already have instead of downloading, checked the same way. A
    chain with usage restrictions, such as UCSC's non-commercial terms, needs
    ``accept_terms=True``.
    """
    chain = _pinned(source_build, target_build, species)
    if chain.restricted and accept_terms is not True:
        raise ValueError(f"{chain.filename} comes from {chain.provider} with these terms: {chain.terms} "
                         f"({chain.terms_url}). Pass accept_terms=True (or --accept-terms) to download it.")
    registry = _chain_registry(chain, cache_root)
    registry.download(chain.name, "1", force=force,
                      source_paths=None if source_path is None else {chain.filename: str(source_path)})
    return registry.local_path(chain.name, "1", asset=chain.filename)


def load_chain(source_build, target_build, *, species=None, cache_root=None):
    """The downloaded pinned chain file from ``source_build`` to ``target_build``, opened offline."""
    chain = _pinned(source_build, target_build, species)
    registry = _chain_registry(chain, cache_root)
    if registry.inspect(chain.name, "1", verify_files=False).status == "missing":
        raise ValueError(f"The {chain.name} chain file is not downloaded; run "
                         f"download_chain({source_build!r}, {target_build!r}) or "
                         f"'orfdata chains download {source_build} {target_build}'")
    path = registry.local_path(chain.name, "1", asset=chain.filename)
    return ChainFile(path, source_build=source_build, target_build=target_build, name=chain.filename)


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
            # The record itself is kept, so placements can report it under the rule that selects it.
            entry = dict(orf_id=record.orf_id, reason=str(error), record=record.to_dict())
            if getattr(error, "contigs", ()):
                entry["contigs"] = list(error.contigs)
            unmapped.append(entry)
    original = catalog.metadata
    metadata = dict(species=original["species"], taxon_id=original["taxon_id"],
                    genome_build=chain.liftover.target_build, coordinates=original["coordinates"],
                    sources=original.get("sources", []))
    # The input catalogue's own description, including any earlier lift, is kept whole.
    metadata["liftover"] = dict(
        source=catalog.source, **asdict(chain.liftover), lifted=len(lifted), unmapped=unmapped,
        input=dict(sha256=hashlib.sha256(json.dumps(catalog.to_dict(), sort_keys=True).encode()).hexdigest(),
                   metadata={k: v for k, v in original.items() if k != "sources"}))
    return ORFCatalog(lifted, metadata)


def _lift(record, chain):
    from .reference import ReferenceContext

    mapped = []
    for a, b in record.coding_blocks:
        hits, partial = chain.locate(record.contig, a, b)
        if len(hits) > 1:
            raise _NotLifted("A coding block aligns to several places", sorted({h[0] for h in hits}))
        if not hits:
            raise _NotLifted("A coding block crosses an alignment gap or edge" if partial
                             else "A coding block is not aligned")
        mapped.append(hits[0])
    if len({(contig, strand) for contig, strand, _, _ in mapped}) != 1:
        raise _NotLifted("Coding blocks land on different contigs or strands",
                         sorted({contig for contig, _, _, _ in mapped}))
    contig, orientation = mapped[0][:2]
    blocks = tuple((a, b) for _, _, a, b in mapped)
    if orientation == "-":
        blocks = blocks[::-1]
    if any(b >= c for (_, b), (c, _) in zip(blocks, blocks[1:])):
        raise _NotLifted("Coding blocks land out of order")
    strand = record.strand if orientation == "+" else {"+": "-", "-": "+"}[record.strand]
    members = tuple(_lift(m, chain) for m in record.members)
    pinned = record.reference_context
    if pinned is not None:  # Keep species and strain; the annotation is gone, the chain pins the coordinates.
        pinned = ReferenceContext(pinned.taxon_id, pinned.species, chain.liftover.target_build, None,
                                  strain=pinned.strain, strain_taxon_id=pinned.strain_taxon_id,
                                  source_sha256=chain.liftover.sha256)
    lifted = replace(record, contig=contig, strand=strand, coding_blocks=blocks,
                     genome_build=chain.liftover.target_build, reference_context=pinned,
                     initiation_context=None, members=members, lifted_from=record, liftover=chain.liftover)
    context = record.initiation_context
    if context is not None:
        positions = [chain.map(record.contig, p, p) for p in context.genomic_positions]
        moved = tuple(m[2] for m in positions if m is not None and m[:2] == (contig, orientation))
        step = 1 if strand == "+" else -1
        if len(moved) == len(positions) and all((b - a) * step > 0 for a, b in zip(moved, moved[1:])) and all(
                (lifted._coding_offset(p) is None) if relative < 0 else lifted._coding_offset(p) == relative - 1
                for relative, p in zip(context.relative_positions, moved)):
            context = replace(context, genomic_positions=moved)
        else:  # The flanks did not move with the coding path; keep their sequence only.
            context = replace(context, genomic_positions=())
        lifted = replace(lifted, initiation_context=context)
    return lifted


def unlifted(metadata):
    """Every record a catalogue's lifts left out, through each hop and merged input: (record, assembly, reason)."""
    from .catalog import _record

    lift = metadata.get("liftover")
    if lift:
        yield from unlifted(lift["input"]["metadata"])
        for entry in lift["unmapped"]:
            if entry.get("record") is not None:  # Written by orfdata 0.23.0 or later.
                yield _record(entry["record"], False), lift["target_build"], entry["reason"]
    for merged in metadata.get("merged_from", ()):
        yield from unlifted(merged["metadata"])
