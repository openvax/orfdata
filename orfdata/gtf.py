"""Write coding transcripts as GTF, so other tools can read ORFdata's models, and read their exons back.

A transcript here is a sequence of genomic exon blocks in transcript order
plus the offset of its first codon. That is enough to emit the feature rows
PyEnsembl indexes — ``transcript``, ``exon``, ``CDS``, ``start_codon`` and
``stop_codon`` — and through PyEnsembl it is what varcode, isovar, topiary and
vaxrank need to give a variant a coding consequence.

Both producers use it: viral CDS features from source GenBank records, and
noncanonical ORFs from an :class:`~orfdata.ORFCatalog`.
"""

import gzip
import io
import os
import re
from collections import defaultdict
from contextlib import contextmanager

GTF_SOURCE = "orfdata"
_CODING = "protein_coding"


_PAIR = re.compile(r'([^\s;"]+)\s+(?:"([^"]*)"|([^;\s"]+))')


def _safe(value):
    """A value as one quoted GTF token.

    GTF has no escape for ``"``, and PyEnsembl's parser reads up to the next
    one, so a quote becomes an apostrophe and a tab or newline becomes a space.
    A semicolon is left alone, since quoted values may contain it.
    """
    return re.sub(r"[\t\r\n]", " ", str(value).replace('"', "'"))


def attributes(pairs):
    """One GTF attribute column; empty values are left out rather than quoted."""
    return " ".join(f'{key} "{_safe(value)}";' for key, value in pairs
                    if value not in (None, ""))


def parse_attribute_pairs(column):
    """Parse one GTF attribute column into ``(key, value)`` pairs in file order.

    This is the inverse of :func:`attributes`: quoted values may contain spaces
    or semicolons, trailing semicolons are optional, and repeated keys (such as
    Ensembl ``tag`` attributes) are preserved in order.

    Parameters
    ----------
    column : str
        Ninth tab-separated GTF column.

    Returns
    -------
    tuple of (str, str)
        Attribute ``(key, value)`` pairs in the order they appear in ``column``.
    """
    return tuple(
        (match.group(1), match.group(2) if match.group(2) is not None else match.group(3))
        for match in _PAIR.finditer(column)
    )


def parse_attributes(column, *, strict=False, repeated_keys=()):
    """Parse one GTF attribute column into a ``{key: value}`` dictionary.

    Parameters
    ----------
    column : str
        Ninth tab-separated GTF column.
    strict : bool, optional
        If ``True``, raise :class:`ValueError` when the same key appears more
        than once with different values (except keys listed in
        ``repeated_keys``). Defaults to ``False``.
    repeated_keys : iterable of str, optional
        Keys permitted to repeat with distinct values when ``strict=True``
        (for example ``("tag",)`` in Ensembl GTFs).

    Returns
    -------
    dict of str to str
        Mapping of attribute names to values.

    Raises
    ------
    ValueError
        If ``strict=True`` and a non-repeated key appears with conflicting
        values.
    """
    if not strict:
        return dict(parse_attribute_pairs(column))
    allowed = set(repeated_keys)
    parsed = {}
    for key, value in parse_attribute_pairs(column):
        if key not in allowed and key in parsed and parsed[key] != value:
            raise ValueError(f"Conflicting GTF attribute: {key}")
        parsed[key] = value
    return parsed


@contextmanager
def open_gtf(source):
    """Open a plain or gzipped GTF file as UTF-8 text, or yield an already-open text stream.

    Compression is detected from the gzip magic bytes (``\\x1f\\x8b``) rather
    than the filename suffix, so cached files without a ``.gz`` extension are
    handled identically.

    Parameters
    ----------
    source : str, os.PathLike, or io.TextIOBase
        Path to a plain or gzipped GTF file, or an already-open text stream.

    Yields
    ------
    io.TextIOBase
        Readable UTF-8 text stream over the GTF lines.
    """
    if isinstance(source, (str, os.PathLike)):
        with open(source, "rb") as raw, io.BufferedReader(raw) as buffered:
            decoded = gzip.GzipFile(fileobj=buffered) if buffered.peek(2).startswith(b"\x1f\x8b") else buffered
            with io.TextIOWrapper(decoded, encoding="utf-8") as text:
                yield text
    else:
        yield source


def read_gtf_rows(source, *, features=None, headers=None, strict_attributes=False,
                  repeated_attribute_keys=(), source_name=None):
    """Stream validated rows from a plain or gzipped GTF file or open text stream.

    Blank lines and ``#`` comments are skipped. When ``headers`` is provided,
    ``#!key value`` metadata lines are recorded into it and conflicting values
    for the same header key raise :class:`ValueError`. Every non-comment row is
    validated to have nine tab-separated fields before feature filtering so
    malformed rows are never silently skipped.

    Parameters
    ----------
    source : str, os.PathLike, or io.TextIOBase
        Path to a plain or gzipped GTF file, or an open text stream.
    features : iterable of str, optional
        Feature types (column 3) to yield, such as ``{"exon"}`` or
        ``{"transcript", "CDS"}``. When ``None``, every feature row is yielded.
    headers : dict, optional
        Mutable mapping populated in place with ``#!key value`` header entries
        encountered while streaming.
    strict_attributes : bool, optional
        Passed to :func:`parse_attributes`. Defaults to ``False``.
    repeated_attribute_keys : iterable of str, optional
        Passed as ``repeated_keys`` to :func:`parse_attributes`.
    source_name : str or os.PathLike, optional
        Display label used in error messages; defaults to ``source`` when a
        path is given or ``getattr(source, "name", "<stream>")`` for a stream.

    Yields
    ------
    tuple of (str, str, int, int, str, dict)
        ``(seqname, feature, start, end, strand, attributes)`` with 1-based
        inclusive integer coordinates.

    Raises
    ------
    ValueError
        On any row without nine tab-separated fields, missing sequence name,
        invalid coordinates, conflicting ``#!`` header, or conflicting strict
        attribute.
    """
    feature_set = None if features is None else set(features)
    label = (
        os.fspath(source_name)
        if source_name is not None
        else (os.fspath(source) if isinstance(source, (str, os.PathLike)) else getattr(source, "name", "<stream>"))
    )
    with open_gtf(source) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if headers is not None and line.startswith("#!"):
                    key, _, value = line[2:].strip().partition(" ")
                    value = value.strip()
                    if key in headers and headers[key] != value:
                        raise ValueError(f"Conflicting GTF header {key!r}: {label}:{number}")
                    headers[key] = value
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF row without nine tab-separated fields: {label}:{number}")
            if feature_set is not None and fields[2] not in feature_set:
                continue
            seqname = fields[0]
            if not seqname:
                raise ValueError(f"GTF row without a sequence name: {label}:{number}")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                raise ValueError(f"GTF row without integer coordinates: {label}:{number}") from None
            if start < 1 or end < start:
                raise ValueError(f"GTF row with invalid coordinates ({start}..{end}): {label}:{number}")
            attrs = parse_attributes(
                fields[8],
                strict=strict_attributes,
                repeated_keys=repeated_attribute_keys,
            )
            yield seqname, fields[2], start, end, fields[6], attrs


def transcript_blocks(exons, offset, length):
    """Genomic blocks covering a transcript interval, with each block's coding offset.

    ``offset`` and ``length`` are transcript coordinates; the blocks come back
    in transcript order as ``(start, end, strand, coding_offset)``, so a GTF
    frame is ``-coding_offset % 3``.
    """
    blocks, position = [], 0
    for start, end, strand in exons:
        span = end - start + 1
        begin, finish = max(offset, position), min(offset + length, position + span)
        if begin < finish:
            lead, tail = begin - position, position + span - finish
            block = ((start + lead, end - tail) if strand == "+" else (start + tail, end - lead))
            blocks.append((block[0], block[1], strand, begin - offset))
        position += span
    return tuple(blocks)


def _row(seqname, feature, start, end, strand, frame, pairs):
    return (f"{seqname}\t{GTF_SOURCE}\t{feature}\t{start}\t{end}\t.\t{strand}\t{frame}\t"
            f"{attributes(pairs)}\n")


def gene_line(seqname, gene_id, gene_name, exons, *, biotype=_CODING):
    """One ``gene`` row spanning every exon given."""
    return _row(seqname, "gene", min(start for start, _, _ in exons),
                max(end for _, end, _ in exons), exons[0][2], ".",
                (("gene_id", gene_id), ("gene_name", gene_name), ("gene_biotype", biotype)))


def transcript_lines(seqname, exons, *, gene_id, gene_name, transcript_id, protein_id,
                     translation_table=1, coding_offset=0, has_stop_codon=True,
                     transcript_name=None, biotype=_CODING, extra=()):
    """Every row placing one coding transcript on ``seqname``.

    ``exons`` are ``(start, end, strand)`` blocks, 1-based and inclusive, in
    transcript order. ``coding_offset`` is how many transcript bases precede
    the first codon, as GenBank's ``codon_start`` records. Following Ensembl,
    the ``CDS`` rows stop before the stop codon, which is its own feature;
    ``has_stop_codon=False`` leaves out both codon features, for a partial or
    frameshifted model whose reading frame does not close.
    """
    strand = exons[0][2]
    shared = (
        ("gene_id", gene_id),
        ("gene_name", gene_name),
        ("gene_biotype", biotype),
        ("transcript_id", transcript_id),
        ("transcript_name", transcript_name or transcript_id),
        ("transcript_biotype", biotype),
        ("transl_table", str(translation_table)),
    ) + tuple(extra)
    length = sum(end - start + 1 for start, end, _ in exons)
    rows = [_row(seqname, "transcript", min(start for start, _, _ in exons),
                 max(end for _, end, _ in exons), strand, ".", shared)]
    for number, (start, end, _) in enumerate(exons, 1):
        rows.append(_row(seqname, "exon", start, end, strand, ".", shared + (
            ("exon_id", f"{transcript_id}.e{number}"), ("exon_number", str(number)))))
    coding_length = length - coding_offset - (3 if has_stop_codon else 0)
    for start, end, _, offset in transcript_blocks(exons, coding_offset, coding_length):
        rows.append(_row(seqname, "CDS", start, end, strand, str(-offset % 3),
                         shared + (("protein_id", protein_id),)))
    if has_stop_codon:
        for feature, begin in (("start_codon", coding_offset), ("stop_codon", length - 3)):
            for start, end, _, offset in transcript_blocks(exons, begin, 3):
                rows.append(_row(seqname, feature, start, end, strand, str(-offset % 3), shared))
    return tuple(rows)


def orf_gtf_lines(records):
    """GTF rows for catalogue ORF records, in their coding-block order.

    Each ORF becomes one transcript named by its ``orf_id`` on the host gene,
    with the host transcript it was placed on kept as an attribute. No ``gene``
    rows are written, so the result can be concatenated with the reference
    annotation of the same assembly — or read on its own, since PyEnsembl
    fills in the genes a GTF leaves out.

    The rows carry ``transcript_biotype "protein_coding"``, because each ORF
    has a complete CDS with start and stop codons and effect predictors only
    translate transcripts marked that way. The source's own category stays in
    ``orf_biotype``.
    """
    lines = []
    for record in records:
        blocks = tuple(record.coding_blocks)
        exons = tuple((start, end, record.strand)
                      for start, end in (blocks if record.strand == "+" else blocks[::-1]))
        lines.extend(transcript_lines(
            record.contig, exons,
            gene_id=record.gene_id or record.orf_id,
            gene_name=record.gene_name or record.gene_id or record.orf_id,
            transcript_id=record.orf_id,
            protein_id=f"{record.orf_id}.p",
            translation_table=record.translation_table,
            extra=(("host_transcript_id", record.transcript_id),
                   ("orf_biotype", record.biotype),
                   ("orf_source", getattr(record, "source", "") or ""))))
    return tuple(lines)


def read_exons(path):
    """Exon blocks of every transcript in a GTF, as ``{transcript_id: (blocks, strand)}``.

    Blocks are ``(start, end)``, 1-based and inclusive, in transcript order as
    ``exon_number`` gives it (which is not coordinate order for a minus-strand
    transcript or one joined across a circular origin). This reads back what
    :func:`transcript_lines` writes, from plain or gzipped GTF files.

    Parameters
    ----------
    path : str, os.PathLike, or io.TextIOBase
        Plain or gzipped GTF file path, or an open text stream.

    Returns
    -------
    dict of str to (tuple of (int, int), str)
        Mapping from ``transcript_id`` to ``(exon_blocks, strand)``.
    """
    rows = {}
    for _, _, start, end, strand, attrs in read_gtf_rows(path, features={"exon"}):
        try:
            row = (int(attrs["exon_number"]), start, end, strand)
            name = attrs["transcript_id"]
        except (KeyError, ValueError):
            raise ValueError(
                f"GTF exon without a readable transcript_id, exon_number and coordinates: {path}"
            ) from None
        rows.setdefault(name, []).append(row)
    exons = {}
    for name, blocks in rows.items():
        numbers = [block[0] for block in blocks]
        if len(set(numbers)) != len(numbers) or len({block[3] for block in blocks}) != 1:
            raise ValueError(f"GTF transcript {name} repeats an exon number or mixes strands: {path}")
        exons[name] = (tuple((start, end) for _, start, end, _ in sorted(blocks)), blocks[0][3])
    return exons


def read_transcript_exons(source, *, transcript_ids=None):
    """Group GTF ``exon`` rows by transcript ID with exons in 5'-to-3' transcript order.

    Parameters
    ----------
    source : str, os.PathLike, or io.TextIOBase
        Plain or gzipped GTF path, or an open text stream.
    transcript_ids : iterable of str, optional
        If provided, retain only exons whose ``transcript_id`` attribute is in
        this set.

    Returns
    -------
    dict of str to dict
        Mapping ``transcript_id -> {"contig", "strand", "version", "gene_id",
        "gene_name", "exons"}`` where ``exons`` is a tuple of 1-based inclusive
        ``(start, end)`` pairs ordered ascending for ``"+"`` and descending for
        ``"-"``.
    """
    wanted = None if transcript_ids is None else set(transcript_ids)
    found = {}
    for seqname, _, start, end, strand, attrs in read_gtf_rows(source, features={"exon"}):
        transcript = attrs.get("transcript_id")
        if transcript is None or (wanted is not None and transcript not in wanted):
            continue
        entry = found.setdefault(
            transcript,
            dict(
                contig=seqname,
                strand=strand,
                version=attrs.get("transcript_version"),
                gene_id=attrs["gene_id"],
                gene_name=attrs.get("gene_name", ""),
                exons=[],
            ),
        )
        entry["exons"].append((start, end))
    for entry in found.values():
        entry["exons"] = tuple(sorted(entry["exons"], reverse=entry["strand"] == "-"))
    return found


def read_transcript_loci(source, *, genome_build=None, transcript_ids=None):
    """Group GTF ``exon`` rows into distinct versioned transcript loci, rejecting overlapping exons.

    Unlike :func:`read_transcript_exons`, this preserves full versioned
    transcript identifiers (appending ``transcript_version`` when present and
    not already part of ``transcript_id``) and keeps multiple annotated loci
    for the same versioned transcript separate (for example pseudoautosomal
    copies).

    Parameters
    ----------
    source : str, os.PathLike, or io.TextIOBase
        Plain or gzipped GTF path, or an open text stream.
    genome_build : str, optional
        When supplied, contig names are normalized with
        :func:`orfdata.model.normalize_contig`.
    transcript_ids : iterable of str, optional
        If provided, retain only loci whose versioned transcript ID is in this
        set.

    Returns
    -------
    dict of str to list of dict
        Mapping ``versioned_transcript_id -> [locus, ...]`` where each locus
        has keys ``"transcript"``, ``"contig"``, ``"strand"``, ``"gene_id"``,
        ``"gene_name"``, and ``"exons"`` (tuple in 5'-to-3' transcript order).

    Raises
    ------
    ValueError
        If a locus contains overlapping exons or any row in ``source`` is
        malformed.
    """
    from .model import normalize_contig

    wanted = None if transcript_ids is None else set(transcript_ids)
    groups = {}
    for seqname, _, start, end, strand, attrs in read_gtf_rows(source, features={"exon"}):
        name = attrs["transcript_id"]
        version = attrs.get("transcript_version")
        if version and "." not in name:
            name += "." + version
        if wanted is not None and name not in wanted:
            continue
        contig = normalize_contig(seqname, genome_build) if genome_build is not None else seqname
        key = name, contig, strand, attrs["gene_id"]
        entry = groups.setdefault(
            key,
            dict(
                transcript=name,
                contig=contig,
                strand=strand,
                gene_id=key[3],
                gene_name=attrs.get("gene_name", ""),
                exons=set(),
            ),
        )
        entry["exons"].add((start, end))
    found = defaultdict(list)
    for entry in groups.values():
        entry["exons"] = tuple(sorted(entry["exons"], reverse=entry["strand"] == "-"))
        ascending = sorted(entry["exons"])
        if any(b >= c for (_, b), (c, _) in zip(ascending, ascending[1:])):
            raise ValueError("Overlapping/contradictory reference exons: " + entry["transcript"])
        found[entry["transcript"]].append(entry)
    return found


def transcript_interval_blocks(exons, strand, start, end):
    """Map a 1-based inclusive transcript interval ``start..end`` to ascending genomic blocks.

    Parameters
    ----------
    exons : sequence of (int, int)
        1-based inclusive genomic ``(start, end)`` exon intervals in 5'-to-3'
        transcript order.
    strand : str
        Transcript strand (``"+"`` or ``"-"``).
    start : int
        1-based inclusive start coordinate in spliced transcript space.
    end : int
        1-based inclusive end coordinate in spliced transcript space.

    Returns
    -------
    tuple of (int, int) or None
        Ascending 1-based inclusive genomic blocks covering ``start..end``, or
        ``None`` if ``end`` extends beyond the spliced exon length.
    """
    blocks, offset = [], 0
    for a, b in exons:
        length = b - a + 1
        lo, hi = max(start, offset + 1), min(end, offset + length)
        if lo <= hi:
            if strand == "+":
                blocks.append((a + lo - offset - 1, a + hi - offset - 1))
            else:
                blocks.append((b - (hi - offset - 1), b - (lo - offset - 1)))
        offset += length
    return tuple(sorted(blocks)) if offset >= end else None
