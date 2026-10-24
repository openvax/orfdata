"""Write coding transcripts as GTF, so other tools can read ORFdata's models.

A transcript here is a sequence of genomic exon blocks in transcript order
plus the offset of its first codon. That is enough to emit the feature rows
PyEnsembl indexes — ``transcript``, ``exon``, ``CDS``, ``start_codon`` and
``stop_codon`` — and through PyEnsembl it is what varcode, isovar, topiary and
vaxrank need to give a variant a coding consequence.

Both producers use it: viral CDS features from source GenBank records, and
noncanonical ORFs from an :class:`~orfdata.ORFCatalog`.
"""

GTF_SOURCE = "orfdata"
_CODING = "protein_coding"


def attributes(pairs):
    """One GTF attribute column; empty values are left out rather than quoted."""
    return " ".join(f'{key} "{value}";' for key, value in pairs if value not in (None, ""))


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
