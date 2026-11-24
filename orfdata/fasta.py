"""Write and read FASTA records the one way, so identifiers stay the first token."""

import gzip
import os
from dataclasses import dataclass

LINE_WIDTH = 80


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record parsed from a file or text stream.

    Attributes
    ----------
    header : str
        Header line after the leading ``>``, with trailing line endings removed.
    sequence : str
        Concatenated sequence lines with surrounding whitespace stripped.
    line_number : int
        1-based line number of the ``>`` header line in the source.
    """

    header: str
    sequence: str
    line_number: int

    @property
    def identifier(self) -> str:
        """Return the first whitespace-delimited token of :attr:`header`."""
        return self.header.split(None, 1)[0]


def write_fasta(stream, identifier, sequence, *, description=""):
    """Write one FASTA record wrapped at :data:`LINE_WIDTH` characters.

    Parameters
    ----------
    stream : io.TextIOBase
        Writable text stream.
    identifier : str
        Record identifier, written as the first header token.
    sequence : str
        Nucleotide or amino-acid sequence.
    description : str, optional
        Optional header description written after a space following
        ``identifier``.
    """
    header = f"{identifier} {description}".rstrip()
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), LINE_WIDTH):
        stream.write(f"{sequence[offset:offset + LINE_WIDTH]}\n")


def iter_fasta_records(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Yield :class:`FastaRecord` entries from a FASTA path or open text stream.

    Strict checks reject headers without an identifier and non-blank sequence
    lines appearing before the first header.

    Parameters
    ----------
    source : str, os.PathLike, or iterable of str
        Path to a FASTA file or an already-open text stream / line iterable.
    upper : bool, optional
        Upper-case each sequence line as it is read.
    gzipped : bool, optional
        When ``source`` is a path, open it with :func:`gzip.open`.
    source_name : str, optional
        Display label used in error messages; defaults to ``str(source)`` for
        paths or ``getattr(source, "name", "<stream>")`` for streams.
    encoding : str, optional
        Text encoding used when opening a file path (default ``"utf-8"``).

    Yields
    ------
    FastaRecord
        Parsed record containing the full header, sequence, and 1-based header
        line number.

    Raises
    ------
    ValueError
        If a ``>`` header has no non-whitespace identifier or if non-blank
        sequence text appears before the first header.
    """
    is_path = isinstance(source, (str, os.PathLike))
    if is_path:
        label = str(source) if source_name is None else source_name
        opener = gzip.open if gzipped else open
        with opener(source, "rt", encoding=encoding) as stream:
            yield from iter_fasta_records(stream, upper=upper, source_name=label)
        return

    label = getattr(source, "name", "<stream>") if source_name is None else source_name
    header, parts, start_line = None, [], None
    for number, line in enumerate(source, 1):
        text = line.rstrip("\r\n")
        if text.startswith(">"):
            if header is not None:
                yield FastaRecord(header, "".join(parts), start_line)
            header = text[1:]
            if not header.strip():
                raise ValueError(f"FASTA header without an identifier: {label}:{number}")
            parts, start_line = [], number
        elif header is not None:
            chunk = text.strip()
            if chunk:
                parts.append(chunk.upper() if upper else chunk)
        elif text.strip():
            raise ValueError(f"FASTA sequence before the first header: {label}:{number}")
    if header is not None:
        yield FastaRecord(header, "".join(parts), start_line)


def iter_fasta(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Yield ``(identifier, sequence)`` for every record of a FASTA file or stream.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one and sequence before the first header are errors rather
    than silently dropped records.

    Parameters
    ----------
    source : str, os.PathLike, or iterable of str
        Path to a FASTA file or an already-open text stream / line iterable.
    upper : bool, optional
        Upper-case each sequence as it is read.
    gzipped : bool, optional
        When ``source`` is a path, open it with :func:`gzip.open`.
    source_name : str, optional
        Display label used in error messages.
    encoding : str, optional
        Text encoding used when opening a file path (default ``"utf-8"``).

    Yields
    ------
    tuple of (str, str)
        ``(identifier, sequence)`` pairs in file order.
    """
    for record in iter_fasta_records(
        source, upper=upper, gzipped=gzipped, source_name=source_name, encoding=encoding
    ):
        yield record.identifier, record.sequence


def read_fasta(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Read every record of a FASTA file or stream as ``{identifier: sequence}``.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one, an identifier used twice, and sequence before the
    first header are errors rather than silently dropped records.

    Parameters
    ----------
    source : str, os.PathLike, or iterable of str
        Path to a FASTA file or an already-open text stream / line iterable.
    upper : bool, optional
        Upper-case each sequence as it is read.
    gzipped : bool, optional
        When ``source`` is a path, open it with :func:`gzip.open`.
    source_name : str, optional
        Display label used in error messages.
    encoding : str, optional
        Text encoding used when opening a file path (default ``"utf-8"``).

    Returns
    -------
    dict of str to str
        Mapping from record identifier to sequence, in file order.

    Raises
    ------
    ValueError
        If a header lacks an identifier, an identifier repeats, or sequence
        text appears before the first header.
    """
    label = (
        str(source)
        if source_name is None and isinstance(source, (str, os.PathLike))
        else (getattr(source, "name", "<stream>") if source_name is None else source_name)
    )
    records = {}
    for name, sequence in iter_fasta(
        source, upper=upper, gzipped=gzipped, source_name=label, encoding=encoding
    ):
        if name in records:
            raise ValueError(f"FASTA identifier repeated: {name} in {label}")
        records[name] = sequence
    return records
