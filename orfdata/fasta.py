"""Write and read FASTA records via ``datacache.fasta``, keeping identifiers as the first token."""

from datacache.fasta import (
    LINE_WIDTH,
    FastaRecord,
    write_fasta,
)
from datacache.fasta import (
    iter_fasta as _iter_fasta,
)
from datacache.fasta import (
    iter_fasta_records as _iter_fasta_records,
)
from datacache.fasta import (
    read_fasta as _read_fasta,
)

__all__ = [
    "LINE_WIDTH",
    "FastaRecord",
    "iter_fasta",
    "iter_fasta_records",
    "read_fasta",
    "write_fasta",
]


def iter_fasta_records(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Yield :class:`FastaRecord` entries from a FASTA path or open text stream.

    Delegates to :func:`datacache.fasta.iter_fasta_records` while preserving
    ORFdata's explicit ``gzipped=False`` default.
    """
    return _iter_fasta_records(
        source,
        upper=upper,
        gzipped=gzipped,
        source_name=source_name,
        encoding=encoding,
    )


def iter_fasta(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Yield ``(identifier, sequence)`` for every record of a FASTA file or stream.

    Delegates to :func:`datacache.fasta.iter_fasta` while preserving ORFdata's
    explicit ``gzipped=False`` default.
    """
    return _iter_fasta(
        source,
        upper=upper,
        gzipped=gzipped,
        source_name=source_name,
        encoding=encoding,
    )


def read_fasta(source, *, upper=False, gzipped=False, source_name=None, encoding="utf-8"):
    """Read every record of a FASTA file or stream as ``{identifier: sequence}``.

    Delegates to :func:`datacache.fasta.read_fasta` while preserving ORFdata's
    explicit ``gzipped=False`` default.
    """
    return _read_fasta(
        source,
        upper=upper,
        gzipped=gzipped,
        source_name=source_name,
        encoding=encoding,
    )
