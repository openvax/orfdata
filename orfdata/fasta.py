"""Write and read FASTA records the one way, so identifiers stay the first token."""

import gzip

LINE_WIDTH = 80


def write_fasta(stream, identifier, sequence, *, description=""):
    """Write one record; ``identifier`` is the lookup key tools read back."""
    header = f"{identifier} {description}".rstrip()
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), LINE_WIDTH):
        stream.write(f"{sequence[offset:offset + LINE_WIDTH]}\n")


def iter_fasta(path, *, upper=False, gzipped=False):
    """Yield ``(identifier, sequence)`` for every record of a FASTA file.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one and sequence before the first header are errors rather
    than silently dropped records. ``upper`` upper-cases each sequence as it is
    read.
    """
    name, parts = None, []
    opener = gzip.open if gzipped else open
    with opener(path, "rt") as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError(f"FASTA header without an identifier: {path}:{number}")
                name, parts = tokens[0], []
            elif name is not None:
                parts.append(line.strip().upper() if upper else line.strip())
            elif line.strip():
                raise ValueError(f"FASTA sequence before the first header: {path}:{number}")
    if name is not None:
        yield name, "".join(parts)


def read_fasta(path, *, upper=False, gzipped=False):
    """Every record of a FASTA file as ``{identifier: sequence}``, in file order.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one, an identifier used twice, and sequence before the
    first header are errors rather than silently dropped records. ``upper``
    upper-cases each sequence as it is read.
    """
    records = {}
    for name, sequence in iter_fasta(path, upper=upper, gzipped=gzipped):
        if name in records:
            raise ValueError(f"FASTA identifier repeated: {name} in {path}")
        records[name] = sequence
    return records
