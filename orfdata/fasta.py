"""Write and read FASTA records the one way, so identifiers stay the first token."""

LINE_WIDTH = 80


def write_fasta(stream, identifier, sequence, *, description=""):
    """Write one record; ``identifier`` is the lookup key tools read back."""
    header = f"{identifier} {description}".rstrip()
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), LINE_WIDTH):
        stream.write(f"{sequence[offset:offset + LINE_WIDTH]}\n")


def read_fasta(path, *, upper=False):
    """Every record of a FASTA file as ``{identifier: sequence}``, in file order.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one, an identifier used twice, and sequence before the
    first header are errors rather than silently dropped records. ``upper``
    upper-cases each sequence as it is read.
    """
    records, name, parts = {}, None, []
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError(f"FASTA header without an identifier: {path}:{number}")
                if tokens[0] in records:
                    raise ValueError(f"FASTA identifier repeated: {tokens[0]} in {path}")
                name, parts = tokens[0], []
            elif name is not None:
                parts.append(line.strip().upper() if upper else line.strip())
            elif line.strip():
                raise ValueError(f"FASTA sequence before the first header: {path}:{number}")
    if name is not None:
        records[name] = "".join(parts)
    return records
