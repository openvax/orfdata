"""Write and read FASTA records the one way, so identifiers stay the first token."""

LINE_WIDTH = 80


def write_fasta(stream, identifier, sequence, *, description=""):
    """Write one record; ``identifier`` is the lookup key tools read back."""
    header = f"{identifier} {description}".rstrip()
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), LINE_WIDTH):
        stream.write(f"{sequence[offset:offset + LINE_WIDTH]}\n")


def read_fasta(path):
    """Every record of a FASTA file as ``{identifier: sequence}``, in file order.

    The identifier is the first header token, as :func:`write_fasta` writes it.
    A header without one, or an identifier used twice, is an error rather than a
    silently dropped record.
    """
    records, name, parts = {}, None, []

    def keep():
        if name is not None:
            records[name] = "".join(parts)

    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith(">"):
                keep()
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError(f"FASTA header without an identifier: {path}:{number}")
                if tokens[0] in records:
                    raise ValueError(f"FASTA identifier repeated: {tokens[0]} in {path}")
                name, parts = tokens[0], []
            else:
                parts.append(line.strip())
    keep()
    return records
