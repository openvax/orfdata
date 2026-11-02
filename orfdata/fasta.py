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
    """
    records, name, parts = {}, None, []
    with open(path) as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records
