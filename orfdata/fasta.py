"""Write FASTA records the one way, so identifiers stay the first token."""

LINE_WIDTH = 80


def write_fasta(stream, identifier, sequence, *, description=""):
    """Write one record; ``identifier`` is the lookup key tools read back."""
    header = f"{identifier} {description}".rstrip()
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), LINE_WIDTH):
        stream.write(f"{sequence[offset:offset + LINE_WIDTH]}\n")
