"""Sequence normalization for BLOSUM62 protein homology searches."""

BLOSUM62_ALPHABET = frozenset("*ABCDEFGHIKLMNPQRSTVWXYZ")
_EQUIVALENTS = str.maketrans({"U": "C", "O": "K", "J": "X"})


def blosum62_sequence(sequence):
    """Return a BLOSUM62-compatible sequence without changing its length.

    Ensembl proteins can contain selenocysteine (U), while Biopython's
    BLOSUM62 matrix has no U row. Use its standard amino-acid analogue C;
    likewise map pyrrolysine O to K and ambiguous leucine/isoleucine J to X.
    """
    normalized = sequence.upper().translate(_EQUIVALENTS)
    unsupported = set(normalized) - BLOSUM62_ALPHABET
    if unsupported:
        raise ValueError(f"sequence contains residues unsupported by BLOSUM62: {sorted(unsupported)}")
    return normalized
