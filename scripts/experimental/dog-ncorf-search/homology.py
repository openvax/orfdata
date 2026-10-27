"""Sequence normalization for BLOSUM62 protein homology searches."""

from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal

BLOSUM62_ALPHABET = frozenset("*ABCDEFGHIKLMNPQRSTVWXYZ")
_EQUIVALENTS = str.maketrans({"U": "C", "O": "K", "J": "X"})
_BLOSUM62 = SubstitutionMatrix.std_protein_matrix()


def blosum62_sequence(sequence):
    """Return a BLOSUM62-compatible sequence without changing its length.

    Ensembl proteins can contain selenocysteine (U), while the standard
    BLOSUM62 matrix has no U row. Use its standard amino-acid analogue C;
    likewise map pyrrolysine O to K and ambiguous leucine/isoleucine J to X.
    """
    normalized = sequence.upper().translate(_EQUIVALENTS)
    unsupported = set(normalized) - BLOSUM62_ALPHABET
    if unsupported:
        raise ValueError(f"sequence contains residues unsupported by BLOSUM62: {sorted(unsupported)}")
    return normalized


def local_alignment(query, target, *, max_alignments=10000):
    """BLOSUM62 local alignment, preserving the previous Gotoh tie order.

    A length-n gap costs -11 - (n-1); local terminal overhangs are unpenalized.
    Equal optima choose the earliest (query end, target end), then traceback
    states in diagonal, query-only, target-only order. This is the ordering of
    the previous PairwiseAligner, independent of Biotite's enumeration order.
    Refuse a truncated set of ties rather than silently changing identity or
    aligned length. No positive-scoring alignment returns None.
    """
    if type(max_alignments) is not int or max_alignments < 1:
        raise ValueError("max_alignments must be a positive integer")
    query, target = blosum62_sequence(query), blosum62_sequence(target)
    if not query or not target:
        return None
    alignments = align_optimal(ProteinSequence(query), ProteinSequence(target), _BLOSUM62,
                               gap_penalty=(-11, -1), terminal_penalty=True,
                               local=True, max_number=max_alignments + 1)
    if not alignments or alignments[0].score <= 0:
        return None
    if len(alignments) > max_alignments:
        raise ValueError("Too many optimal local alignments to preserve complete tie handling")

    def key(alignment):
        trace = alignment.trace
        states = tuple(0 if q >= 0 and t >= 0 else 1 if q >= 0 else 2
                       for q, t in trace[::-1])
        return tuple(trace[-1]), states

    return min(alignments, key=key)
