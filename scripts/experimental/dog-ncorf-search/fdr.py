"""Small target-decoy FDR helpers shared by the dog ncORF search checks."""

import numpy as np


def qvalues(score, is_decoy):
    """Target-decoy q-values, assigning one value to every tied score.

    A score threshold includes the entire tie. Treating tied rows in their
    input order would otherwise make the result depend on whether a target or
    decoy happened to occur first.
    """
    score = np.asarray(score, dtype=float)
    is_decoy = np.asarray(is_decoy, dtype=bool)
    if score.ndim != 1 or is_decoy.ndim != 1 or len(score) != len(is_decoy):
        raise ValueError("score and is_decoy must be one-dimensional arrays of equal length")
    if len(score) == 0:
        return np.empty(0, dtype=float)
    if not np.isfinite(score).all():
        raise ValueError("scores must be finite")

    order = np.argsort(-score, kind="mergesort")
    ordered_score = score[order]
    decoys = np.cumsum(is_decoy[order])
    targets = np.cumsum(~is_decoy[order])
    fdr = decoys / np.maximum(targets, 1)

    # The FDR at a threshold is the value after every row with that score has
    # entered the accepted set. Broadcast that value across the whole tie.
    ends = np.r_[np.flatnonzero(ordered_score[:-1] != ordered_score[1:]), len(score) - 1]
    tied_fdr = np.empty(len(score), dtype=float)
    start = 0
    for end in ends:
        tied_fdr[start:end + 1] = fdr[end]
        start = end + 1
    ordered_q = np.minimum.accumulate(tied_fdr[::-1])[::-1]
    result = np.empty(len(score), dtype=float)
    result[order] = ordered_q
    return result
