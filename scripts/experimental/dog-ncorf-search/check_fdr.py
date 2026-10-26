"""Self-test of class-specific target-decoy q-values; needs numpy only.

usage: python -I check_fdr.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdr import qvalues  # noqa: E402, I001


scores = np.array([10.0, 9.0, 9.0, 8.0])
decoys = np.array([False, False, True, True])
expected = np.array([0.0, 0.5, 0.5, 1.0])
np.testing.assert_allclose(qvalues(scores, decoys), expected)

# Reordering a target/decoy tie must not change either row's value.
permuted = np.array([0, 2, 1, 3])
restored = np.empty(4)
restored[permuted] = qvalues(scores[permuted], decoys[permuted])
np.testing.assert_allclose(restored, expected)
np.testing.assert_array_equal(qvalues([], []), np.array([]))

for bad_score, bad_decoy in (([1], []), ([[1]], [[False]]), ([float("nan")], [False])):
    try:
        qvalues(bad_score, bad_decoy)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid q-value input was accepted")

print("FDR checks passed")
