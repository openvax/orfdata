"""Self-test of BLOSUM62 sequence normalization; needs no external data.

usage: python -I check_homology.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homology import blosum62_sequence  # noqa: E402, I001


assert blosum62_sequence("ACUOJX*") == "ACCKXX*"
try:
    blosum62_sequence("AC-")
except ValueError as error:
    assert "'-'" in str(error)
else:
    raise AssertionError("an unsupported alignment residue was accepted")

print("homology checks passed")
