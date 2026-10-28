import importlib.util
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts" / "experimental" / "dog-ncorf-search"
pytestmark = pytest.mark.skipif(not SCRIPTS.is_dir(), reason="experimental scripts are not in the sdist")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dog_search_qvalues_include_complete_score_ties():
    qvalues = _load("fdr").qvalues
    scores = np.array([10.0, 9.0, 9.0, 8.0])
    decoys = np.array([False, False, True, True])
    expected = np.array([0.0, 0.5, 0.5, 1.0])
    np.testing.assert_allclose(qvalues(scores, decoys), expected)
    np.testing.assert_allclose(qvalues(scores[[0, 2, 1, 3]], decoys[[0, 2, 1, 3]]), expected)


def test_dog_search_normalizes_non_blosum_residues():
    normalize = _load("homology").blosum62_sequence
    assert normalize("ACUOJX*") == "ACCKXX*"
    with pytest.raises(ValueError, match="unsupported"):
        normalize("AC-")


@pytest.mark.parametrize("query,target", [
    ("AAAA", "AAAAAAAA"),
    ("ACDEFGACDEFG", "ACDEFG"),
    ("ACDEFWWWWGHIK", "ACDEFGHIK"),
    ("ACDEFAAAAGHIK", "ACDEFAAAAAGHIK"),
    ("AX*", "AX*"),
    ("WY", "AA"),
    ("TRFSATV*TYS*XAEDEPIRPDFQCANIMSGKAARMKAXRKAL",
     "TRFSATV*TIRYS*XAEDEPXACQANIMAARMKAXRKAL"),
    ("MARYAGGSNDGHHRXRRIYSSLSSNVFWDAFW*MARAPXAVVAWCVQRELAAQT",
     "MARYAGGSNDGHHRXRRIYSSLSSNVLFALWAFW*MARAPXAVVAWCVQRELAAQT"),
])
def test_biotite_local_alignment_preserves_scores_gaps_and_ties(query, target):
    from Bio import Align
    from Bio.Align import substitution_matrices

    previous = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                     open_gap_score=-11, extend_gap_score=-1)
    expected = previous.align(query, target)
    actual = _load("homology").local_alignment(query, target)
    if not expected:
        assert actual is None
        return
    expected = expected[0]
    assert actual.score == expected.score
    assert [(q, t) for q, t in actual.trace if q >= 0 and t >= 0] == [
        (qs + k, ts + k) for (qs, qe), (ts, te) in zip(*expected.aligned) for k in range(qe - qs)]


def test_biotite_alignment_matches_independent_oracle_for_seeded_sequence_pairs():
    randomizer = random.Random(131)
    alphabet = "AAAAACDEFGHIKLMNPQRSTVWYX*"
    for index in range(200):
        query = "".join(randomizer.choices(alphabet, k=randomizer.randrange(5, 60)))
        target = list(query if index % 2 else "".join(randomizer.choices(alphabet, k=randomizer.randrange(5, 60))))
        for _ in range(randomizer.randrange(1, 7)):
            if not target:
                break
            position = randomizer.randrange(len(target))
            target[position:position + randomizer.randrange(0, 5)] = randomizer.choices(alphabet, k=randomizer.randrange(0, 5))
        if target:
            test_biotite_local_alignment_preserves_scores_gaps_and_ties(query, "".join(target))


def test_biotite_alignment_does_not_silently_truncate_optimal_ties():
    align = _load("homology").local_alignment
    with pytest.raises(ValueError, match="complete tie handling"):
        align("AAAA", "AAAAAAAA", max_alignments=1)
    assert align("", "AAAA") is None
    assert align("AAAA", "") is None


def test_biotite_blosum62_matrix_matches_previous_matrix():
    from Bio.Align import substitution_matrices

    homology = _load("homology")
    previous = substitution_matrices.load("BLOSUM62")
    for first in homology.BLOSUM62_ALPHABET:
        for second in homology.BLOSUM62_ALPHABET:
            assert homology._BLOSUM62.get_score(first, second) == previous[first, second]


def test_dog_search_driver_fails_when_a_run_has_a_failure_marker(tmp_path):
    work = tmp_path / "search"
    (work / "db").mkdir(parents=True)
    (work / "work").mkdir()
    (work / "db" / "dog_db.fasta").write_text(">CANON000000\nPEPTIDE\n")
    (work / "runs.tsv").write_text(
        "failed\tfixture\thttps://example.invalid/archive.zip\tfailed.mzML\t1\t1\n")
    (work / "work" / "failed.fetch_failed").write_text("intentional test marker\n")

    sage = tmp_path / "sage"
    sage.write_text("#!/bin/sh\n[ \"${1:-}\" = --version ] && echo 'sage test'\n")
    sage.chmod(0o755)
    python = tmp_path / "python"
    python.write_text(
        f"#!/bin/sh\n"
        "case \"$*\" in\n"
        "  *'import pandas, pyarrow, inflate64'*) exit 0;;\n"
        "esac\n"
        f"exec {sys.executable} \"$@\"\n")
    python.chmod(0o755)
    env = dict(os.environ, SAGE=str(sage), PYTHON=str(python), MIN_FREE_GB="0")
    result = subprocess.run(
        ["bash", str(SCRIPTS / "driver.sh"), str(work)], env=env,
        capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert "INCOMPLETE searched=0 failed=1" in result.stdout
