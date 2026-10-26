import importlib.util
import os
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
