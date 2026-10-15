"""Pinned annotations PyEnsembl's EnsemblRelease doesn't cover, such as T2T-CHM13 (#72)."""

import gzip
import hashlib
from pathlib import Path

import pytest

from orfdata import install_genome, load_genome, pinned_genomes
from orfdata.cli import run

FIXTURES = Path(__file__).parent / "data/openprot"


def test_the_t2t_annotation_is_pinned():
    t2t, = pinned_genomes(species="human")
    assert (t2t.name, t2t.reference_name, t2t.annotation_name, t2t.annotation_version) == (
        "t2t-chm13", "T2T-CHM13v2.0", "ensembl-rapid", "2022_07")
    assert sorted(t2t.assets) == ["gtf", "proteins", "transcripts"]
    assert all(len(a["sha256"]) == 64 and a["url"].startswith("https://ftp.ensembl.org/") for a in t2t.assets.values())


def slice_pins(tmp_path):
    """Real Ensembl 106 GTF and cDNA slices pinned as an annotation."""
    proteins = tmp_path / "slice.pep.fa.gz"
    proteins.write_bytes(gzip.compress(b">ENSMUSP00000105230.2 pep\nMSSMNPEYDYLFK\n", mtime=0))
    paths = dict(gtf=FIXTURES / "ensembl106.gtf.gz", transcripts=FIXTURES / "ensembl106.cdna.fa.gz",
                 proteins=proteins)
    assets = {role: dict(url=f"https://example.invalid/{path.name}", sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         size=path.stat().st_size) for role, path in paths.items()}
    pins = {"mouse-slice": dict(species="mus_musculus", reference_name="GRCm39", annotation_name="ensembl-slice",
                                annotation_version="106", title="Ensembl 106 slice", assets=assets,
                                provider="Ensembl", terms="none", terms_url="https://example.invalid")}
    return paths, pins


def test_a_pinned_annotation_installs_from_local_copies_and_opens_offline(tmp_path, monkeypatch):
    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    with pytest.raises(ValueError, match="not installed; install it once with install_genome.*orfdata genomes install"):
        load_genome("mouse-slice")
    with pytest.raises(ValueError, match="mismatch"):
        install_genome("mouse-slice", source_paths=dict(paths, gtf=paths["transcripts"]))
    installed = install_genome("mouse-slice", source_paths=paths)
    assert (installed.reference_name, installed.species.latin_name) == ("GRCm39", "mus_musculus")
    opened = load_genome("mouse-slice")
    assert opened.transcript_by_id("ENSMUST00000109601").sequence.startswith(
        installed.transcript_by_id("ENSMUST00000109601").sequence[:20])
    with pytest.raises(ValueError, match="No pinned genome 'mouse-sliced'"):
        load_genome("mouse-sliced")
    # The genome keeps its species and assembly accession through pickling, as multiprocessing and Varcode need.
    import pickle

    copied = pickle.loads(pickle.dumps(opened))
    assert copied.species.latin_name == "mus_musculus" and copied.reference_name == "GRCm39"


def test_a_damaged_index_or_file_says_how_to_reinstall(tmp_path, monkeypatch):
    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    install_genome("mouse-slice", source_paths=paths, cache_root=tmp_path / "root")
    index, = (tmp_path / "root").rglob("*.db")
    index.write_bytes(b"not a database")
    with pytest.raises(ValueError, match="index is damaged .*'orfdata genomes install mouse-slice --force'"):
        load_genome("mouse-slice", cache_root=tmp_path / "root")
    assert install_genome("mouse-slice", force=True, source_paths=paths, cache_root=tmp_path / "root")
    assert load_genome("mouse-slice", cache_root=tmp_path / "root").transcript_by_id("ENSMUST00000109601")


def test_the_proteome_cli_takes_a_pinned_genome_or_a_release(capsys):
    for argv, message in ((["proteome"], "Pass --species and --release, or --genome"),
                          (["proteome", "--genome", "t2t-chm13", "--species", "human"], "not both"),
                          (["proteome", "--genome", "t2t-chm13"], "t2t-chm13 annotation is not installed")):
        with pytest.raises(SystemExit):
            run(argv)
        assert message in capsys.readouterr().err
