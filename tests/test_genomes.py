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
    # The command to run names the same cache.
    with pytest.raises(ValueError, match=r"index is damaged .*'orfdata genomes --cache-root \S+root install "
                                         r"mouse-slice --force'"):
        load_genome("mouse-slice", cache_root=tmp_path / "root")
    assert install_genome("mouse-slice", force=True, source_paths=paths, cache_root=tmp_path / "root")
    assert load_genome("mouse-slice", cache_root=tmp_path / "root").transcript_by_id("ENSMUST00000109601")
    # An empty FASTA index, too, says how to repair it rather than raising EOFError.
    pickled = next(p for p in (tmp_path / "root").rglob("*.pickle") if "cdna" in p.name)
    pickled.write_bytes(b"")
    with pytest.raises(ValueError, match="index is damaged"):
        load_genome("mouse-slice", cache_root=tmp_path / "root").transcript_by_id("ENSMUST00000109601").sequence


def test_a_genome_installed_under_a_relative_cache_works_from_anywhere(tmp_path, monkeypatch):
    import os
    import pickle

    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    monkeypatch.chdir(tmp_path)
    blob = pickle.dumps(install_genome("mouse-slice", source_paths=paths, cache_root="relative"))
    monkeypatch.chdir(os.sep)
    assert pickle.loads(blob).transcript_by_id("ENSMUST00000109601").sequence
    assert len({hash(g) for g in pinned_genomes()}) == 1  # Hashable, as a frozen dataclass should be.


def test_the_proteome_cli_takes_a_pinned_genome_or_a_release(capsys):
    for argv, message in ((["proteome"], "Pass --species and --release, or --genome"),
                          (["proteome", "--genome", "t2t-chm13", "--species", "human"], "not both"),
                          (["proteome", "--genome", "t2t-chm13"], "t2t-chm13 annotation is not installed")):
        with pytest.raises(SystemExit):
            run(argv)
        assert message in capsys.readouterr().err


def test_nonstandard_assembly_keeps_species_genetic_codes_and_saved_overrides():
    import pickle

    from orfdata._pinned_genome import PinnedEnsemblGenome

    genome = PinnedEnsemblGenome("T2T-CHM13v2.0", "ensembl-rapid", "2022_07",
                                species="human", assembly_accession="GCA_009914755.4")
    for copied in (genome, pickle.loads(pickle.dumps(genome)),
                   PinnedEnsemblGenome.from_json(genome.to_json())):
        assert copied.taxon_id == 9606
        assert copied.reference_metadata["translation_tables"] == {
            "nuclear": 1, "mitochondrial": 2, "plastid": None}
        assert copied.assembly_accession == "GCA_009914755.4"
    # Saved provenance must not be silently refreshed from today's catalogue.
    state = genome.to_dict()
    state["reference_metadata"]["translation_tables"]["nuclear"] = 2
    state["reference_metadata"]["translation_table_sources"]["nuclear"] = "reference"
    copied = PinnedEnsemblGenome.from_dict(state)
    assert copied.reference_metadata == state["reference_metadata"]


def test_explicit_chromosome_dna_install_and_offline_load(tmp_path, monkeypatch):
    from pyensembl.genome_fasta import MissingGenomeFastaError

    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    dna = tmp_path / "chromosomes.fa"
    dna.write_text(">MT\nATATGAAGG\n")
    installed = install_genome("mouse-slice", source_paths=paths, genome_fasta_path_or_url=dna)
    assert installed.sequence("MT", 2, 6) == "TATGA"
    assert installed.taxon_id == 10090
    opened = load_genome("mouse-slice", genome_fasta_path_or_url=dna)
    assert opened.sequence("MT", 1, 9) == "ATATGAAGG"
    assert opened.transcript_by_id("ENSMUST00000109601").sequence
    assert not load_genome("mouse-slice").requires_genome_fasta
    missing = load_genome("mouse-slice", genome_fasta_path_or_url=tmp_path / "missing.fa")
    with pytest.raises(MissingGenomeFastaError):
        missing.sequence("MT", 1, 3)
    monkeypatch.setattr(type(opened), "download_genome_fasta", lambda *a, **k: pytest.fail("load downloaded DNA"))
    monkeypatch.setattr(type(opened), "index_genome_fasta", lambda *a, **k: pytest.fail("load opened DNA"))
    assert load_genome("mouse-slice", genome_fasta_path_or_url="https://example.invalid/chromosomes.fa").requires_genome_fasta


def test_genome_cli_passes_explicit_dna_source(tmp_path, monkeypatch, capsys):
    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    dna = tmp_path / "chromosomes.fa"
    dna.write_text(">MT\nATATGAAGG\n")
    argv = ["genomes", "install", "mouse-slice", "--genome-fasta", str(dna)]
    for role, path in paths.items():
        argv += ["--source-path", f"{role}={path}"]
    run(argv)
    assert "reference name: GRCm39" in capsys.readouterr().out
    run(argv + ["--json"])
    assert '"reference_name": "GRCm39"' in capsys.readouterr().out
    assert load_genome("mouse-slice", genome_fasta_path_or_url=dna).sequence("MT", 1, 3) == "ATA"
