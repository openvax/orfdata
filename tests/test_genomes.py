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


# T2T-CHM13v2.0's names, contigs and mitochondria (#204). Checked 2026-11-09 against NCBI's assembly reports
# (tests/data/t2t), UCSC's hs1.chromAlias.txt and the marbl/CHM13 READMEs.
T2T = Path(__file__).parent / "data/t2t"
OWN, RCRS = "CP068254.1", "NC_012920.1"
T2T_NAMES = {
    "t2t-chm13": OWN, "T2T-CHM13v2.0": OWN, "GCA_009914755.4": OWN, "GCF_009914755.1": None, "hs1": OWN,
    "chm13v2.0": OWN, "chm13v2.0_maskedY": OWN, "chm13v2.0_noY": OWN,
    "chm13v2.0_maskedY_rCRS": RCRS, "chm13v2.0_maskedY.rCRS": RCRS,
    "T2T-CHM13v1.1": OWN, "GCA_009914755.3": OWN, "chm13.draft_v1.1": OWN,
}
OLDER_T2T = ("T2T-CHM13v1.0", "GCA_009914755.2", "chm13.draft_v1.0", "chm13.draft_v0.9", "T2T-CHM13v0.7",
             "GCA_009914755.1", "chm13.draft_v0.7")


def assembly_report(path):
    return [line.rstrip("\n").split("\t") for line in path.read_text().splitlines() if not line.startswith("#")]


def test_every_t2t_name_opens_the_pin_with_its_own_mitochondrion():
    from orfdata import reference_names
    from orfdata.genomes import _pinned

    names = reference_names()
    assert list(names.columns) == ["name", "genome", "mitochondrion", "note"]
    found = {row.name: (row.genome, row.mitochondrion) for row in names.itertuples()}
    assert {name: found[name] for name in T2T_NAMES} == {name: ("t2t-chm13", mt) for name, mt in T2T_NAMES.items()}
    assert {name: found[name] for name in OLDER_T2T} == {name: (None, None) for name in OLDER_T2T}
    assert len(names) == len(T2T_NAMES) + len(OLDER_T2T) and names["name"].is_unique
    for name in T2T_NAMES:
        assert _pinned(name).name == _pinned(name.upper()).name == "t2t-chm13"
    # Bare "chm13" names every version, so it selects none.
    with pytest.raises(ValueError, match="No pinned genome 'chm13'.*orfdata genomes names"):
        _pinned("chm13")


def test_older_t2t_versions_say_why_their_coordinates_differ():
    from orfdata import contig_aliases, reference_names
    from orfdata.genomes import _pinned

    v10 = {row[0]: int(row[8]) for row in assembly_report(T2T / "GCA_009914755.2_T2T-CHM13v1.0_assembly_report.txt")}
    v20 = {row[0]: int(row[8]) for row in assembly_report(T2T / "GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt")}
    assert (v10["1"], v20["1"]) == (248387497, 248387328)
    assert sum(v10[name] != v20[name] for name in v10) == 23 and len(v10) == 24
    notes = reference_names().set_index("name")["note"]
    for name in OLDER_T2T:
        for call in (_pinned, contig_aliases, load_genome):
            with pytest.raises(ValueError, match=f"{name} has no pinned genome. ") as error:
                call(name)
            assert str(error.value).endswith(notes[name])
    assert "248,387,497" in notes["T2T-CHM13v1.0"] and "draft" in notes["chm13.draft_v0.9"]


def test_t2t_contigs_match_ncbis_assembly_report_and_ucscs_aliases():
    from orfdata import contig_aliases

    report = T2T / "GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt"
    frame = contig_aliases("t2t-chm13")
    assert frame.attrs["source"]["sha256"] == hashlib.sha256(report.read_bytes()).hexdigest()
    assert frame.attrs["source"]["url"].endswith(report.name) and frame.attrs["genome"] == "t2t-chm13"
    expected = [(row[0], row[9], row[4], None if row[6] == "na" else row[6], int(row[8]))
                for row in assembly_report(report)]
    assert len(expected) == 25
    assert list(frame.itertuples(index=False, name=None)) == expected
    ucsc = {row[0]: tuple(row[1:4]) for row in (line.split("\t") for line in
                                                (T2T / "hs1.chromAlias.txt").read_text().splitlines()[1:])}
    assert {row.ucsc: (row.genbank, row.refseq or "", row.name) for row in frame.itertuples()} == ucsc
    # The chromosome a GRCh38 file would be caught by, and the mitochondrion header lengths can't tell apart.
    lengths = dict(zip(frame["ucsc"], frame["length"]))
    assert lengths["chr1"] == 248387328 and lengths["chrM"] == 16569


def test_each_t2t_name_lists_only_the_contigs_its_files_share_with_the_pin():
    from orfdata import contig_aliases

    everything = set(contig_aliases("t2t-chm13")["name"])
    for name in T2T_NAMES:
        contigs = set(contig_aliases(name)["name"])
        missing = everything - contigs
        if name in ("chm13v2.0_noY", "T2T-CHM13v1.1", "GCA_009914755.3", "chm13.draft_v1.1"):
            assert missing == {"Y"}, name
        elif T2T_NAMES[name] != OWN:
            # The _rCRS set's chrM is the rCRS, a different molecule; RefSeq's assembly has no mitochondrion.
            assert missing == {"MT"}, name
        else:
            assert not missing, name
    # The _rCRS set's chrM is the rCRS, which the genome annotates as NC_012920.1, never as MT.
    rcrs = contig_aliases("chm13v2.0_maskedY_rCRS").set_index("ucsc")["name"]
    assert rcrs["chrM"] == "NC_012920.1" and set(rcrs.index) >= {"chr1", "chrX", "chrY"}


def test_a_t2t_name_that_is_not_installed_names_the_pins_install_command(tmp_path):
    from orfdata import NotInstalledError

    with pytest.raises(NotInstalledError) as error:
        load_genome("chm13v2.0_maskedY_rCRS", cache_root=tmp_path)
    assert str(error.value).startswith("The t2t-chm13 annotation (for chm13v2.0_maskedY_rCRS) is not installed")
    assert f"'orfdata genomes --cache-root {tmp_path} install t2t-chm13'" in str(error.value)
    assert isinstance(error.value, FileNotFoundError) and isinstance(error.value, ValueError)


def test_an_installed_genome_opens_under_each_of_its_names(tmp_path, monkeypatch):
    from orfdata import contig_aliases, genomes, reference_names

    paths, pins = slice_pins(tmp_path)
    pins["mouse-slice"].update(
        names={"mm39-slice": dict(mitochondrion="AY172335.1", note="A test name")},
        other_versions={"mm10-slice": "A different assembly."},
        contigs=[dict(name="MT", ucsc="chrM", genbank="AY172335.1", refseq="NC_005089.1", length=16299)])
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    install_genome("mm39-slice", source_paths=paths)
    for name in ("mouse-slice", "GRCm39", "MM39-SLICE"):
        assert load_genome(name).transcript_by_id("ENSMUST00000109601").sequence
    assert list(reference_names()["name"]) == ["mouse-slice", "GRCm39", "mm39-slice", "mm10-slice"]
    assert contig_aliases("mm39-slice")["ucsc"].tolist() == ["chrM"]
    with pytest.raises(ValueError, match="mm10-slice has no pinned genome. A different assembly."):
        load_genome("mm10-slice")


def test_catalogue_hints_use_the_pinned_names():
    from orfdata import load_catalog

    with pytest.raises(ValueError, match="orfdata names that assembly T2T-CHM13v2.0"):
        load_catalog(species="human", genome_build="hs1")
    with pytest.raises(ValueError, match="chm13.draft_v1.0 has no pinned genome. T2T-CHM13v1.0's coordinates"):
        load_catalog(species="human", genome_build="chm13.draft_v1.0")
    with pytest.raises(ValueError) as error:
        load_catalog(species="human", genome_build="chm13")
    assert "T2T-CHM13v2.0, as Ensembl does" not in str(error.value)


def test_genome_names_and_contigs_commands(capsys):
    run(["genomes", "names"])
    output = capsys.readouterr().out
    assert "chm13v2.0_maskedY_rCRS  t2t-chm13  NC_012920.1" in output and "T2T-CHM13v1.0 " in output
    run(["genomes", "names", "--json"])
    import json

    rows = {row["name"]: row for row in json.loads(capsys.readouterr().out)}
    assert rows["GCF_009914755.1"]["mitochondrion"] is None and rows["hs1"]["genome"] == "t2t-chm13"
    run(["genomes", "contigs", "hs1", "--json"])
    contigs = json.loads(capsys.readouterr().out)
    assert len(contigs) == 25 and contigs[0] == dict(name="1", ucsc="chr1", genbank="CP068277.2",
                                                     refseq="NC_060925.1", length=248387328)
    with pytest.raises(SystemExit):
        run(["genomes", "contigs", "T2T-CHM13v1.0"])
    assert "has no pinned genome" in capsys.readouterr().err


# The rCRS mitochondrion t2t-chm13 carries as NC_012920.1 (openvax/varcode#722).
RCRS = "NC_012920.1"


def rcrs_sequence():
    text = (T2T / "NC_012920.1.fasta").read_text()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def packaged_rcrs():
    """The packaged GTF rows, transcripts and proteins, read only when they are the pinned bytes."""
    from orfdata.reference_table import packaged_bytes

    t2t, = pinned_genomes(species="human")
    assets = t2t.mitochondria[RCRS]["assets"]
    texts = {role: gzip.decompress(packaged_bytes(a["resource"], a["sha256"])).decode() for role, a in assets.items()}

    def fasta(text, key=lambda header: header.split()[0]):
        records = {}
        for block in text.split(">")[1:]:
            header, *lines = block.splitlines()
            records[key(header).split(".")[0]] = "".join(lines)
        return records

    rows = [line.split("\t") for line in texts["gtf"].splitlines() if not line.startswith("#")]
    # Proteins by the transcript their header names.
    return rows, fasta(texts["transcripts"]), fasta(texts["proteins"], lambda h: h.split("transcript:")[1].split()[0])


def test_the_packaged_rcrs_annotation_spells_ncbis_rcrs():
    import re

    from orfdata.sequence import genetic_code
    from orfdata.viral_annotation import extract, translate

    dna = rcrs_sequence()
    assert len(dna) == 16569 and dna.startswith("GATCACAGGTCTATCACCC")
    rows, transcripts, proteins = packaged_rcrs()
    assert {row[0] for row in rows} == {RCRS}
    exons, cds = {}, {}
    for row in rows:
        transcript = re.search(r'transcript_id "([^"]+)"', row[8])
        if row[2] == "exon":
            exons.setdefault(transcript[1], []).append((int(row[3]), int(row[4]), row[6]))
        elif row[2] == "CDS":
            cds.setdefault(transcript[1], []).append((int(row[3]), int(row[4]), row[6]))
    genes = {re.search(r'gene_name "([^"]+)"', row[8])[1] for row in rows if row[2] == "gene"}
    assert len(genes) == 37 and {"MT-ND1", "MT-CO1", "MT-ATP6", "MT-CYB", "MT-TF", "MT-RNR2"} <= genes
    assert set(exons) == set(transcripts) and len(cds) == len(proteins) == 13
    for transcript, spans in exons.items():
        assert extract(dna, sorted(spans, reverse=spans[0][2] == "-")) == transcripts[transcript], transcript
    starts = genetic_code(2).start_codons
    for transcript, spans in cds.items():
        coding = extract(dna, sorted(spans, reverse=spans[0][2] == "-"))
        protein = translate(coding, table=2)
        assert ("M" + protein[1:] if coding[:3] in starts else protein) == proteins[transcript], transcript
        # The vertebrate mitochondrial code matters: TGA is tryptophan, not a stop.
        assert translate(coding, table=1) != protein


def test_the_rcrs_files_contig_table_names_their_chrm_by_its_accession():
    from orfdata import contig_aliases

    for name in ("chm13v2.0_maskedY_rCRS", "chm13v2.0_maskedY.rCRS"):
        frame = contig_aliases(name)
        row = frame[frame["ucsc"] == "chrM"].iloc[0].to_dict()
        assert row == dict(name=RCRS, ucsc="chrM", genbank="J01415.2", refseq=RCRS, length=16569)
        assert "MT" not in set(frame["name"]) and len(frame) == 25
    # Files with CHM13's own mitochondrion keep MT, CP068254.1; none lists the rCRS.
    for name in ("t2t-chm13", "hs1", "chm13v2.0", "GCF_009914755.1"):
        assert RCRS not in set(contig_aliases(name)["name"]), name
    assert contig_aliases("hs1").set_index("ucsc").loc["chrM", "name"] == "MT"


def test_a_genome_with_a_packaged_mitochondrion_annotates_it_on_its_own_contig(tmp_path, monkeypatch):
    import pickle

    from orfdata import genomes

    paths, pins = slice_pins(tmp_path)
    t2t, = pinned_genomes(species="human")
    pins["mouse-slice"]["mitochondria"] = {RCRS: dict(t2t.mitochondria[RCRS])}
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    root = tmp_path / "root"
    installed = install_genome("mouse-slice", source_paths=paths, cache_root=root)
    for genome in (installed, load_genome("mouse-slice", cache_root=root),
                   pickle.loads(pickle.dumps(installed))):
        assert RCRS in genome.contigs() and "17" in genome.contigs()
        assert genome.transcript_by_id("ENSMUST00000109601").sequence  # The slice is still there.
        nd1, = genome.transcripts_at_locus(RCRS, 3307, 3307)
        assert (nd1.gene_name, nd1.biotype, nd1.strand) == ("MT-ND1", "protein_coding", "+")
        assert nd1.translation_table == 2 and nd1.protein_sequence.startswith("MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEP")
        assert genome.genes_by_name("MT-TF")[0].start == 577
    # The packaged files are copied beside the index, so the genome doesn't depend on where ORFdata is installed.
    assert {p.name for p in root.rglob("rcrs-ensembl115.*")} == {"rcrs-ensembl115.cdna.fa.gz",
                                                                 "rcrs-ensembl115.pep.fa.gz"}
    merged, = root.rglob(f"*+{RCRS}.gtf.gz")
    assert merged.stat().st_size == paths["gtf"].stat().st_size + t2t.mitochondria[RCRS]["assets"]["gtf"]["size"]
    # A forced install rebuilds the index and the merged GTF.
    assert install_genome("mouse-slice", source_paths=paths, cache_root=root, force=True).genes_by_name("MT-CO1")


def test_adding_a_mitochondrion_reindexes_without_downloading_again(tmp_path, monkeypatch):
    from orfdata import genomes

    t2t, = pinned_genomes(species="human")
    paths, pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: pins)
    root = tmp_path / "root"
    install_genome("mouse-slice", source_paths=paths, cache_root=root)
    assert not load_genome("mouse-slice", cache_root=root).transcripts_at_locus(RCRS, 3307, 3307)
    pins["mouse-slice"]["mitochondria"] = {RCRS: dict(t2t.mitochondria[RCRS])}
    # The downloaded files' store is unchanged, so loading builds the new index from them, offline.
    opened = load_genome("mouse-slice", cache_root=root)
    assert opened.transcripts_at_locus(RCRS, 3307, 3307)[0].gene_name == "MT-ND1"
