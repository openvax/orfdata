"""Dog peptides: DLA ligands (PXD074485) and the Canine PeptideAtlas, from real source slices."""

import csv
import hashlib
import io
import json
import re
import shutil
import zipfile
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import DatasetManifest, DatasetRegistry, ReferenceContext, ScopedTarget, build_dataset
from orfdata.cli import run
from orfdata.sources.kaabinejadian import import_ligands, irt_rules, irt_standard, peptide_sequence
from orfdata.sources.peptideatlas import import_build

DLA = Path(__file__).parent / "data/kaabinejadian"
ATLAS = Path(__file__).parent / "data/peptideatlas"


def _source(root):
    manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    return manifest, {asset.asset_id: root / asset.asset_id for asset in manifest.sources}


def _build(root, importer, tmp_path, manifest=None):
    found, inputs = _source(root)
    manifest = manifest or found
    return build_dataset(manifest, lambda i, o: importer(i, o, manifest), source_paths=inputs,
                         cache_root=tmp_path, download_options={"show_progress": False})


def _tables():
    with zipfile.ZipFile(DLA / "Results.zip") as archive:
        for name in archive.namelist():
            rows = list(csv.reader(io.TextIOWrapper(archive.open(name), encoding="utf-8-sig", newline="")))
            yield name[len("Results/"):-len(".csv")], rows[0], list(enumerate(rows[1:], 2))


def test_packaged_dog_peptide_releases_are_discoverable_offline(monkeypatch):
    import requests

    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("Discovery contacted network"))
    registry = DatasetRegistry()
    ligands = registry.info("kaabinejadian-dog-dla-ligands")
    atlas = registry.info("peptideatlas-dog-peptides")
    assert ligands in registry.list(species="dog", assay="mhc_ms")
    assert ligands in registry.list(species="human")  # Its monoallelic lines are human cells.
    assert atlas in registry.list(species="dog", assay="digest_ms") and atlas not in registry.list(species="human")
    assert registry.published("kaabinejadian-dog-dla-ligands")
    assert not registry.published("peptideatlas-dog-peptides")
    assert registry.download_metadata("peptideatlas-dog-peptides")["distribution"]["redistributable"] is False
    with pytest.raises(ValueError, match="build it from its pinned sources"):
        registry.download("peptideatlas-dog-peptides")
    # The Ensembl context equals PyEnsembl's dog release 110, so mapped queries need no conversion.
    _, ensembl = atlas.references
    assert ensembl.key == ReferenceContext(9615, "canis_familiaris", "ROS_Cfam_1.0", "ensembl", "110").key


def test_peaks_rows_keep_values_identification_type_and_modified_form(tmp_path):
    dataset = _build(DLA, import_ligands, tmp_path)
    manifest = dataset.manifest
    values = {m.measurement_id: m for m in dataset.measurements()}
    rows = 0
    for member, header, lines in _tables():
        for line, row in lines:
            rows += 1
            area, intensity, spectra = (values[f"{member}:{line}:{column}"] for column in ("Area", "Intensity", "#Spec"))
            assert area.target.kind == "peptide" and area.target.namespace == "amino-acid-sequence"
            assert area.target.identifier == re.sub(r"\([^)]*\)", "", row[0])
            assert area.target == intensity.target == spectra.target
            assert area.value == (float(row[14]) if row[14] else None)
            assert area.status == ("measured" if row[14] else "missing")
            assert (area.unit, intensity.unit, spectra.unit) == ("PEAKS area", "PEAKS intensity", "PSM")
            assert spectra.quantity == "spectrum_count" and spectra.value == int(row[17])
            assert spectra.provenance["peptide"] == row[0] and spectra.provenance["found_by"] == row[25]
            rule = spectra.provenance["irt_rule"]
            assert rule in (None, "irt_sequence", "kit_only", "irt_precursor", "de_novo_irt")
            assert spectra.provenance["irt_standard"] is (rule is not None)
            assert spectra.provenance["raw_file_deposited"] is True
            assert spectra.provenance["column"] == header[17]
    assert len(values) == rows * 3 == 99
    dog, human = manifest.references
    assert (dog.taxon_id, human.taxon_id) == (9615, 9606) and not dog.is_genomic
    # The same peptide in two searches of Lola's spectra is one target, measured in two assays.
    shared = ScopedTarget(dog, "peptide", "AAARLTH", "amino-acid-sequence")
    assert {m.assay_id for m in dataset.measurements(target=shared, quantity="spectrum_count")} == {
        "Lola-BB76:UniProt9615", "Lola-BB76:MergedCanFams"}
    assert [m.sample_id for m in dataset.measurements(target=ScopedTarget(
        human, "peptide", "AAARLTH", "amino-acid-sequence"))] == ["HCT116-DLA8850101"] * 3


def test_samples_and_assays_keep_dla_typing_antibody_and_search(tmp_path):
    dataset = _build(DLA, import_ligands, tmp_path)
    samples = {s.sample_id: s for s in dataset.samples()}
    lola, cells = samples["Lola"], samples["HCT116-DLA8850101"]
    assert lola.reference.taxon_id == 9615 and cells.reference.taxon_id == 9606
    assert (lola.tissue, lola.condition) == ("bone tissue", "osteosarcoma")
    assert lola.attributes["breed"] == "Black Russian Terrier"
    assert lola.attributes["mhc_typing"] == ("DLA-88*501:01",)
    assert cells.attributes["mhc_typing"] == ("DLA-88*501:01",) and cells.cell_type == "HCT 116"
    assert cells.attributes["mhc_typing_source"] == "deposit SDRF; the article's transduced allele"
    assays = {a.assay_id: a for a in dataset.assays()}
    assert set(assays) == {"Lola-BB76:UniProt9615", "Lola-BB76:MergedCanFams", "Lola-H58A:UniProt9615",
                           "HCT116-DLA8850101-BB76:SwissProt9606", "Tumor163828A-BB76:UniProt9615",
                           "Tumor163828A-BB76:MergedCanFams"}
    assert all(a.kind == "mhc_ms" for a in assays.values())
    # One tumour, two immunoprecipitations: one sample, with an assay per antibody and search.
    assert {a.sample_id for a in assays.values() if a.sample_id == "Lola"} == {"Lola"}
    assert assays["Lola-BB76:UniProt9615"].attributes["antibody"] == "BB7.6 (anti-MHC class I)"
    assert assays["Lola-H58A:UniProt9615"].attributes["antibody"] == "H58A (anti-MHC class I)"
    assert assays["Lola-BB76:MergedCanFams"].run_id == "20250310-LolaTumor-BB76-022725"
    assert assays["Lola-H58A:UniProt9615"].run_id == "23-0425-Lola Canine Tumor-H58A-040423"


def test_accessions_are_kept_per_search_with_standards_and_failures_explicit(tmp_path):
    dataset = _build(DLA, import_ligands, tmp_path)
    mappings = {m.mapping_id: m for m in dataset.mappings()}
    unique = mappings["MergedCanFams:DB Search:AERHGYIKGIVKDIIHD"]
    assert unique.status == "unique" and unique.source_namespace == "PXD074485:MergedCanFams:DB Search"
    assert [(t.kind, t.identifier, t.namespace) for t in unique.targets] == [
        ("protein", "ENSCAFP00805008970.1", "PXD074485:MergedCanFams")]
    assert unique.method == "PEAKS-reported proteins (DB Search)"
    ambiguous = mappings["UniProt9615:DB Search:AIKKGLEEL"]
    assert ambiguous.status == "ambiguous" and {t.identifier for t in ambiguous.targets} == {"F1PFM6", "J9P889"}
    homolog = mappings["MergedCanFams:Homolog:GKVVSETNDTKVLRH"]
    assert homolog.status == "unique" and homolog.method == "PEAKS-reported proteins (Homolog)"
    # The kit peptide and its acetylated and oxidized forms share one sequence.
    standard = mappings["UniProt9615:DB Search:ADVTPADFSEWSK"]
    assert standard.status == "unmapped" and standard.reason == "Biognosys iRT retention-time standard"
    assert standard.evidence["irt_standards"] == ("Biognosys|iRT-Kit_peptide_10",)
    # A homolog of a kit peptide, which PEAKS attributes only to the kit.
    assert mappings["UniProt9615:Homolog:ADTTPADFSEWSK"].reason == "Biognosys iRT retention-time standard"
    # The merged database lacks the kit, so PEAKS names dog proteins for a kit fragment.
    fragment = mappings["MergedCanFams:DeepNovo:AASYYAPVR"]
    assert fragment.reason == "Biognosys iRT retention-time standard" and fragment.targets == ()
    assert fragment.evidence["reported_proteins"]
    # A de novo read that PEAKS matches to the kit among dog proteins.
    novo = mappings["UniProt9615:DeepNovo:AAFSKDPAAVLR"]
    assert novo.status == "unmapped" and novo.reason == "De novo read of a Biognosys iRT standard"
    assert "J9P3Z7" in novo.evidence["reported_proteins"] and novo.evidence["irt_standards"]
    assert novo.evidence["irt_rule"] == "de_novo_irt"
    assert mappings["UniProt9615:DeepNovo:ADVTMEAPDAVESK"].reason.startswith("De novo sequence")
    # Homology and kit-less database-search reads of iRT spectra lose their dog proteins.
    for key in ("UniProt9615:Homolog:TGFLLDPVLNR", "MergedCanFams:DB Search:IFLAFQQKAYVAP"):
        misread = mappings[key]
        assert misread.status == "unmapped" and misread.evidence["irt_rule"] == "irt_precursor"
        assert misread.reason.startswith("Observed at a Biognosys iRT standard's precursor")
        assert misread.evidence["reported_proteins"]
    assert mappings["MergedCanFams:DB Search:TNAENEFVTIK"].status == "ambiguous"
    assert "irt_rule" not in mappings["UniProt9615:DB Search:TNAENEFVTIK"].evidence
    assert mappings["UniProt9615:DeepNovo:AAARLTH"].reason == "PEAKS reported no protein"
    human = mappings["SwissProt9606:DB Search:AAALGIAFVEV"]
    assert human.reference.taxon_id == 9606 and [t.identifier for t in human.targets] == ["P08962"]
    assert {m.evidence["found_by"] for m in mappings.values()} == {"DB Search", "DeepNovo", "Homolog"}
    assert all(m.status == "unmapped" for m in mappings.values() if m.evidence["found_by"] == "DeepNovo")


def test_irt_standards_and_misread_irt_spectra_are_flagged(tmp_path):
    dataset = _build(DLA, import_ligands, tmp_path)
    dog = dataset.manifest.references[0]

    def rules(sequence):
        target = ScopedTarget(dog, "peptide", sequence, "amino-acid-sequence")
        return {m.assay_id: m.provenance["irt_rule"] for m in dataset.measurements(target=target, quantity="spectrum_count")}

    assert rules("ADVTPADFSEWSK") == {"Lola-BB76:UniProt9615": "irt_sequence"}
    assert rules("AASYYAPVR") == {"Lola-BB76:MergedCanFams": "irt_sequence"}  # A kit fragment, read de novo.
    assert rules("ADTTPADFSEWSK") == {"Lola-BB76:UniProt9615": "kit_only"}
    # De novo reads: matched to the kit, or sharing seven residues (DPAAVLR) with a kit peptide.
    assert rules("AAFSKDPAAVLR") == {"Lola-BB76:UniProt9615": "de_novo_irt", "Lola-BB76:MergedCanFams": "de_novo_irt"}
    # Reads at the precursor of a kit peptide that a search of the same immunoprecipitation found.
    assert rules("EGVTPANFSQ") == {"Lola-BB76:MergedCanFams": "irt_precursor"}  # ADVTPADFSE, de novo
    assert rules("TGFLLDPVLNR") == {"Tumor163828A-BB76:UniProt9615": "irt_precursor"}  # GTFIIDPGGVIR's mass
    assert rules("IFLAFQQKAYVAP") == {"Tumor163828A-BB76:MergedCanFams": "irt_precursor"}  # No kit to choose.
    # Database search that could choose the kit peptide keeps its match, and so does the kit-less
    # search's match when a search with the kit found the same peptide.
    assert rules("AAAGQAKDEL") == {"Tumor163828A-BB76:UniProt9615": None}
    assert rules("TNAENEFVTIK") == {"Tumor163828A-BB76:UniProt9615": None, "Tumor163828A-BB76:MergedCanFams": None}
    assert rules("AIKKGLEEL") == {"Lola-BB76:UniProt9615": None}
    assert irt_standard("GTFIIDPAAVIR") and irt_standard("IIDPAAVIR") and not irt_standard("SIINFEKL")
    assert irt_standard("ADTTPADFSEWSK", ["Biognosys|iRT-Kit_peptide_10"])


def test_irt_rules_hold_their_tolerances_and_scope():
    kit = ["UniProt9615"]
    anchor = ("A", "UniProt9615", "DB Search", "GTFIIDPGGVIR", ["Biognosys|iRT-Kit_peptide_3"], 2, 622.8533, 54.0)

    def rule(row, *rows):
        return irt_rules([anchor, *rows, row], kit, 20, 1.0, 6)[-1]

    novo = ("A", "MergedCanFams", "DeepNovo", "SIINFEKL", [], 2, 622.8533, 54.0)
    assert rule(novo) == "irt_precursor"
    assert rule(novo[:6] + (622.8533 * (1 + 21e-6), 54.0)) is None  # Outside 20 ppm.
    assert rule(novo[:7] + (55.1,)) is None  # More than a minute apart.
    assert rule(novo[:5] + (3,) + novo[6:]) is None  # Another charge.
    assert rule(("B",) + novo[1:]) is None  # Another immunoprecipitation.
    assert rule(("A", "UniProt9615", "DB Search", "SIINFEKL", ["P1"], 2, 622.8533, 54.0)) is None  # Kit search.
    merged = ("A", "MergedCanFams", "DB Search", "SIINFEKL", ["P1"], 2, 622.8533, 54.0)
    assert rule(merged) == "irt_precursor"
    assert rule(merged, ("A", "UniProt9615", "DB Search", "SIINFEKL", ["P1"], 2, 500.0, 30.0)) is None
    assert rule(("A", "MergedCanFams", "DeepNovo", "GTFLLDPGG", [], 2, 400.0, 20.0)) == "irt_sequence"  # I = L.
    assert rule(("A", "MergedCanFams", "DeepNovo", "KATLFLQFGA", [], 2, 400.0, 20.0)) == "de_novo_irt"
    assert rule(("A", "MergedCanFams", "DeepNovo", "KATLFLQFG", [], 2, 400.0, 20.0)) == "de_novo_irt"
    assert rule(("A", "MergedCanFams", "DeepNovo", "KATLFLQF", [], 2, 400.0, 20.0)) is None  # Five shared.
    assert rule(("A", "MergedCanFams", "DeepNovo", "KATSIINFEKL", ["Biognosys|x", "P1"], 2, 400.0, 20.0)) == "de_novo_irt"
    assert rule(("A", "MergedCanFams", "Homolog", "KATLFLQFGAQ", ["P1"], 2, 400.0, 20.0)) is None


def test_the_first_release_rebuilds_with_its_own_flags(tmp_path):
    manifest = DatasetManifest.from_dict(json.loads((DLA / "manifest-v1.json").read_text()))
    assert manifest.processing["importer"] == "kaabinejadian-dla-v1"
    dataset = _build(DLA, import_ligands, tmp_path, manifest)
    flags = {m.measurement_id: m.provenance for m in dataset.measurements(quantity="spectrum_count")}
    for member, _, lines in _tables():
        for line, row in lines:
            provenance = flags[f"{member}:{line}:#Spec"]
            listed = [value for value in row[24].split(":") if value]
            assert provenance["irt_standard"] is irt_standard(peptide_sequence(row[0]), listed)
            assert "irt_rule" not in provenance
    mappings = {m.mapping_id: m for m in dataset.mappings()}
    assert mappings["UniProt9615:Homolog:TGFLLDPVLNR"].status == "ambiguous"
    assert mappings["UniProt9615:DeepNovo:AAFSKDPAAVLR"].reason.startswith("De novo sequence")
    assert "irt_rule" not in mappings["MergedCanFams:DeepNovo:AASYYAPVR"].evidence


def _copy_with_sdrf(tmp_path, change):
    root = tmp_path / "source"
    shutil.copytree(DLA, root)
    (root / "sdrf.tsv").write_text(change((root / "sdrf.tsv").read_text()))
    manifest = json.loads((root / "manifest.json").read_text())
    data = (root / "sdrf.tsv").read_bytes()
    for asset in manifest["sources"]:
        if asset["asset_id"] == "sdrf.tsv":
            asset.update(sha256=hashlib.sha256(data).hexdigest(), size=len(data))
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root


def test_one_sample_must_have_the_same_sdrf_facts_in_every_assay(tmp_path):
    def older(text):  # Make Lola five in one row and six in the other.
        lines = text.split("\n")
        index = next(i for i, line in enumerate(lines) if line.startswith("Lola\t") and "H58A" in line)
        lines[index] = lines[index].replace("\t5Y\t", "\t6Y\t")
        return "\n".join(lines)
    root = _copy_with_sdrf(tmp_path, older)
    with pytest.raises(ValueError, match="disagree between assays"):
        _build(root, import_ligands, tmp_path / "cache")


def test_rows_from_undeposited_runs_are_marked(tmp_path):
    manifest, _ = _source(DLA)
    tables = manifest.to_dict()["processing"]["tables"]
    tables["Lola-H58A"]["runs"] = {run: False for run in tables["Lola-H58A"]["runs"]}
    dataset = _build(DLA, import_ligands, tmp_path / "marked",
                     replace(manifest, processing=dict(manifest.processing, tables=tables)))
    deposited = {m.assay_id: m.provenance["raw_file_deposited"] for m in dataset.measurements(quantity="spectrum_count")}
    assert deposited["Lola-H58A:UniProt9615"] is False and deposited["Lola-BB76:UniProt9615"] is True
    tables["Lola-H58A"]["runs"] = {}
    with pytest.raises(ValueError, match="undeclared run"):
        _build(DLA, import_ligands, tmp_path / "undeclared",
               replace(manifest, processing=dict(manifest.processing, tables=tables)))


def test_ligand_import_reconciles_against_the_deposit(tmp_path):
    manifest, _ = _source(DLA)
    counts = manifest.to_dict()["coverage"]["source_counts"]
    counts["Lola-BB76-UniProt9615"]["irt_standard_rows"] += 1
    changed = replace(manifest, coverage=dict(manifest.coverage, source_counts=counts))
    with pytest.raises(ValueError, match="reconciliation failed"):
        _build(DLA, import_ligands, tmp_path / "counts", changed)
    counts = manifest.to_dict()["coverage"]["source_counts"]
    counts["Lola-BB76-UniProt9615"]["irt_rules"] = {"irt_sequence": 4, "kit_only": 1}
    changed = replace(manifest, coverage=dict(manifest.coverage, source_counts=counts))
    with pytest.raises(ValueError, match="reconciliation failed"):
        _build(DLA, import_ligands, tmp_path / "rules", changed)
    tables = manifest.to_dict()["processing"]["tables"]
    tables["Lola-BB76"]["antibody"] = "H58A (anti-MHC class I)"
    changed = replace(manifest, processing=dict(manifest.processing, tables=tables))
    with pytest.raises(ValueError, match="SDRF"):
        _build(DLA, import_ligands, tmp_path / "sdrf", changed)


@pytest.mark.parametrize("peptide, sequence", [("A(+42.01)DVTPADFSEWSK", "ADVTPADFSEWSK"),
                                               ("ADVTPAD(+21.98)FSEW(+15.99)SK", "ADVTPADFSEWSK"),
                                               ("Q(-17.03)LLEEL", "QLLEEL")])
def test_peaks_modifications_are_stripped(peptide, sequence):
    assert peptide_sequence(peptide) == sequence
    with pytest.raises(ValueError):
        peptide_sequence("K.PEPTIDE.R")


def test_atlas_keeps_psms_per_sample_and_protein_calls(tmp_path):
    dataset = _build(ATLAS, import_build, tmp_path)
    source, ensembl = dataset.manifest.references
    assert source.genome_build == "ROS_Cfam_1.0" and source.annotation_name == "PeptideAtlas biosequence set"
    samples = {s.sample_id: s for s in dataset.samples()}
    assert {s: samples[s].tissue for s in samples if s != "build-573"} == {
        "PAe009170": "brain", "PAe009174": "liver", "PAe009533": "plasma"}
    assert samples["PAe009533"].attributes["preparation"] == "depleted timsTOF"
    assert samples["PAe009533"].attributes["breed_source"] == "the article"
    assert samples["PAe009170"].attributes["breed"] == "Labrador Retriever"
    assert samples["build-573"].pool_members == ("PAe009170", "PAe009174", "PAe009533")
    peptide = ScopedTarget(source, "peptide", "AAHEEICTTNEGVMYR", "amino-acid-sequence")
    counts = {(m.sample_id, m.unit): m.value for m in dataset.measurements(target=peptide)}
    assert counts == {("PAe009174", "PSM"): 17, ("PAe009533", "PSM"): 675}  # Liver and timsTOF plasma.
    calls = {m.target.identifier: m for m in dataset.measurements(sample_id="build-573", reference=source)}
    assert len(calls) == 7 and calls["ENSCAFP00845001027.1"].target.namespace == "ensembl"
    assert calls["ENSCAFP00845001027.1"].provenance["presence_level"] == "noncore-canonical"
    assert calls["DECOY_A0A1K0GV43"].target.namespace == "PeptideAtlas:decoy"
    assert calls["CONTAM_sp|Q9BQ66|KR412_HUMAN"].target.namespace == "PeptideAtlas:contaminant"
    assert dataset.counts["measurement"] == 4 + 7


def test_atlas_peptides_keep_every_location_and_flag_decoys(tmp_path):
    dataset = _build(ATLAS, import_build, tmp_path)
    mappings = {m.mapping_id: m for m in dataset.mappings()}
    assert len(mappings) == 6
    shared = mappings["PAp00000001"]
    assert shared.source_identifier == "AAHEEICTTNEGVMYR" and shared.status == "ambiguous"
    assert {t.namespace for t in shared.targets} == {"uniprot", "ensembl", "PeptideAtlas:contaminant"}
    exons = [place for place in shared.evidence["locations"] if place[0] == "ENSCAFP00845030549.1"]
    # An exon-spanning peptide on the minus strand: its first residues lie at the higher coordinates.
    assert exons == [("ENSCAFP00845030549.1", 467, 482, "37", "-", 22534513, 22534516, True),
                     ("ENSCAFP00845030549.1", 467, 482, "37", "-", 22533401, 22533444, True)]
    assert shared.evidence["exon_spanning"] is True and shared.evidence["psms"] == 768
    decoy = mappings["PAp16560428"]
    assert decoy.status == "unmapped" and decoy.reason == "Maps only to decoy sequences"
    assert decoy.evidence["decoy_locations"] >= 1
    assert any(place[7] is False for place in mappings["PAp07488729"].evidence["locations"])  # I/L-equivalent.
    assert mappings["PAp16545961"].status == "unique"


def test_atlas_protein_calls_reach_ensembl_110_by_identical_sequence(tmp_path):
    dataset = _build(ATLAS, import_build, tmp_path)
    _, ensembl = dataset.manifest.references
    mappings = {m.source.identifier: m for m in dataset.reference_mappings()}
    assert set(mappings) == {"A0A0A0RBT6", "A0A5F4CZ72", "Q867B0", "ENSCAFP00845001027.1", "A0A8I3PHG0"}
    assert mappings["Q867B0"].status == "unmapped"
    assert mappings["A0A5F4CZ72"].status == "ambiguous"
    native = mappings["ENSCAFP00845001027.1"]
    assert native.status == "unique" and native.sequence_validated
    target = native.targets[0]
    assert (target.reference, target.identifier, target.namespace) == (ensembl, "ENSCAFP00845001027", "ensembl")
    mapped = list(dataset.mapped_measurements(target, require_sequence=True))
    assert [(m.measurement.value, m.measurement.target.identifier) for m in mapped] == [
        (3124, "ENSCAFP00845001027.1")]


def _build_only(manifest):
    return {"peptideatlas-dog-peptides": dict(default_version="fixture", versions={"fixture": dict(
        manifest=manifest.to_dict(), assets={},
        distribution=dict(access_terms="Build locally; not redistributed", redistributable=False))})}


def test_build_only_release_builds_opens_and_inspects_locally(tmp_path, capsys):
    manifest, inputs = _source(ATLAS)
    definitions = _build_only(manifest)
    registry = DatasetRegistry(definitions, cache_root=tmp_path)
    assert registry.inspect("peptideatlas-dog-peptides").status == "missing"
    with pytest.raises(FileNotFoundError, match="hasn't been built.*orfdata datasets build peptideatlas-dog-peptides"):
        registry.open("peptideatlas-dog-peptides")
    config = tmp_path / "registry.json"
    config.write_text(json.dumps(definitions))
    cli = ["datasets", "--registry", str(config), "--cache-root", str(tmp_path)]
    with pytest.raises(SystemExit):
        run(cli + ["download", "peptideatlas-dog-peptides"] + ["--json"])
    assert "build it from its pinned sources" in capsys.readouterr().err
    built = registry.build("peptideatlas-dog-peptides", source_paths=inputs,
                           download_options={"show_progress": False})
    assert registry.inspect("peptideatlas-dog-peptides").status == "available"
    assert registry.open("peptideatlas-dog-peptides").path == built.path
    run(cli + ["inspect", "peptideatlas-dog-peptides"] + ["--json"])
    inspected = json.loads(capsys.readouterr().out)
    assert inspected["status"] == "available" and inspected["origin"] == "source_build"
    source = manifest.references[0]
    run(cli + ["query", "peptideatlas-dog-peptides", "--reference", source.key, "--kind", "peptide",
               "--namespace", "amino-acid-sequence", "--target", "AAHEEICTTNEGVMYR"] + ["--format", "jsonl"])
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert sorted(row["value"] for row in rows) == [17, 675]
    unknown = json.loads(json.dumps(definitions))
    unknown["peptideatlas-dog-peptides"]["versions"]["fixture"]["manifest"]["processing"]["importer"] = "none"
    with pytest.raises(ValueError, match="No packaged importer"):
        DatasetRegistry(unknown)


def test_a_build_only_default_can_sit_beside_published_versions(tmp_path):
    manifest, _ = _source(ATLAS)
    definitions = _build_only(manifest)
    published = replace(manifest, version="published")
    definitions["peptideatlas-dog-peptides"]["versions"]["published"] = dict(
        manifest=published.to_dict(), distribution=dict(access_terms="Test artifact", redistributable=True),
        assets={"dataset.sqlite": dict(url="https://example.invalid/fixture.sqlite", sha256="0" * 64, size=1)})
    registry = DatasetRegistry(definitions, cache_root=tmp_path)
    assert not registry.published("peptideatlas-dog-peptides")
    assert registry.published("peptideatlas-dog-peptides", "published")
    assert registry.inspect("peptideatlas-dog-peptides", "published").status == "missing"
    with pytest.raises(ValueError, match="no published artifact, because its source terms"):
        registry.download("peptideatlas-dog-peptides")
