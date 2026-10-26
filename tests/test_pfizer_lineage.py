"""Source-backed pool roles and native acquisition identities stay distinct."""

import base64
import copy
import hashlib
import json
import struct
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import QuantitativeDataset, SourceAsset
from orfdata.sources.pfizer import protein_manifest
from orfdata.sources.pfizer_lineage import (
    lineage_records,
    native_assay_id,
    native_sample_id,
    reconcile_dog_lineage,
)
from orfdata.sources.pfizer_vendor import DistributionDecoder, NativeQuantification
from tools.export_dog_lineage import (
    DOCREL,
    PANEL_SHA256,
    REL,
    SHEET,
    SML,
    dog_panel,
    study_metadata,
)

DIRECTORY = Path(__file__).parent / "data/pfizer/vendor"


def sources():
    return (study_metadata(DIRECTORY / "study-slice.xml"),
            json.loads((DIRECTORY / "dog-panel.json").read_text()),
            json.loads((DIRECTORY / "summary-samples.json").read_text())["samples"],
            json.loads((DIRECTORY / "native.json").read_text())["tables"]["WorkflowInputFiles"])


def test_native_channels_runs_and_summary_pools_reconcile_without_individual_dog_claims():
    result = reconcile_dog_lineage(*sources())
    assert result["counts"] == dict(channels=44, runs=96, tissue_pool=38, bridge=4, empty=2)
    channels = {(c["plex"], c["channel"]): c for c in result["channels"]}
    adrenal = channels["Plex1", "126"]
    assert adrenal["file_set_number"] == 5
    assert adrenal["file_set_id"] == "a1baa969-e535-4cc8-80bc-ba874d9e4cd4"
    assert adrenal["tissue"] == "Adrenal" and adrenal["tissue_pool_id"] == "RBI-0000137273"
    assert adrenal["summary_sample"]["Run"] == "Pool1_1"
    assert adrenal["publication_cell"] == "N2"
    assert adrenal["pool_members"] == [] and adrenal["pool_membership"] == "unresolved"
    for plex in ("Plex1", "Plex4"):
        empty = channels[plex, "131N"]
        assert empty["role"] == "empty" and empty["vendor_sample_type"] == "Sample"
        assert empty["tissue"] is empty["summary_sample_id"] is empty["tissue_pool_id"] is None
    for plex in ("Plex1", "Plex2", "Plex3", "Plex4"):
        bridge = channels[plex, "131C"]
        assert bridge["role"] == "bridge" and bridge["tissue"] is None
        assert bridge["summary_sample"] is None and bridge["pool_members"] == []
    runs = {r["study_file_id"]: r for r in result["runs"]}
    assert (runs["F5.1"]["workflow_id"], runs["F5.1"]["file_id"]) == (-160, 2793)
    assert runs["F5.1"]["file_name"] == "IQP1-19058_1269D-P1-F1-RTS.raw"
    assert runs["F5.1"]["file_size"] == 945301680
    assert runs["F5.1"]["fraction"] == 1 and runs["F5.1"]["technical_replicate"] is None
    reports = result["instrument_reconciliation"]
    assert reports["native_run_counts"] == {"Orbitrap Fusion Lumos": 96}
    assert reports["publication"]["instrument"] == "Orbitrap Fusion Eclipse"
    assert reports["deposit"]["instrument"] == "Orbitrap Fusion Lumos"
    json.dumps(result, allow_nan=False)


def test_summary_mixtures_match_entire_profiles_independently_of_labels():
    study, panel, samples, files = sources()
    for sample in samples:
        sample["Mixture"] = "arbitrary-label-" + sample["Mixture"][-1]
    result = reconcile_dog_lineage(study, panel, samples, files)
    assert result["counts"]["tissue_pool"] == 38
    first = result["channels"][0]
    assert first["summary_sample"]["Mixture"] == "arbitrary-label-1"
    assert "profile equality" in first["summary_to_publication_basis"]
    samples[0]["Condition"] = "Bladder"
    with pytest.raises(ValueError, match="complete channel/tissue profile"):
        reconcile_dog_lineage(study, panel, samples, files)


def test_common_api_keeps_pool_identity_across_24_fractions_and_separates_controls():
    lineage = reconcile_dog_lineage(*sources())
    reference, = protein_manifest("dog").references
    samples, assays = lineage_records(lineage, reference)
    assert len(samples) == 44 and len(assays) == 44 + 96 * 11
    assert len({s.sample_id for s in samples}) == 44
    assert len({a.assay_id for a in assays}) == len(assays)
    adrenal = next(c for c in lineage["channels"] if c["tissue"] == "Adrenal")
    sample = next(s for s in samples if s.sample_id == native_sample_id(adrenal))
    assert sample.biological_unit_id == "MSV000099744:pool:RBI-0000137273"
    assert sample.pool_members == () and sample.biological_replicate is None
    related = [a for a in assays if a.sample_id == sample.sample_id]
    assert len(related) == 25
    aggregate = next(a for a in related if a.assay_id == native_assay_id(adrenal))
    assert aggregate.run_id is None and aggregate.fraction is None
    runs = [a for a in related if a.run_id is not None]
    assert {a.fraction for a in runs} == {str(i) for i in range(1, 25)}
    assert all(a.technical_replicate is None for a in related)
    assert all(s.tissue is None and s.biological_unit_id is None for s in samples
               if s.attributes["role"] in ("empty", "bridge"))
    assert all(a.kind == "digest_ms" and a.bridge == "131C" for a in assays)


def quantification():
    lineage = reconcile_dog_lineage(*sources())
    tables = json.loads((DIRECTORY / "native.json").read_text())["tables"]
    decoder = DistributionDecoder(*(tables[t] for t in (
        "DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
        "EnumDataTypes", "EnumDataTypeValues")))
    reference, = protein_manifest("dog").references
    converter = NativeQuantification(decoder, lineage, reference, "native-quantification-fixture", "test",
                                     source_asset="native.json")
    for table in ("TargetPsms", "TargetPeptideGroups"):
        tables[table] = [{key: base64.b64decode(value["base64"]) if isinstance(value, dict)
                         and value.get("type") == "blob" else value for key, value in row.items()}
                        for row in tables[table]]
    return lineage, tables, converter


def test_real_native_observations_roundtrip_through_common_dataset_api(tmp_path):
    lineage, tables, converter = quantification()
    psms = [m for row in tables["TargetPsms"] for m in converter.peptide_abundances("TargetPsms", row)]
    group = list(converter.peptide_abundances("TargetPeptideGroups", tables["TargetPeptideGroups"][0]))
    assert len(psms) == 33 and len(group) == 44
    first = psms[0]
    assert first.value == 7496.059512428747 and first.target.identifier == "GHEEEHHHR"
    assert first.unit == "PD2.4 TargetPsms.Abundances"
    assert first.provenance["native_record"] == {"WorkflowID": -160, "PeptideID": 1}
    assert first.provenance["native_quality"]["PercolatorqValue"] == 0.0003468
    assert first.provenance["native_quality"]["PercolatorPEP"] == 0.002091
    assert first.provenance["native_quality"]["MatchConfidence"] == 3
    assert first.provenance["native_modifications"] == "N-Term(TMT6plex)"
    assert all(m.value is None and m.status == "missing" for m in psms[1:11])
    assert psms[9].provenance["channel_role"] == "empty"
    assert psms[10].provenance["channel_role"] == "bridge"
    assert all(m.status == "missing" and m.provenance["native_distribution"]["missing_reason"] == "null_vector"
               for m in group)
    assert {m.unit for m in group} == {"PD2.4 TargetPeptideGroups.Abundances"}
    assets = tuple(SourceAsset(p.name, p.resolve().as_uri(), hashlib.sha256(p.read_bytes()).hexdigest(),
                               p.stat().st_size, "Factual test slice; see vendor/NOTICE.md")
                   for p in (DIRECTORY / name for name in
                             ("native.json", "study-slice.xml", "dog-panel.json", "summary-samples.json")))
    manifest = replace(protein_manifest("dog"), dataset_id=converter.dataset_id, version=converter.dataset_version,
                       title="Explicit native quantification test slice", sources=assets,
                       coverage=dict(scope="Three native PSMs and one unquantified peptide group"),
                       processing=dict(importer="native quantification test"))
    samples, assays = lineage_records(lineage, converter.reference)
    dataset = QuantitativeDataset.create(tmp_path / "native.sqlite", manifest, samples=samples, assays=assays,
                                         measurements=psms + group)
    saved = {m.measurement_id: m for m in dataset.measurements()}
    assert saved == {m.measurement_id: m for m in psms + group}
    observed = list(dataset.measurements(target=first.target, unit=first.unit, status="measured"))
    assert observed == [first]
    source_sample = next(s for s in dataset.samples() if s.sample_id == first.sample_id)
    assert source_sample.tissue == "Adrenal" and source_sample.pool_members == ()


def test_quantification_refuses_wrong_run_or_channel_instead_of_guessing():
    _, tables, converter = quantification()
    row = copy.deepcopy(tables["TargetPsms"][0])
    row["StudyFileId"] = "F6.1"
    with pytest.raises(ValueError, match="acquisition identity"):
        list(converter.peptide_abundances("TargetPsms", row))
    row["SpectrumFileId"] = 99999999
    with pytest.raises(ValueError, match="unresolved native acquisition"):
        list(converter.peptide_abundances("TargetPsms", row))
    with pytest.raises(ValueError, match="Unsupported native peptide abundance"):
        list(converter.peptide_abundances("TargetPsms", row, "Intensity"))
    del converter.channels[5, "126"]
    with pytest.raises(ValueError, match="lacks reconciled lineage"):
        list(converter.peptide_abundances("TargetPsms", tables["TargetPsms"][0]))


def test_common_nonfinite_status_and_native_normalizations_stay_explicit():
    _, tables, converter = quantification()
    row = copy.deepcopy(tables["TargetPsms"][0])
    # Synthetic numeric boundary cases in the real native eleven-channel layout.
    row["Abundances"] = struct.pack("<dB", float("inf"), 1) + struct.pack("<dB", 0.0, 1) * 10
    values = list(converter.peptide_abundances("TargetPsms", row))
    assert values[0].value is None and values[0].status == "missing"
    assert values[0].provenance["native_distribution"]["status"] == "nonfinite"
    assert values[0].provenance["native_distribution"]["native_value"] == "Infinity"
    assert values[1].value == 0.0 and values[1].status == "measured"
    group = tables["TargetPeptideGroups"][0]
    native = {column: list(converter.peptide_abundances("TargetPeptideGroups", group, column))
              for column in ("Abundances", "AbundancesNormalized", "AbundancesScaled")}
    assert {len(rows) for rows in native.values()} == {44}
    assert len({rows[0].unit for rows in native.values()}) == 3
    assert {rows[0].provenance["native_distribution"]["distribution_id"] for rows in native.values()} == {39, 45, 47}


def test_panel_reader_preserves_cells_without_recalculation(tmp_path):
    # Minimal synthetic OOXML wrapper of the factual dog cells. The publisher's
    # complete workbook/layout is not redistributed in the MIT source tree.
    source = json.loads((DIRECTORY / "dog-panel.json").read_text())
    assert source["source_sha256"] == PANEL_SHA256
    workbook = ET.Element(SML + "workbook")
    ET.SubElement(ET.SubElement(workbook, SML + "sheets"), SML + "sheet",
                  {"name": SHEET, DOCREL + "id": "rId2"})
    relations = ET.Element(REL + "Relationships")
    ET.SubElement(relations, REL + "Relationship", {"Id": "rId2", "Target": "worksheets/sheet2.xml"})
    worksheet, strings = ET.Element(SML + "worksheet"), ET.Element(SML + "sst")
    data = ET.SubElement(worksheet, SML + "sheetData")
    for address, value in source["cells"].items():
        if address == "M2":
            cell = ET.SubElement(data, SML + "c", {"r": address})
            ET.SubElement(cell, SML + "v").text = value
        else:
            cell = ET.SubElement(data, SML + "c", {"r": address, "t": "s"})
            ET.SubElement(cell, SML + "v").text = str(len(strings))
            ET.SubElement(ET.SubElement(strings, SML + "si"), SML + "t").text = value
    path = tmp_path / "panel.xlsx"
    with zipfile.ZipFile(path, "w") as archive:
        for name, root in (("xl/workbook.xml", workbook), ("xl/_rels/workbook.xml.rels", relations),
                           ("xl/worksheets/sheet2.xml", worksheet), ("xl/sharedStrings.xml", strings)):
            archive.writestr(name, ET.tostring(root))
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    panel = dog_panel(path)
    assert before == panel["source_sha256"]
    assert panel["cells"] == source["cells"]
    assert panel["cells"]["N11"] == panel["cells"]["Q11"] == "Empty"
    assert panel["cells"]["O11"] == "Uterus_Cervix"
    assert panel["cells"]["P11"] == "Endothelium"
    assert panel["cells"]["M2"] == "126"
    assert panel["searched_fasta_note"] == (
        "FASTA: Canis lupus familiaris from Uniprot, downloaded July 2023, 20,649 sequences")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before


@pytest.mark.parametrize("change", ["duplicate", "size", "basename", "missing", "native_key"])
def test_run_join_rejects_inconsistent_or_incomplete_original_ids(change):
    study, panel, samples, files = sources()
    raw = [r for r in files if r["FileType"] == "XcaliburRawfile"]
    if change == "duplicate":
        files.append(copy.deepcopy(raw[0]))
    elif change == "size":
        raw[0]["FileSize"] += 1
    elif change == "basename":
        raw[0]["FileName"] = "different.raw"
    elif change == "missing":
        files.remove(raw[0])
    else:
        raw[0]["WorkflowID"], raw[0]["FileID"] = raw[1]["WorkflowID"], raw[1]["FileID"]
    with pytest.raises(ValueError, match="Duplicate|identity/size|run coverage"):
        reconcile_dog_lineage(study, panel, samples, files)


def test_missing_native_channel_and_changed_publication_role_are_not_guessed():
    study, panel, samples, files = sources()
    study["samples"] = [s for s in study["samples"] if s["Name"] != "IQP_1269_Dog_P1 - [131N]"]
    with pytest.raises(ValueError, match="eleven distinct channels"):
        reconcile_dog_lineage(study, panel, samples, files)
    study, panel, samples, files = sources()
    next(r for r in panel["channels"] if r["cell"] == "N11")["source_label"] = "Adrenal"
    with pytest.raises(ValueError, match="empty channels"):
        reconcile_dog_lineage(study, panel, samples, files)
