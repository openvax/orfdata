"""Native PD metadata governs values, labels and missingness."""

import base64
import copy
import hashlib
import json
import struct
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata.search_reference import SearchReference
from orfdata.sources.pfizer import protein_manifest, search_references
from orfdata.sources.pfizer_vendor import DistributionDecoder, peptide_candidates

NATIVE = json.loads((Path(__file__).parent / "data/pfizer/vendor/native.json").read_text())
TABLES = NATIVE["tables"]


def decoder(tables=TABLES):
    return DistributionDecoder(*(tables[t] for t in (
        "DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
        "EnumDataTypes", "EnumDataTypeValues")))


def blob(table, field, row=0):
    value = TABLES[table][row][field]
    return base64.b64decode(value["base64"]) if value is not None else None


def test_native_psm_reporters_keep_missing_values_and_channel_order():
    psm = TABLES["TargetPsms"][0]
    assert (psm["WorkflowID"], psm["PeptideID"], psm["Sequence"]) == (-160, 1, "GHEEEHHHR")
    values = list(decoder().decode("TargetPsms", "Abundances", blob("TargetPsms", "Abundances")))
    assert len(values) == 11
    assert values[0]["value"] == 7496.059512428747
    assert values[0]["status"] == "measured" and values[0]["distribution_id"] == 43
    assert all(v["value"] is None and v["status"] == "missing" for v in values[1:])
    assert all(v["native_value"] == 0 and v["native_present"] is False for v in values[1:])
    assert all(v["missing_reason"] == "null_element" for v in values[1:])
    assert "126" in values[0]["label"] and "131C" in values[-1]["label"]


def test_native_protein_sample_vectors_and_counts_are_distinct():
    dec = decoder()
    abundance = list(dec.decode("TargetProteins", "Abundances", blob("TargetProteins", "Abundances")))
    normalized = list(dec.decode("TargetProteins", "AbundancesNormalized", blob("TargetProteins", "AbundancesNormalized")))
    counts = list(dec.decode("TargetProteins", "AbundancesCounts", blob("TargetProteins", "AbundancesCounts")))
    assert len(abundance) == len(normalized) == len(counts) == 44
    assert abundance[0]["value"] == 32552.227340917958
    assert normalized[0]["value"] == 77329.4733908357
    assert counts[0]["value"] == 2 and isinstance(counts[0]["value"], int)
    assert len({abundance[0]["column_id"], normalized[0]["column_id"], counts[0]["column_id"]}) == 3
    assert abundance[0]["label"] == "F5: 126, Sample"
    assert abundance[-1]["label"] == "F8: 131C, Sample"
    # An 11-channel PSM array cannot be interpreted as the 44-sample protein array.
    with pytest.raises(ValueError, match="Distribution length"):
        list(dec.decode("TargetProteins", "Abundances", blob("TargetPsms", "Abundances")))


def test_native_unquantified_peptide_is_a_missing_vector():
    assert TABLES["TargetPeptideGroups"][0]["PeptideGroupID"] == 1
    assert blob("TargetPeptideGroups", "Abundances") is None
    values = list(decoder().decode("TargetPeptideGroups", "Abundances", None))
    assert len(values) == 44 and all(v["status"] == "missing" for v in values)
    assert all(v["missing_reason"] == "null_vector" for v in values)
    assert all(v["native_value"] is None for v in values)


def test_present_nonfinite_values_and_missing_payloads_remain_distinct():
    data = (struct.pack("<dB", float("nan"), 1) + struct.pack("<dB", float("inf"), 1)
            + struct.pack("<dB", 7.0, 0) + struct.pack("<dB", 0.0, 1) * 8)
    values = list(decoder().decode("TargetPsms", "Abundances", data))
    assert [v["status"] for v in values[:4]] == ["nonfinite", "nonfinite", "missing", "measured"]
    assert [v["native_value"] for v in values[:4]] == ["NaN", "Infinity", 7.0, 0.0]
    assert [v["value"] for v in values[:4]] == [None, None, None, 0.0]
    json.dumps(values, allow_nan=False)
    bad = data[:8] + b"\x02" + data[9:]
    with pytest.raises(ValueError, match="presence flag"):
        list(decoder().decode("TargetPsms", "Abundances", bad))
    with pytest.raises(ValueError, match="Distribution length"):
        list(decoder().decode("TargetPsms", "Abundances", data[:-1]))


def test_inconsistent_metadata_is_rejected():
    tables = copy.deepcopy(TABLES)
    tables["DataDistributionBoxes"].append(tables["DataDistributionBoxes"][0])
    with pytest.raises(ValueError, match="positions"):
        decoder(tables)
    tables = copy.deepcopy(TABLES)
    tables["DataTypesColumns"].append(tables["DataTypesColumns"][0])
    with pytest.raises(ValueError, match="column identity"):
        decoder(tables)


def fixture_search(path, count):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    original, = search_references(protein_manifest("dog"))
    source = replace(original.source, url=path.as_uri(), size=path.stat().st_size, sha256=digest)
    reference = replace(original.reference, sequence_sha256=digest,
                        annotation_name="explicit test FASTA slice")
    spec = replace(original, dataset_id="pfizer-dog-candidate-fixture", source=source,
                   reference=reference, expected_sequences=count)
    return SearchReference(path, spec)


def test_native_shared_peptide_has_more_search_candidates_than_reported_members():
    directory = Path(__file__).parent / "data/pfizer/vendor"
    audit = json.loads((directory / "search-slice.json").read_text())
    fasta = directory / "search-slice.fasta"
    assert hashlib.sha256(fasta.read_bytes()).hexdigest() == audit["fixture_sha256"]
    reference = fixture_search(fasta, audit["fixture_sequences"])
    peptides = {row["peptide"] for row in audit["candidate_placements"]}
    assert list(peptide_candidates(peptides, reference)) == audit["candidate_placements"]
    shared = next(row for row in TABLES["SharedPeptideExamples"] if row["PeptideGroupID"] == 8)
    reported = {accession.strip() for accession in shared["ParentProteinAccessions"].split(";")}
    candidates = {row["accession"] for row in peptide_candidates([shared["Sequence"]], reference)}
    assert reported == {"A0A5F4BP72", "E2RMY1"}
    assert candidates == reported | {"E2RPC3", "J9NWC2"}
    assert shared["ParentProteinCount"] == 2


def test_candidate_matching_keeps_repeated_sites_il_ambiguity_and_no_matches(tmp_path):
    fasta = tmp_path / "synthetic.fasta"
    fasta.write_text(">sp|ONE|ONE_DOG synthetic OX=9615\nMIIILLLAAA\n"
                     ">tr|TWO|TWO_DOG synthetic OX=9615\nMAAAA\n")
    reference = fixture_search(fasta, 2)
    matches = list(peptide_candidates(["II", "AA", "NOPE", "II"], reference))
    ii = [r for r in matches if r["peptide"] == "II"]
    assert [(r["start"], r["end"]) for r in ii] == [(1, 3), (2, 4), (3, 5), (4, 6), (5, 7)]
    assert [r["basis"] for r in ii] == ["exact_residues"] * 2 + ["isoleucine_leucine_equivalent"] * 3
    assert len(list(peptide_candidates(["II"], reference, equate_il=False))) == 2
    aa = [r for r in matches if r["peptide"] == "AA"]
    assert len(aa) == 5 and {r["accession"] for r in aa} == {"ONE", "TWO"}
    assert matches[-1] == dict(peptide="NOPE", accession=None, start=None, end=None,
                              status="no_search_match", basis="exact_or_il_equivalent")
    assert not list(peptide_candidates([], reference))
    with pytest.raises(ValueError, match="residue budget"):
        list(peptide_candidates(["AAAA"], reference, max_residues=3))
    with pytest.raises(ValueError, match="unmodified uppercase"):
        list(peptide_candidates(["M[Oxidation]AA"], reference))
