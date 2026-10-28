"""Native labels require their own domains and preserve undefined source codes."""

import copy
import json
from pathlib import Path

import pytest

from orfdata.sources.pfizer_definitions import NativeDefinitions

SOURCE = json.loads((Path(__file__).parent / "data/pfizer/vendor/definitions.json").read_text())


def definitions(tables=None):
    tables = tables or SOURCE["tables"]
    return NativeDefinitions(*(tables[t] for t in ("DataTypes", "DataTypesColumns", "EnumDataTypes", "EnumDataTypeValues")))


def test_confidence_ambiguity_master_and_quantification_codes_are_separate_domains():
    native = definitions()
    examples = [
        ("TargetPsms", "MatchConfidence", 3, 12, "High"),
        ("TargetPsms", "PsmAmbiguity", 3, 11, "Ambiguous"),
        ("TargetProteins", "IsMasterProtein", 3, 14, "Master Protein Rejected"),
        ("TargetPsms", "QuanInfo", 3, 19, "Redundant"),
        ("TargetPeptideGroups", "Confidence", 0, 12, "Unknown"),
    ]
    for table, column, value, domain, label in examples:
        result = native.scalar(table, column, value)
        assert result["value"] == value and result["enum_id"] == domain
        assert result["labels"][0]["DisplayName"] == label and result["status"] == "defined"
    assert native.scalar("TargetPsms", "ExcludedBy", -1)["status"] == "native"
    assert "node number" in native.field("TargetPsms", "ExcludedBy")["Property_Description"]
    with pytest.raises(ValueError, match="own native map"):
        native.scalar("TargetProteins", "ProteinFDRConfidence", 3)


def test_missing_unknown_and_combined_flags_do_not_gain_invented_labels():
    native = definitions()
    missing = native.scalar("TargetPsms", "QuanInfo", None)
    assert missing["status"] == "missing" and missing["labels"] == [] and missing["value"] is None
    unknown = native.scalar("TargetPsms", "MatchConfidence", 999)
    assert unknown["status"] == "undefined_enum_code" and unknown["value"] == 999 and unknown["labels"] == []
    # ActivationTypes is a source-declared flag enum; combine two known single bits.
    domain = native.field("TargetPsms", "ActivationTypes")["SpecialValueTypeID"]
    codes = [code for enum, code in native.values if enum == domain and code > 0 and code & (code - 1) == 0]
    value = codes[0] | codes[1]
    flags = native.scalar("TargetPsms", "ActivationTypes", value)
    assert flags["status"] == "defined_flags" and flags["unknown_bits"] == 0
    assert {row["Value"] for row in flags["labels"]} == set(codes[:2])
    partial = native.scalar("TargetPsms", "ActivationTypes", value | (1 << 60))
    assert partial["status"] == "partly_defined" and partial["unknown_bits"] == 1 << 60
    with pytest.raises(TypeError, match="integers"):
        native.scalar("TargetPsms", "MatchConfidence", 3.0)


def test_native_modification_specificity_does_not_reinterpret_its_position():
    native = definitions()
    result = native.modification(SOURCE["modification"], SOURCE["modification_link"])
    assert result["native_modification"]["Name"] == "Acetyl"
    assert result["native_modification"]["UnimodAccession"] == 1
    assert result["native_modification"]["DeltaMonoisotopicMass"] == 42.010565
    assert result["source_position"] == 1
    assert result["native_link"]["TargetPsmsWorkflowID"] == -160
    assert result["native_link"]["TargetPsmsPeptideID"] == 360
    assert result["specificity"]["labels"][0]["DisplayName"] == "Protein N-Terminus"
    wrong = dict(SOURCE["modification_link"], FoundModificationsModificationID=3)
    with pytest.raises(ValueError, match="disagree"):
        native.modification(SOURCE["modification"], wrong)


def test_definition_integrity_and_nonfinite_scalars_remain_explicit():
    tables = copy.deepcopy(SOURCE["tables"])
    tables["EnumDataTypes"] = [r for r in tables["EnumDataTypes"] if r["EnumID"] != 12]
    with pytest.raises(ValueError, match="domain"):
        definitions(tables)
    tables = copy.deepcopy(SOURCE["tables"])
    tables["EnumDataTypeValues"].append(tables["EnumDataTypeValues"][0])
    with pytest.raises(ValueError, match="Duplicate"):
        definitions(tables)
    result = definitions().scalar("TargetPsms", "PercolatorPEP", float("inf"))
    assert result["status"] == "nonfinite" and result["value"] is None and result["native_value"] == "Infinity"
    json.dumps(result, allow_nan=False)
