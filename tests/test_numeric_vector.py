"""Normalized vectors retain source bits, missing payloads and null-vector identity."""

import base64
import json
import math
import struct
from pathlib import Path

import pytest

from orfdata.numeric_vector import NumericVector
from orfdata.sources.pfizer_vendor import DistributionDecoder


def test_every_native_fixture_numeric_vector_reconstructs_exact_original_bytes():
    tables = json.loads((Path(__file__).parent / "data/pfizer/vendor/native.json").read_text())["tables"]
    decoder = DistributionDecoder(*(tables[t] for t in (
        "DataTypes", "DataTypesColumns", "DataDistributionMaps", "DataDistributionBoxes",
        "EnumDataTypes", "EnumDataTypeValues")))
    observed_types = set()
    for table in ("TargetPsms", "TargetPeptideGroups", "TargetProteins"):
        columns = [r["DBColumnName"] for r in tables["DataTypesColumns"]
                   if r["DataTypeID"] == decoder.tables[table]["DataTypeID"]]
        for row in tables[table]:
            for column in columns:
                if column not in row:
                    continue
                encoded = row[column]
                original = base64.b64decode(encoded["base64"]) if encoded is not None else None
                vector = decoder.vector(table, column, original)
                assert decoder.native_bytes(vector) == original
                observed_types.add(vector.value_type)
                decoded = list(decoder.decode(table, column, original))
                assert [cell["value"] for cell in vector.cells()] == [cell["value"] for cell in decoded]
    assert observed_types == {"int32", "float64"}


def test_nan_payloads_infinities_missing_values_and_signed_int64_are_lossless():
    # Non-default IEEE NaN payload remains byte-for-byte unchanged even when
    # its common JSON value is represented as nonfinite rather than a number.
    payloads = struct.pack("<Qdd", 0x7FF80000000000A5, float("-inf"), 7.5)
    vector = NumericVector("float64", 3, payloads, b"\x03")
    values = list(vector.elements())
    assert math.isnan(values[0][0]) and values[0][1]
    assert values[1] == (float("-inf"), True) and values[2] == (7.5, False)
    assert DistributionDecoder.native_bytes(vector) == payloads[:8] + b"\x01" + payloads[8:16] + b"\x01" + payloads[16:] + b"\x00"
    cells = list(vector.cells())
    assert [c["status"] for c in cells] == ["nonfinite", "nonfinite", "missing"]
    assert [c["native_value"] for c in cells] == ["NaN", "-Infinity", 7.5]
    json.dumps(cells, allow_nan=False)
    original = -9222281927842617047
    integer = NumericVector("int64", 1, struct.pack("<q", original), b"\x01")
    assert list(integer.elements()) == [(original, True)]
    assert list(integer.cells())[0]["value"] == original


def test_null_vector_is_distinct_from_all_missing_and_measured_zero():
    null = NumericVector("float64", 11, b"", b"", True)
    absent = NumericVector("float64", 11, bytes(88), bytes(2))
    zero = NumericVector("float64", 11, bytes(88), b"\xff\x07")
    assert DistributionDecoder.native_bytes(null) is None
    assert DistributionDecoder.native_bytes(absent) == bytes(99)
    assert all(c["missing_reason"] == "null_vector" and c["native_value"] is None for c in null.cells())
    assert all(c["missing_reason"] == "null_element" and c["native_value"] == 0.0 for c in absent.cells())
    assert all(c["status"] == "measured" and c["value"] == 0.0 for c in zero.cells())


@pytest.mark.parametrize("arguments", [
    ("float32", 1, bytes(4), b"\x01"),
    ("int32", True, bytes(4), b"\x01"),
    ("int32", 1, bytes(3), b"\x01"),
    ("int32", 1, bytes(4), b""),
    ("int32", 1, bytes(4), b"\x81"),
    ("int32", 1, bytes(4), b"\x01", True),
])
def test_malformed_normalized_vectors_are_rejected(arguments):
    with pytest.raises((ValueError, TypeError)):
        NumericVector(*arguments)
