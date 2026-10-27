"""Lossless numeric vectors with separate presence bits for dataset storage.

This encoding has no vendor dependency. Values are little-endian fixed-width
numbers; bit i (least-significant first) in ``present`` describes element i.
Missing payload values and IEEE-754 bits are retained, including NaN payloads.
A null source vector is distinct from a present vector of missing elements.
"""

import math
import struct
from dataclasses import dataclass


@dataclass(frozen=True)
class NumericVector:
    value_type: str
    count: int
    values: bytes
    present: bytes
    source_null: bool = False

    FORMATS = {"int32": "<i", "int64": "<q", "float64": "<d"}
    ENCODING = "orfdata.numeric-vector.v1"

    def __post_init__(self):
        if self.value_type not in self.FORMATS or type(self.count) is not int or self.count < 0:
            raise ValueError("Invalid numeric vector type/count")
        if type(self.source_null) is not bool or not isinstance(self.values, bytes) or not isinstance(self.present, bytes):
            raise TypeError("Numeric vector payloads must be bytes with an explicit null flag")
        width = struct.calcsize(self.FORMATS[self.value_type])
        if self.source_null:
            if self.values or self.present:
                raise ValueError("A null numeric vector has no element payloads")
        elif len(self.values) != width * self.count or len(self.present) != (self.count + 7) // 8:
            raise ValueError("Numeric vector length differs from its declared type/count")
        elif self.count % 8 and self.present[-1] >> (self.count % 8):
            raise ValueError("Numeric vector has nonzero padding bits")

    def elements(self):
        """Yield raw numeric values and their independent presence flags."""
        if self.source_null:
            yield from ((None, False) for _ in range(self.count))
            return
        for index, (value,) in enumerate(struct.iter_unpack(self.FORMATS[self.value_type], self.values)):
            yield value, bool(self.present[index // 8] & (1 << (index % 8)))

    def cells(self):
        """Yield JSON-safe values with explicit missing/nonfinite source states."""
        for native, present in self.elements():
            nonfinite = isinstance(native, float) and not math.isfinite(native)
            status = "missing" if not present else "nonfinite" if nonfinite else "measured"
            original = ("NaN" if math.isnan(native) else "Infinity" if native > 0 else "-Infinity") if nonfinite else native
            yield dict(value=native if status == "measured" else None, status=status,
                       native_value=original, native_present=present,
                       missing_reason="null_vector" if self.source_null else "null_element" if not present else None)
