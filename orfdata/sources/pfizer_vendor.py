"""Source-defined numeric distributions in the Pfizer Proteome Discoverer result.

Types and ordered labels come from the deposited database metadata. No file
acquisition occurs here. This module does not infer protein attribution or turn
an absent reporter into a measured zero.
"""

import math
import struct


class DistributionDecoder:
    """Decode PD numeric value/presence vectors against their exact metadata.

    PD 2.4 stores each numeric element followed by a presence byte. Int values
    occupy four bytes, Double and Int64 eight; the database declares the type.
    Element type must not be guessed from the vector length. These storage facts
    are independently documented by proteinDiscover's transformRaw.R; native
    fixture bytes and distribution metadata supply this reader's validation.
    """

    def __init__(self, data_types, columns, maps, boxes, enum_types, enum_values):
        self.tables = self._unique(data_types, "TableName")
        self.maps = self._unique(maps, "ID")
        enums = [r for r in enum_types if r["EnumType"].split(",")[0].endswith(".CustomDataType")]
        if len(enums) != 1:
            raise ValueError("Expected one native CustomDataType enum")
        self.numeric_types = {r["Value"]: r["DisplayName"] for r in enum_values
                              if r["EnumID"] == enums[0]["EnumID"]}
        self.columns = {}
        for column in columns:
            key = column["DataTypeID"], column["DBColumnName"]
            if key in self.columns:
                raise ValueError("Duplicate native column identity")
            self.columns[key] = column
        self.boxes = {}
        for box in boxes:
            self.boxes.setdefault(box["DataDistributionMapID"], []).append(box)
        for values in self.boxes.values():
            values.sort(key=lambda b: b["Position"])
            if [b["Position"] for b in values] != list(range(1, len(values) + 1)):
                raise ValueError("Distribution positions must be complete and unique")
            if len({b["BoxID"] for b in values}) != len(values):
                raise ValueError("Duplicate native distribution box identity")

    @staticmethod
    def _unique(rows, field):
        result = {}
        for row in rows:
            if row[field] in result:
                raise ValueError("Duplicate native metadata identity: " + field)
            result[row[field]] = row
        return result

    def decode(self, table, column, blob):
        """Yield labelled cells with finite values or explicit missing/nonfinite status."""
        definition = self.columns[(self.tables[table]["DataTypeID"], column)]
        if definition["SpecialValueType"] != "DataDistribution":
            raise ValueError("Column is not a native data distribution")
        distribution = self.maps[definition["SpecialValueTypeID"]]
        boxes = self.boxes[distribution["ID"]]
        kind = self.numeric_types[distribution["CustomDataType"]]
        formats = {"Int": "<iB", "Int64": "<qB", "Double": "<dB"}
        if kind not in formats:
            raise ValueError("Unsupported declared distribution type: " + kind)
        record = struct.Struct(formats[kind])
        if blob is not None and (not isinstance(blob, bytes) or len(blob) != record.size * len(boxes)):
            raise ValueError("Distribution length differs from its declared type and ordered boxes")
        elements = [(None, 0)] * len(boxes) if blob is None else record.iter_unpack(blob)
        for box, (native, present) in zip(boxes, elements):
            if present not in (0, 1):
                raise ValueError("Invalid native distribution presence flag")
            nonfinite = isinstance(native, float) and not math.isfinite(native)
            status = "missing" if not present else "nonfinite" if nonfinite else "measured"
            original = ("NaN" if math.isnan(native) else "Infinity" if native > 0 else "-Infinity") if nonfinite else native
            yield dict(column_id=definition["ColumnID"], distribution_id=distribution["ID"],
                       box_id=box["BoxID"], position=box["Position"], label=box["Name"],
                       value=native if status == "measured" else None, status=status,
                       native_value=original, native_present=bool(present),
                       missing_reason="null_vector" if blob is None else "null_element" if not present else None)
