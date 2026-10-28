"""Source-native field and enum definitions; confidence domains stay separate."""

import math


class NativeDefinitions:
    def __init__(self, data_types, columns, enum_types, enum_values):
        self.tables = self._unique(data_types, lambda row: row["TableName"])
        self.columns = self._unique(columns, lambda row: (row["DataTypeID"], row["DBColumnName"]))
        self.enums = self._unique(enum_types, lambda row: row["EnumID"])
        self.values = self._unique(enum_values, lambda row: (row["EnumID"], row["Value"]))
        if any(key[0] not in self.enums for key in self.values):
            raise ValueError("Enum value has no native domain definition")
        for row in self.columns.values():
            if row["SpecialValueType"] == "Enum" and row["SpecialValueTypeID"] not in self.enums:
                raise ValueError("Native field references a missing enum domain")

    @staticmethod
    def _unique(rows, key):
        result = {}
        for row in rows:
            identifier = key(row)
            if identifier in result:
                raise ValueError("Duplicate native definition: " + repr(identifier))
            result[identifier] = row
        return result

    def field(self, table, column):
        """Return the deposited column definition, including its native domain ID."""
        return dict(self.columns[self.tables[table]["DataTypeID"], column])

    def scalar(self, table, column, value):
        """Describe a scalar without translating unrelated native confidence scales.

        Unknown enum codes are explicit. A numeric distribution such as protein
        FDR confidence must be decoded with its map metadata; a peptide/PSM
        enum must not be borrowed to label its integer values.
        """
        definition = self.field(table, column)
        if definition["SpecialValueType"] == "DataDistribution" or isinstance(value, bytes):
            raise ValueError("Numeric distributions require their own native map, not a scalar enum")
        if not (value is None or type(value) in (int, float, str)):
            raise TypeError("Expected a native SQLite scalar")
        nonfinite = isinstance(value, float) and not math.isfinite(value)
        native = ("NaN" if math.isnan(value) else "Infinity" if value > 0 else "-Infinity") if nonfinite else value
        result = dict(column_id=definition["ColumnID"], value=None if nonfinite else value,
                      native_value=native, status="missing" if value is None else "nonfinite" if nonfinite else "native",
                      definition=definition)
        if definition["SpecialValueType"] != "Enum":
            return result
        domain_id = definition["SpecialValueTypeID"]
        domain = self.enums[domain_id]
        result.update(enum_id=domain_id, enum_type=domain["EnumType"], enum_definition=dict(domain), labels=[])
        if value is None:
            return result
        if type(value) is not int:
            raise TypeError("Native enum codes must be integers")
        exact = self.values.get((domain_id, value))
        if exact is not None:
            result.update(status="defined", labels=[dict(exact)])
        elif domain["IsFlagsEnum"] and value >= 0:
            known, labels = 0, []
            for (enum, code), row in sorted(self.values.items()):
                if enum == domain_id and code > 0 and code & (code - 1) == 0 and value & code:
                    known |= code
                    labels.append(dict(row))
            unknown = value & ~known
            result.update(status="partly_defined" if unknown else "defined_flags", labels=labels, unknown_bits=unknown)
        else:
            result.update(status="undefined_enum_code")
        return result

    def modification(self, record, link):
        """Keep native modification specificity and position without guessing offsets.

        The link's Position is retained verbatim. Terminal specificity comes from
        its own enum; this does not declare a residue-index coordinate system or
        infer a localized modification from an ambiguous display string.
        """
        if record["ModificationID"] != link["FoundModificationsModificationID"]:
            raise ValueError("Modification link and definition disagree")
        return dict(native_modification=dict(record), native_link=dict(link), source_position=link["Position"],
                    specificity=self.scalar("FoundModifications", "PositionType", record["PositionType"]),
                    position_interpretation="source Position retained without residue-offset conversion")
