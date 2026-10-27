"""Source-defined numeric distributions in the Pfizer Proteome Discoverer result.

Types and ordered labels come from the deposited database metadata. No file
acquisition occurs here. This module does not infer protein attribution or turn
an absent reporter into a measured zero.
"""

import math
import re
import struct

from ..measurements import Measurement
from ..numeric_vector import NumericVector
from ..reference import ScopedTarget
from .pfizer_lineage import native_assay_id, native_sample_id


class NativeQuantification:
    """Attach decoded native abundances to checked pools/runs in the common API.

    Peptide groups and PSMs retain independent source record IDs and native
    quality fields. No confidence threshold, reassignment, re-normalization or
    inference of protein/individual-dog abundance occurs here.
    """

    COLUMNS = {"TargetPeptideGroups": ("Abundances", "AbundancesNormalized", "AbundancesScaled"),
               "TargetPsms": ("Abundances",)}

    def __init__(self, decoder, lineage, reference, dataset_id, dataset_version, *, source_asset):
        if reference.taxon_id != 9615:
            raise ValueError("Canine native quantification requires a canine reference")
        if not isinstance(source_asset, str) or not source_asset.strip():
            raise ValueError("Native quantification requires an explicit source asset ID")
        self.decoder, self.reference = decoder, reference
        self.dataset_id, self.dataset_version = dataset_id, dataset_version
        self.source_asset = source_asset
        self.channels = {(c["file_set_number"], c["channel"]): c for c in lineage["channels"]}
        self.runs = {(r["workflow_id"], r["file_id"]): r for r in lineage["runs"]}
        if len(self.channels) != len(lineage["channels"]) or len(self.runs) != len(lineage["runs"]):
            raise ValueError("Duplicate native quantification lineage identity")

    def peptide_axes(self, table, row, column="Abundances"):
        """Resolve ordered distribution boxes to checked sample/assay identities.

        This validates run identity without expanding measurements or interpreting
        numeric values. It also supplies the shared axes of compact storage.
        """
        if table not in self.COLUMNS or column not in self.COLUMNS[table]:
            raise ValueError("Unsupported native peptide abundance table/column")
        if not isinstance(row["Sequence"], str) or not re.fullmatch("[A-Z]+", row["Sequence"]):
            raise ValueError("Native peptide sequence must be unmodified uppercase residues")
        run = None
        if table == "TargetPsms":
            key = row["WorkflowID"], row["SpectrumFileId"]
            if key not in self.runs:
                raise ValueError("PSM references an unresolved native acquisition")
            run = self.runs[key]
            if (row["StudyFileId"] != run["study_file_id"]
                    or row["SpectrumFileName"] != run["file_name"]):
                raise ValueError("PSM acquisition identity disagrees with native run metadata")
        _, _, boxes = self.decoder.distribution(table, column)
        axes = []
        for box in boxes:
            if run is not None:
                channel_key = run["file_set_number"], box["Name"]
            else:
                label = re.fullmatch(r"F(\d+): (\d+(?:N|C)?), Sample", box["Name"])
                if label is None:
                    raise ValueError("Unrecognized native peptide-group sample label")
                channel_key = int(label[1]), label[2]
            if channel_key not in self.channels:
                raise ValueError("Native quantification channel lacks reconciled lineage")
            channel = self.channels[channel_key]
            axes.append(dict(box_id=box["BoxID"], sample_id=native_sample_id(channel),
                             assay_id=native_assay_id(channel, run), channel_role=channel["role"]))
        return tuple(axes)

    def peptide_abundances(self, table, row, column="Abundances"):
        """Yield every native cell, including null vectors, missing and nonfinite cells.

        All units name the original table and column. In particular PSM
        Abundances is not ReporterQuanResults.Height or the released iBAQ.
        A nonfinite native value has common status missing, with its original
        spelling and nonfinite state preserved explicitly in provenance.
        """
        axes = self.peptide_axes(table, row, column)
        identifiers = (dict(WorkflowID=row["WorkflowID"], PeptideID=row["PeptideID"]) if table == "TargetPsms"
                       else dict(PeptideGroupID=row["PeptideGroupID"]))
        source_id = ":".join(str(value) for value in identifiers.values())
        quality_fields = ("Confidence", "MatchConfidence", "Qvalityqvalue", "QvalityPEP", "SVM_Score",
                          "PercolatorqValue", "PercolatorPEP", "PercolatorSVMScore", "ExcludedBy", "PsmAmbiguity")
        quality = {key: row[key] for key in quality_fields if key in row}
        # SQLite scalar infinities are legal; common metadata must remain JSON-safe.
        for key, value in quality.items():
            if isinstance(value, float) and not math.isfinite(value):
                quality[key] = dict(status="nonfinite", native_value=
                                    "NaN" if math.isnan(value) else "Infinity" if value > 0 else "-Infinity")
        target = ScopedTarget(self.reference, "peptide", row["Sequence"], "amino-acid-sequence")
        definition = self.decoder.columns[self.decoder.tables[table]["DataTypeID"], column]
        for axis, cell in zip(axes, self.decoder.decode(table, column, row[column]), strict=True):
            yield Measurement(
                f"{table}:{source_id}:{column}:{cell['box_id']}", target,
                axis["sample_id"], axis["assay_id"],
                "peptide_abundance", cell["value"], f"PD2.4 {table}.{column}",
                self.dataset_id, self.dataset_version,
                status="measured" if cell["status"] == "measured" else "missing",
                provenance=dict(
                    source_asset=self.source_asset, native_table=table, native_record=identifiers,
                    native_column=column, native_partition="target",
                    native_quality=quality, native_modifications=row.get("Modifications"),
                    native_modified_sequence=row.get("ModifiedSequence", row.get("AnnotatedSequence")),
                    native_distribution=cell, channel_role=axis["channel_role"],
                    source_column_definition={k: definition[k] for k in
                                              ("ColumnID", "Property_DisplayName", "Property_Description",
                                               "SpecialValueType", "SpecialValueTypeID")},
                    protein_attribution="separate reported memberships and searched-reference candidates",
                    transformation="none"))


def peptide_candidates(peptides, search_reference, *, max_residues=250000, equate_il=True):
    """Stream all placements in an already verified, complete searched FASTA.

    Inputs are a bounded batch of unmodified sequence strings. Every accession
    and repeated occurrence is retained. The optional I/L comparison keeps its
    result distinct from exact residue equality; it does not modify sequences.
    Missing candidates are explicit. This establishes sequence compatibility,
    not detection, unique attribution, or coverage of unsearched canine isoforms.
    The native vendor-reported assignments must be retained separately.
    """
    from ..search_reference import SearchReference

    if not isinstance(search_reference, SearchReference):
        raise TypeError("Peptide matching requires a checksum-verified SearchReference")
    if type(max_residues) is not int or max_residues < 1 or type(equate_il) is not bool:
        raise ValueError("Invalid peptide matching resource limit or I/L policy")
    sequences, residues = set(), 0
    for peptide in peptides:
        if not isinstance(peptide, str) or not re.fullmatch("[A-Z]+", peptide):
            raise ValueError("Peptide must be a nonempty unmodified uppercase sequence")
        if peptide not in sequences:
            residues += len(peptide)
            if residues > max_residues:
                raise ValueError("Peptide matching batch exceeds its explicit residue budget")
            sequences.add(peptide)
    if not sequences:
        return
    try:
        import ahocorasick
    except ImportError:
        raise ImportError("Peptide extraction requires orfdata[canine-extraction]") from None
    normalized = {}
    for peptide in sorted(sequences):
        key = peptide.replace("I", "L") if equate_il else peptide
        normalized.setdefault(key, []).append(peptide)
    matcher = ahocorasick.Automaton()
    for key, originals in sorted(normalized.items()):
        matcher.add_word(key, originals)
    matcher.make_automaton()
    found = set()
    for protein in sorted(search_reference.entries(), key=lambda record: record.target.identifier):
        sequence = protein.sequence.replace("I", "L") if equate_il else protein.sequence
        for end, originals in matcher.iter(sequence):
            for peptide in originals:
                start = end + 1 - len(peptide)
                exact = protein.sequence[start:end + 1] == peptide
                found.add(peptide)
                yield dict(peptide=peptide, accession=protein.target.identifier,
                           start=start, end=end + 1, status="candidate",
                           basis="exact_residues" if exact else "isoleucine_leucine_equivalent")
    for peptide in sorted(sequences - found):
        yield dict(peptide=peptide, accession=None, start=None, end=None,
                   status="no_search_match", basis="exact_or_il_equivalent" if equate_il else "exact_residues")


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
        definition, distribution, boxes = self.distribution(table, column)
        vector = self.vector(table, column, blob)
        for box, cell in zip(boxes, vector.cells()):
            yield dict(column_id=definition["ColumnID"], distribution_id=distribution["ID"],
                       box_id=box["BoxID"], position=box["Position"], label=box["Name"], **cell)

    def distribution(self, table, column):
        """Return the exact native field, numeric map and ordered box definitions."""
        definition = self.columns[(self.tables[table]["DataTypeID"], column)]
        if definition["SpecialValueType"] != "DataDistribution":
            raise ValueError("Column is not a native data distribution")
        distribution = self.maps[definition["SpecialValueTypeID"]]
        boxes = self.boxes[distribution["ID"]]
        return definition, distribution, boxes

    def vector(self, table, column, blob):
        """Normalize a PD vector without interpreting or discarding numeric payload bits."""
        _, distribution, boxes = self.distribution(table, column)
        kind = self.numeric_types[distribution["CustomDataType"]]
        types = {"Int": "int32", "Int64": "int64", "Double": "float64"}
        if kind not in types:
            raise ValueError("Unsupported declared distribution type: " + kind)
        value_type = types[kind]
        width = struct.calcsize(NumericVector.FORMATS[value_type])
        count = len(boxes)
        if blob is None:
            return NumericVector(value_type, count, b"", b"", True)
        if not isinstance(blob, bytes) or len(blob) != (width + 1) * count:
            raise ValueError("Distribution length differs from its declared type and ordered boxes")
        values, presence = bytearray(), bytearray((count + 7) // 8)
        for index in range(count):
            start = index * (width + 1)
            values.extend(blob[start:start + width])
            present = blob[start + width]
            if present not in (0, 1):
                raise ValueError("Invalid native distribution presence flag")
            presence[index // 8] |= present << (index % 8)
        return NumericVector(value_type, count, bytes(values), bytes(presence))

    @staticmethod
    def native_bytes(vector):
        """Reconstruct original PD bytes for independent hash/round-trip checks."""
        if not isinstance(vector, NumericVector):
            raise TypeError("Expected a normalized numeric vector")
        if vector.source_null:
            return None
        width = struct.calcsize(NumericVector.FORMATS[vector.value_type])
        output = bytearray()
        for index in range(vector.count):
            output.extend(vector.values[index * width:(index + 1) * width])
            output.append((vector.present[index // 8] >> (index % 8)) & 1)
        return bytes(output)
