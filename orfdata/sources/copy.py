"""CoPY yeast SRM, label-free estimates and RNA, with native calibration scope."""

import csv
import gzip
import hashlib
import math
import re
from collections import Counter

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample, TargetMapping
from ..pairing import AssayPair, TargetRelation
from ..reference import ScopedTarget
from ..reference_mapping import ReferenceMapping

HEADERS = {
    "ProteinQuant": ("Function", "Protein", "Name", "Type", "Peptide X", "Quant 1", "Quant 2", "Quant 3", "Quant 4", "Median Quant X", "Peptide Y", "Quant 1", "Quant 2", "Quant 3", "Quant 4", "Median Quant Y", "Log2 Ratio Peptide X / Peptide Y", "Final Quant", "MAD", "rSD", "rCV"),
    "PeptideQuant": ("CopyCAT", "FunctionalGroup", "Protein", "ProteinName", "ProteinClassification", "PeptideSequence", "Index within CoPYCAT", "PeptideClassification", "Quant1", "Quant2", "Quant3", "Quant4", "MedianQuant", "MAD", "rSD", "rCV"),
    "RNA": ("Transcript", "FPKM", "cpc"),
    "LabelFree": ("Name", "Qexactive (Copies Per Cell)", "SAX (Copies per Cell)"),
}
RNA_COLUMNS = ("tracking_id", "class_code", "nearest_ref_id", "gene_id", "gene_short_name", "tss_id", "locus", "length", "coverage", "FPKM", "FPKM_conf_lo", "FPKM_conf_hi", "FPKM_status")
RETARGETED = "Not Used. Protein re-targetted by another QconCAT"


def _numeric(value):
    if isinstance(value, bool):
        raise ValueError("CoPY abundance must be a number, not a boolean")
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError("CoPY abundance must be finite and nonnegative")
    return number


def _quantity(value):
    if value is None:
        return dict(value=None, status="missing")
    if value in ("-", "--", RETARGETED):
        return dict(value=None, status="not_measured")
    if isinstance(value, str) and value.startswith("<"):
        return dict(value=None, status="below_detection", detection_limit=_numeric(value[1:]))
    return dict(value=_numeric(value), status="measured")


def _workbook(path, expected_counts):
    from openpyxl import load_workbook
    workbook = load_workbook(path, read_only=True, data_only=True)
    tables = {}
    try:
        for sheet, header in HEADERS.items():
            rows = tuple(workbook[sheet].values)
            if not rows or rows[0] != header:
                raise ValueError(f"Unexpected CoPY {sheet} columns")
            tables[sheet] = tuple((line, row) for line, row in enumerate(rows[1:], 2)
                                  if any(value is not None for value in row))
            if len(tables[sheet]) != expected_counts[sheet]:
                raise ValueError(f"CoPY {sheet} row count changed")
    finally:
        workbook.close()
    for sheet, index in (("ProteinQuant", 1), ("RNA", 0), ("LabelFree", 0)):
        identifiers = [row[index] for _, row in tables[sheet]]
        if any(not isinstance(value, str) or not value for value in identifiers):
            raise ValueError(f"CoPY {sheet} target IDs must be nonempty strings")
        if sheet != "LabelFree" and len(set(identifiers)) != len(identifiers):
            raise ValueError(f"CoPY {sheet} target IDs must be nonempty and unique")
    return tables


def _standards(path):
    sequences, identifier = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split("|")[0]
                if identifier in sequences or not identifier:
                    raise ValueError("Duplicate/empty QconCAT standard ID")
                sequences[identifier] = ""
            elif identifier is None or not re.fullmatch("[A-Za-z]+", line):
                raise ValueError("Invalid QconCAT standard FASTA")
            else:
                sequences[identifier] += line.upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError("Empty QconCAT standard")
    return sequences


def _standard_record(row, standards):
    identifier = "OxidativeMopUp" if row[0] == "OxMop" else f"CopyCat{int(row[0]):03d}"
    sequence = standards.get(identifier)
    peptide = row[5]
    if not isinstance(peptide, str) or not re.fullmatch("[A-Z]+", peptide):
        raise ValueError("Invalid CoPY peptide sequence")
    positions = tuple(match.start() + 1 for match in re.finditer(f"(?={peptide})", sequence or ""))
    return dict(standard_id=identifier, source_asset="qconcat.fasta", peptide=peptide,
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest() if sequence else None,
                peptide_positions_1based=positions,
                validation="exact_subsequence" if positions else "peptide_not_in_standard" if sequence else "standard_missing")


def _rna_rows(path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != RNA_COLUMNS:
            raise ValueError("Unexpected CoPY Cufflinks RNA columns")
        seen = set()
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()) or not row["tracking_id"]:
                raise ValueError("Incomplete CoPY RNA row")
            if row["tracking_id"] in seen:
                raise ValueError("Duplicate RNA tracking ID within a replicate")
            seen.add(row["tracking_id"])
            if row["FPKM_status"] not in ("OK", "LOWDATA", "HIDATA", "FAIL"):
                raise ValueError("Unknown Cufflinks quantification status")
            for name in ("FPKM", "FPKM_conf_lo", "FPKM_conf_hi"):
                _numeric(row[name])
            yield line, row


def import_copy(inputs, output, manifest):
    """Transform DataCache-verified source files without inferring replicate pairs."""
    if manifest.processing.get("importer") != "copy-yeast-v1":
        raise ValueError("Not a CoPY source manifest")
    references = {reference.key: reference for reference in manifest.references}
    roles = manifest.processing["reference_roles"]
    reference = references[roles["workbook"]]
    rna_reference = references[roles["rna_replicates"]]
    if any(value.taxon_id != 4932 for value in (reference, rna_reference)):
        raise ValueError("CoPY sources must remain yeast-scoped")
    assets = {asset.asset_id: asset for asset in manifest.sources}
    if (reference.source_sha256 != assets["workbook.xlsx"].sha256
            or rna_reference.source_sha256 != assets["geo-series.soft"].sha256):
        raise ValueError("CoPY native references disagree with their pinned sources")
    tables = _workbook(inputs["workbook.xlsx"], manifest.coverage["sheet_rows"])
    standards = _standards(inputs["qconcat.fasta"])
    if len(standards) != manifest.coverage["standard_sequences"]:
        raise ValueError("QconCAT standard count changed")
    standard_records = {line: _standard_record(row, standards) for line, row in tables["PeptideQuant"]}
    if Counter(value["validation"] for value in standard_records.values()) != dict(manifest.coverage["standard_validation"]):
        raise ValueError("QconCAT peptide/standard reconciliation changed")
    if Counter(row[3] for _, row in tables["ProteinQuant"]) != dict(manifest.coverage["protein_classes"]):
        raise ValueError("Protein classification counts changed")
    rna_inputs = manifest.processing["rna_inputs"]
    for entry in rna_inputs:
        rows = tuple(_rna_rows(inputs[entry["asset_id"]]))
        if len(rows) != entry["rows"] or Counter(row["FPKM_status"] for _, row in rows) != dict(entry["statuses"]):
            raise ValueError("RNA replicate count/quality reconciliation changed")
    total_fpkm = sum(_numeric(row[1]) for _, row in tables["RNA"])
    if total_fpkm <= 0:
        raise ValueError("Aggregate RNA calibration requires a positive FPKM total")
    for _, row in tables["RNA"]:
        if not math.isclose(_numeric(row[2]), _numeric(row[1]) / total_fpkm * 60000, rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError("Author RNA copies/cell disagree with the documented 60,000-mRNA calibration")

    sample_attributes = dict(strain="EUROSCARF Y11335 BY4742", genotype="MATalpha his3delta1 leu2delta0 lys2delta0 ura3delta0 YJL088w::kanMX4",
                             medium="defined C-limiting F1; glucose 10 g/L", dilution_rate_per_hour=0.1,
                             cross_assay_replicate_correspondence="unresolved")
    cohort = "CoPY:chemostat-cohort"
    samples = [Sample("cohort", reference, biological_unit_id=cohort, condition="glucose-limited chemostat",
                      attributes=dict(sample_attributes, level="source cohort aggregate", biological_cultures=4),
                      provenance={"source_asset": "workbook.xlsx", "methods_asset": "paper.xml"})]
    for replicate in range(1, 5):
        samples.append(Sample(f"SRM:Quant{replicate}", reference, condition="glucose-limited chemostat",
                              biological_replicate=str(replicate), attributes=sample_attributes,
                              provenance={"source_asset": "workbook.xlsx", "column": f"Quant{replicate}"}))
    for entry in rna_inputs:
        samples.append(Sample(entry["sample_id"], rna_reference, condition="glucose-limited chemostat",
                              biological_replicate=entry["biological_replicate"], attributes=sample_attributes,
                              provenance={"source_asset": entry["metadata_asset"], "source_title": entry["title"]}))

    srm_calibration = dict(method="protein-specific stable-isotope QconCAT SRM; heavy/light ratio and cell count",
                           standards_asset="qconcat.fasta", source_standard_sequences=len(standards),
                           sample_cell_equivalents_per_injection=200000,
                           light_heavy_ratio_per_transition="not supplied in processed workbook")
    srm_acquisition = dict(instrument="Waters nanoAcquity UPLC / Xevo TQMS tandem quadrupole",
                           transitions_per_peptide=3, transition_identities="not supplied in processed workbook",
                           scheduled_window_minutes=3, source_srm_accession="PASS00717")
    assays = [Assay("SRM:aggregate", "cohort", "targeted_ms", reference, calibration=srm_calibration,
                    attributes=dict(srm_acquisition, level="source aggregate"),
                    provenance={"source_asset": "workbook.xlsx"}),
              Assay("RNA:aggregate", "cohort", "rna_seq", reference,
                    normalization={"method": "source counts aggregated over four biological replicates before FPKM estimation"},
                    calibration={"RNA_copies_per_cell": "derived assumption", "assumed_total_mrna_per_cell": 60000,
                                 "source_fpkm_sum": total_fpkm}, provenance={"source_asset": "workbook.xlsx"}),
              Assay("LabelFree:Qexactive", "cohort", "digest_ms", reference,
                    search_reference="UniProt yeast 6560 proteins plus P00489 spike; release unspecified",
                    calibration={"method": "Top3 with glycogen phosphorylase B internal standard", "standard_uniprot": "P00489",
                                 "standard_fmol": 50, "cell_equivalents": 400000},
                    attributes={"replicate_correspondence": "unresolved", "source_ms_accession": "PXD002694"},
                    provenance={"source_asset": "workbook.xlsx", "methods_asset": "methods.pdf", "pages": [5, 6]}),
              Assay("LabelFree:SAX", "cohort", "digest_ms", reference,
                    search_reference="UniProt yeast 6560 proteins; release unspecified",
                    calibration={"method": "Top3 scaled to assumed total protein", "assumed_total_protein_molecules_per_cell": 60000000},
                    attributes={"replicate_correspondence": "unresolved", "source_ms_accession": "PXD002694"},
                    provenance={"source_asset": "workbook.xlsx", "methods_asset": "methods.pdf", "pages": [6, 7]})]
    for replicate in range(1, 5):
        assays.append(Assay(f"SRM:Quant{replicate}", f"SRM:Quant{replicate}", "targeted_ms", reference,
                            calibration=srm_calibration, attributes=dict(srm_acquisition, biological_replicate_column=f"Quant{replicate}"),
                            provenance={"source_asset": "workbook.xlsx"}))
    for entry in rna_inputs:
        assays.append(Assay(entry["sample_id"], entry["sample_id"], "rna_seq", rna_reference,
                            search_reference="Bowtie yeast prebuilt genome and SGD GTF; exact releases unspecified",
                            normalization={"method": "source Cufflinks 2.0.0 FPKM", "quantification_status_required": "OK"},
                            provenance={"source_asset": entry["asset_id"]}))

    def protein_target(identifier, *, label_free=False):
        grouped = label_free or "_" in identifier
        return ScopedTarget(reference, "protein_group" if grouped else "protein", identifier,
                            "CoPY:LabelFree" if label_free else "CoPY:protein-group" if grouped else "sgd")

    relations = []
    for line, row in tables["ProteinQuant"]:
        target = protein_target(row[1])
        for gene in row[1].split("_"):
            grouped = target.kind == "protein_group"
            relations.append(TargetRelation(f"ProteinQuant:{line}:{gene}", manifest.dataset_id, manifest.version,
                                            ScopedTarget(reference, "gene", gene, "sgd"), target,
                                            "gene_protein_group_context" if grouped else "canonical_translation",
                                            "many_to_one" if grouped else "one_to_one",
                                            {"source_asset": "workbook.xlsx", "source_systematic_id": row[1],
                                             "method": "source SGD systematic gene/protein labels", "sequence_validation": "not performed"}))
    pair = AssayPair("RNA-SRM:cohort", reference, manifest.dataset_id, manifest.version,
                     "RNA:aggregate", "SRM:aggregate", "condition_aggregate",
                     {"source_asset": "paper.xml", "source": "same-cell transcriptome and protein study, aggregate tables",
                      "granularity": "four-culture source cohort", "individual_replicate_crosswalk": "unresolved"}, "concurrent")

    def mappings():
        for line, row in tables["PeptideQuant"]:
            yield TargetMapping(f"PeptideQuant:{line}", "CoPY:PeptideQuant-row", f"row:{line}", reference,
                                (protein_target(row[2]),), "unique", "source-reported target label",
                                evidence={"resolution": "assignment to one reported target, not proof of unique peptide/protein attribution",
                                          "peptide_sequence": row[5], "standard": standard_records[line]},
                                provenance={"source_asset": "workbook.xlsx", "sheet": "PeptideQuant", "row": line})

    def reference_mappings():
        sources = {protein_target(row[1]) for _, row in tables["ProteinQuant"] if "_" in row[1]}
        sources.update(protein_target(row[0], label_free=True) for _, row in tables["LabelFree"])
        for source in sorted(sources, key=lambda target: target.key):
            members = tuple(ScopedTarget(reference, "protein", value, "sgd") for value in source.identifier.split("_"))
            yield ReferenceMapping(f"{source.namespace}:{source.identifier}", manifest.dataset_id, manifest.version,
                                   source, reference, members, "ambiguous" if len(members) > 1 else "unresolved",
                                   "source-reported systematic identifiers",
                                   {"sequence_assignment": "not tested", "membership": "reported labels; no full-background uniqueness"},
                                   reason="Group attribution retained; individual abundance unresolved",
                                   provenance={"source_asset": "workbook.xlsx"})

    def observation(identifier, target, sample, assay, quantity, value, *, provenance, unit="copies/cell", **kwargs):
        return Measurement(identifier, target, sample, assay, quantity, unit=unit,
                           dataset_id=manifest.dataset_id, dataset_version=manifest.version,
                           provenance=provenance, **_quantity(value), **kwargs)

    def uncertainty(value):
        if isinstance(value, (int, float)):
            return dict(uncertainty=_numeric(value), uncertainty_type="median_absolute_deviation")
        return {}

    def measurements():
        for line, row in tables["ProteinQuant"]:
            provenance = dict(source_asset="workbook.xlsx", sheet="ProteinQuant", row=line, column="Final Quant",
                              classification=row[3], native_row=row, source_headers=HEADERS["ProteinQuant"],
                              estimate="author-reported calibrated protein summary; not raw instrument intensity")
            yield observation(f"ProteinQuant:{line}", protein_target(row[1]), "cohort", "SRM:aggregate",
                              "protein_abundance", row[17], provenance=provenance, **uncertainty(row[18]))
        for line, row in tables["PeptideQuant"]:
            target = ScopedTarget(reference, "peptide", row[5], "amino-acid-sequence")
            provenance = dict(source_asset="workbook.xlsx", sheet="PeptideQuant", row=line,
                              native_row=row, source_headers=HEADERS["PeptideQuant"],
                              protein_label=row[2], classification=row[7], standard=standard_records[line])
            for index in range(4):
                assay = f"SRM:Quant{index + 1}"
                yield observation(f"PeptideQuant:{line}:Quant{index + 1}", target, assay, assay,
                                  "peptide_abundance", row[8 + index], provenance=dict(provenance, column=f"Quant{index + 1}"))
            yield observation(f"PeptideQuant:{line}:Median", target, "cohort", "SRM:aggregate", "peptide_abundance",
                              row[12], provenance=dict(provenance, column="MedianQuant"), **uncertainty(row[13]))
        for line, row in tables["RNA"]:
            target = ScopedTarget(reference, "gene", row[0], "sgd")
            provenance = dict(source_asset="workbook.xlsx", sheet="RNA", row=line, native_row=row,
                              source_feature_column="Transcript", target_resolution="SGD systematic gene identifier")
            identifier = f"RNA:{line}:FPKM"
            yield observation(identifier, target, "cohort", "RNA:aggregate", "rna_abundance", row[1],
                              unit="FPKM", provenance=dict(provenance, column="FPKM"))
            yield observation(f"RNA:{line}:cpc", target, "cohort", "RNA:aggregate", "rna_abundance", row[2],
                              provenance=dict(provenance, column="cpc"), origin="derived", input_measurements=(identifier,),
                              derivation={"method": "author FPKM / aggregate FPKM sum * assumed total mRNA",
                                          "normalization_fpkm_sum": total_fpkm, "assumed_total_mrna_per_cell": 60000,
                                          "normalization_source": "workbook.xlsx/RNA; all source rows",
                                          "absolute_RNA_directly_measured": False})
        for line, row in tables["LabelFree"]:
            for index, method in ((1, "Qexactive"), (2, "SAX")):
                yield observation(f"LabelFree:{line}:{method}", protein_target(row[0], label_free=True),
                                  "cohort", f"LabelFree:{method}", "protein_abundance", row[index],
                                  provenance=dict(source_asset="workbook.xlsx", sheet="LabelFree", row=line,
                                                  column=HEADERS["LabelFree"][index], native_row=row,
                                                  estimate="source calibrated label-free estimate; distinct from targeted SRM"))
        for entry in rna_inputs:
            sample = entry["sample_id"]
            for line, row in _rna_rows(inputs[entry["asset_id"]]):
                accepted = row["FPKM_status"] == "OK"
                target = ScopedTarget(rna_reference, "transcript", row["tracking_id"], "CoPY:Cufflinks")
                yield Measurement(f"{sample}:{line}", target, sample, sample, "rna_abundance",
                                  _numeric(row["FPKM"]) if accepted else None, "FPKM", manifest.dataset_id, manifest.version,
                                  status="measured" if accepted else "not_measured",
                                  provenance=dict(source_asset=entry["asset_id"], line=line, native_row=row,
                                                  confidence_interval={"level": 0.95, "lower": _numeric(row["FPKM_conf_lo"]),
                                                                       "upper": _numeric(row["FPKM_conf_hi"])},
                                                  unresolved_reason=None if accepted else "source_quantification_" + row["FPKM_status"]))

    return QuantitativeDataset.create(output, manifest, samples=samples, assays=assays,
                                      measurements=measurements(), mappings=mappings(), pairings=(pair,),
                                      relations=relations, reference_mappings=reference_mappings())
