"""Paired canine mammary-stroma measurements from GSE135183 / PXD023023.

Same archived specimens do not imply identical aliquots. Native protein groups,
peptide rows, controls and unresolved search references remain source-scoped.
"""

import csv
import gzip
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample, TargetMapping
from ..pairing import AssayPair
from ..reference import ReferenceContext, ScopedTarget

IMPORTER = "canine-mammary-stroma-v1"
SOFT_FIELDS = {
    "geo_accession",
    "title",
    "organism_ch1",
    "taxid_ch1",
    "characteristics_ch1",
    "description",
    "data_processing",
    "relation",
    "instrument_model",
}
QUANTITIES = {
    "proteinGroups.txt": (
        ("Intensity ", "intensity"),
        ("LFQ intensity ", "LFQ intensity"),
        ("iBAQ ", "iBAQ"),
    ),
    "peptides.txt": (("Intensity ", "intensity"), ("LFQ intensity ", "LFQ intensity")),
}


def _identity_text(value):
    value = re.sub(r"[^a-z0-9]", "", value.lower())
    return {"bassett": "basset", "viszla": "vizsla", "nd": "notdisclosed"}.get(value, value)


def clinical_cases(path):
    """Read factual Table 1 cells, preserving source spellings and case IDs."""
    table = ET.parse(path).find(".//table-wrap[@id='tbl0001']")
    if table is None:
        raise ValueError("Clinical Table 1 is missing")
    rows = [
        [" ".join("".join(cell.itertext()).split()) for cell in row] for row in table.iter("tr")
    ]
    expected = [
        "Case #",
        "Gender",
        "Breed",
        "Age [years]",
        "Simple carcinoma",
        "Sub type of simple carcinoma",
    ]
    if not rows or rows[0] != expected:
        raise ValueError("Clinical Table 1 schema changed")
    result = {}
    for cells in rows[1:]:
        if len(cells) != 6 or not cells[0].isdigit() or cells[0] in result:
            raise ValueError("Invalid/duplicate clinical case")
        result[cells[0]] = dict(zip(expected, cells))
    return result


def geo_samples(path):
    samples, current = {}, None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line, text in enumerate(handle, 1):
            key, separator, value = text.rstrip("\r\n").partition(" = ")
            if not separator:
                continue
            if key == "^SAMPLE":
                if value in samples:
                    raise ValueError("Duplicate GEO sample")
                current = dict(accession=value, source_line=line, fields={})
                samples[value] = current
            elif key.startswith("^"):
                current = None
            elif current is not None and key.startswith("!Sample_") and key[8:] in SOFT_FIELDS:
                current["fields"].setdefault(key[8:], []).append(value)
    result = {}
    for accession, sample in samples.items():
        fields = sample["fields"]
        if fields.get("geo_accession") != [accession] or fields.get("taxid_ch1") != ["9615"]:
            raise ValueError("GEO identity/taxon mismatch")
        if len(fields.get("title", ())) != 1:
            raise ValueError("Missing GEO title")
        match = re.search(r"\b(\d+)_(norm|tumour)$", fields["title"][0], re.I)
        if not match:
            raise ValueError("Unrecognized source sample title")
        case, condition = match.groups()
        column = case + ("_Norm" if condition.lower() == "norm" else "_Tumour")
        if column in result:
            raise ValueError("Duplicate RNA matrix sample")
        traits = {}
        for item in fields["characteristics_ch1"]:
            key, separator, value = item.partition(": ")
            if not separator or key in traits:
                raise ValueError("Invalid/duplicate sample characteristic")
            traits[key] = value
        if traits.get("Sex", "").lower() != "female":
            raise ValueError("Unexpected cohort sex")
        result[column] = dict(sample, case=case, condition=condition.lower(), traits=traits)
    return result


def specimen_crosswalk(clinical, rna):
    """Unique metadata agreement within the author-declared shared cohort."""
    by_case = {}
    for row in rna.values():
        traits = row["traits"]
        key = tuple(
            _identity_text(traits[k])
            for k in ("breed", "age of patient (years)", "subtype of simple carcinoma")
        )
        if row["case"] in by_case and by_case[row["case"]] != key:
            raise ValueError("RNA compartments disagree about donor identity")
        by_case[row["case"]] = key
    crosswalk = {}
    for case, row in clinical.items():
        key = tuple(
            _identity_text(row[k]) for k in ("Breed", "Age [years]", "Sub type of simple carcinoma")
        )
        candidates = [identifier for identifier, identity in by_case.items() if identity == key]
        if len(candidates) != 1 or candidates[0] in crosswalk.values():
            raise ValueError("Clinical/RNA case crosswalk is absent or ambiguous")
        crosswalk[case] = candidates[0]
    return crosswalk


def _table(path, expected):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = tuple(reader.fieldnames or ())
        digest = hashlib.sha256(json.dumps(columns, separators=(",", ":")).encode()).hexdigest()
        if len(set(columns)) != len(columns) or digest != expected["header_sha256"]:
            raise ValueError("Native MS table schema changed")
        rows, seen = [], set()
        for line, row in enumerate(reader, 2):
            if (
                None in row
                or any(value is None for value in row.values())
                or not row["id"].isdigit()
            ):
                raise ValueError("Malformed native MS row")
            if row["id"] in seen:
                raise ValueError("Duplicate native MS row ID")
            seen.add(row["id"])
            if any(row[k] not in ("", "+") for k in ("Reverse", "Potential contaminant")):
                raise ValueError("Unknown native control/decoy flag")
            rows.append((line, row))
    if len(rows) != expected["rows"]:
        raise ValueError("Native MS row-count reconciliation failed")
    if any(
        dict(Counter(row[key] for _, row in rows)) != dict(counts)
        for key, counts in expected["qc_counts"].items()
    ):
        raise ValueError("Native MS control/decoy counts changed")
    return columns, rows


def _value(raw, *, zero_missing=False, count=False):
    if zero_missing and raw == "NaN":
        return None, "missing"
    value = float(raw)
    if not math.isfinite(value) or value < 0 or (count and not value.is_integer()):
        raise ValueError(
            "Native quantities must be finite nonnegative numbers in the reported scale"
        )
    return (
        (None, "missing")
        if zero_missing and value == 0
        else (int(value) if count else value, "measured")
    )


def import_stroma(inputs, output, manifest):
    """Import pinned processed files through the common quantitative dataset API."""
    if manifest.processing.get("importer") != IMPORTER or len(manifest.references) != 2:
        raise ValueError("Not a paired mammary-stroma source manifest")
    assets = {a.asset_id: a for a in manifest.sources}
    if set(inputs) != set(assets):
        raise ValueError("All pinned study inputs are required")
    for name, path in inputs.items():
        with open(path, "rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != assets[name].sha256:
                raise ValueError("Study input checksum mismatch: " + name)
    rna_reference, protein_reference = manifest.references
    if tuple(manifest.references) != (
        ReferenceContext(
            9615, "dog", "CanFam3.1", "ensembl", source_sha256=assets["geo.soft.gz"].sha256
        ),
        ReferenceContext(
            9615, "dog", None, "uniprot", source_sha256=assets["proteinGroups.txt"].sha256
        ),
    ):
        raise ValueError("Native incomplete reference identities disagree with the study sources")
    rna = geo_samples(inputs["geo.soft.gz"])
    clinical = clinical_cases(inputs["paper.xml"])
    crosswalk = specimen_crosswalk(clinical, rna)
    if crosswalk != dict(manifest.processing["specimen_crosswalk"]):
        raise ValueError("Specimen crosswalk differs from the audited cohort")
    if len(rna) != manifest.coverage["rna_samples"]:
        raise ValueError("RNA sample count changed")
    if any(
        {row["condition"] for row in rna.values() if row["case"] == case} != {"norm", "tumour"}
        for case in {row["case"] for row in rna.values()}
    ):
        raise ValueError("Every RNA case must retain both stromal compartments")
    tables = {name: _table(inputs[name], manifest.coverage["tables"][name]) for name in QUANTITIES}
    runs = sorted(
        c.removeprefix("LFQ intensity ")
        for c in tables["proteinGroups.txt"][0]
        if c.startswith("LFQ intensity ")
    )
    if len(runs) != manifest.coverage["protein_samples"]:
        raise ValueError("Protein sample count changed")
    for name, (columns, _) in tables.items():
        for prefix, _ in QUANTITIES[name]:
            if sorted(c[len(prefix) :] for c in columns if c.startswith(prefix)) != runs:
                raise ValueError("Native quantitative columns disagree about sample coverage")
    samples, assays, protein_samples = [], [], {}
    condition_names = {
        "norm": "normal stroma in a cancer-bearing mammary specimen",
        "tumour": "cancer-associated mammary stroma",
    }
    for column, row in sorted(rna.items()):
        identifier = "rna:" + row["accession"]
        samples.append(
            Sample(
                identifier,
                rna_reference,
                biological_unit_id="GSE135183:case:" + row["case"],
                tissue="mammary gland",
                cell_type="stroma",
                condition=condition_names[row["condition"]],
                attributes=dict(
                    source_case=row["case"], pooling="single", source_metadata=row["fields"]
                ),
                provenance=dict(
                    source_asset="geo.soft.gz", source_line=row["source_line"], matrix_column=column
                ),
            )
        )
        assays.append(
            Assay(
                identifier,
                identifier,
                "rna_seq",
                rna_reference,
                search_reference="Ensembl CanFam3.1 transcripts; exact release and sequence bytes unresolved",
                normalization=dict(
                    method="Kallisto 0.44.0 estimated counts, gene summarization by tximport 1.8.0",
                    source_scale="integer estimated gene counts; no TPM conversion",
                ),
                attributes=dict(source_accession=row["accession"]),
                provenance=dict(source_asset="geo.soft.gz"),
            )
        )
    for run in runs:
        match = re.fullmatch(r"\d+_\d+_S\d+_(\d+)_([NC])_(Control|Diseased)", run)
        if not match or match.group(1) not in crosswalk:
            raise ValueError("Unrecognized source proteomics run")
        case, label, state = match.groups()
        if (label, state) not in (("N", "Control"), ("C", "Diseased")):
            raise ValueError("Conflicting native tissue-compartment labels")
        condition = "norm" if label == "N" else "tumour"
        column = crosswalk[case] + ("_Norm" if condition == "norm" else "_Tumour")
        if column not in rna or column in protein_samples:
            raise ValueError("Protein/RNA compartment crosswalk is missing or duplicated")
        identifier = "protein:" + run
        protein_samples[column] = identifier
        samples.append(
            Sample(
                identifier,
                protein_reference,
                biological_unit_id="GSE135183:case:" + crosswalk[case],
                tissue="mammary gland",
                cell_type="stroma",
                condition=condition_names[condition],
                attributes=dict(
                    source_case=case,
                    pooling="single",
                    clinical=clinical[case],
                    preparation="Separate LCM preparation from the same archived FFPE specimen",
                ),
                provenance=dict(source_asset="paper.xml", table="tbl0001", source_case=case),
            )
        )
        assays.append(
            Assay(
                identifier,
                identifier,
                "digest_ms",
                protein_reference,
                run_id=run,
                search_reference="Canonical UniProt dog proteome dated 2019-07-29, reversed decoys and contaminants; FASTA unavailable",
                normalization=dict(
                    method="Source MaxQuant 1.6.2.3 quantities; no re-normalization or imputation"
                ),
                attributes=dict(
                    search_engine="Andromeda",
                    peptide_fdr=0.01,
                    protein_fdr=0.05,
                    match_between_runs_minutes=2,
                ),
                provenance=dict(source_asset="paper.xml"),
            )
        )
    pairs = []
    inverse = {value: key for key, value in crosswalk.items()}
    for column, row in sorted(rna.items()):
        pairs.append(
            AssayPair(
                "rna-protein:" + column,
                rna_reference,
                manifest.dataset_id,
                manifest.version,
                "rna:" + row["accession"],
                protein_samples.get(column),
                "same_specimen",
                dict(
                    basis="Author-declared shared FFPE cohort; unique breed/age/histology reconciliation",
                    source_assets=["paper.xml", "geo.soft.gz"],
                    source_table="tbl0001",
                    rna_case=row["case"],
                    protein_case=inverse.get(row["case"]),
                    crosswalk_status="metadata_derived; not author-provided ID mapping",
                    preparation="Separate LCM preparations; aliquot and section adjacency not established",
                    missing_mate_reason=None
                    if column in protein_samples
                    else "No matching proteomics case; study reports case 7 excluded because its block was unavailable",
                ),
                "concurrent" if column in protein_samples else "unknown",
                right_reference=protein_reference,
            )
        )

    def target(name, identifier):
        return ScopedTarget(
            protein_reference,
            "protein_group" if name == "proteinGroups.txt" else "peptide",
            identifier,
            "PXD023023:" + name,
        )

    groups = {row["id"] for _, row in tables["proteinGroups.txt"][1]}

    def mappings():
        for name, (_, rows) in tables.items():
            for line, row in rows:
                group_ids = (
                    row.get("Protein group IDs", "").split(";") if name == "peptides.txt" else []
                )
                members = tuple(
                    target("proteinGroups.txt", group) for group in group_ids if group in groups
                )
                unresolved = any(group not in groups for group in group_ids)
                yield TargetMapping(
                    name + ":" + row["id"],
                    "PXD023023:" + name,
                    row["id"],
                    protein_reference,
                    members,
                    "unresolved"
                    if name == "proteinGroups.txt" or unresolved
                    else "unique"
                    if len(members) == 1
                    else "ambiguous"
                    if members
                    else "unmapped",
                    "Native peptide-to-group labels; protein reference assignment not established",
                    "Complete searched FASTA unavailable; no unique protein attribution"
                    if name == "proteinGroups.txt"
                    else "Source group endpoint absent"
                    if unresolved
                    else None,
                    evidence=dict(
                        unique_protein_detection=False,
                        source_partition="decoy"
                        if row["Reverse"]
                        else "contaminant"
                        if row["Potential contaminant"]
                        else "target",
                        sequence_validation="not performed; no complete searched FASTA",
                    ),
                    provenance=dict(source_asset=name, line=line, native_record=row),
                )

    def measurements():
        with gzip.open(inputs["rna.csv.gz"], "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = tuple(reader.fieldnames or ())
            if (
                not columns
                or columns[0] != ""
                or len(set(columns)) != len(columns)
                or set(columns[1:]) != set(rna)
            ):
                raise ValueError("RNA matrix and GEO sample identities differ")
            seen = set()
            for line, row in enumerate(reader, 2):
                gene = row[""]
                if (
                    None in row
                    or any(value is None for value in row.values())
                    or not re.fullmatch(r"ENSCAFG\d+", gene)
                    or gene in seen
                ):
                    raise ValueError("RNA gene identity/row is invalid or duplicated")
                seen.add(gene)
                for column in columns[1:]:
                    value, status = _value(row[column], count=True)
                    assay = "rna:" + rna[column]["accession"]
                    yield Measurement(
                        f"rna:{gene}:{column}",
                        ScopedTarget(rna_reference, "gene", gene, "ensembl"),
                        assay,
                        assay,
                        "rna_abundance",
                        value,
                        "estimated_count",
                        manifest.dataset_id,
                        manifest.version,
                        status=status,
                        provenance=dict(
                            source_asset="rna.csv.gz",
                            line=line,
                            column=column,
                            native_value=row[column],
                        ),
                    )
            if len(seen) != manifest.coverage["rna_genes"]:
                raise ValueError("RNA gene-count reconciliation failed")
        for name, (_, rows) in tables.items():
            counts = Counter()
            for line, row in rows:
                qc = {key: row[key] for key in ("Reverse", "Potential contaminant")}
                if "Only identified by site" in row:
                    qc["Only identified by site"] = row["Only identified by site"]
                for prefix, unit in QUANTITIES[name]:
                    for run in runs:
                        column = prefix + run
                        value, status = _value(row[column], zero_missing=True)
                        counts[unit + ":" + status] += 1
                        identifier = "protein:" + run
                        yield Measurement(
                            f"{name}:{row['id']}:{column}",
                            target(name, row["id"]),
                            identifier,
                            identifier,
                            "protein_abundance"
                            if name == "proteinGroups.txt"
                            else "peptide_abundance",
                            value,
                            unit,
                            manifest.dataset_id,
                            manifest.version,
                            status=status,
                            provenance=dict(
                                source_asset=name,
                                line=line,
                                column=column,
                                native_value=row[column],
                                source_qc=qc,
                                missing_reason=(
                                    "native NaN"
                                    if row[column] == "NaN"
                                    else "zero intensity sentinel"
                                )
                                if status == "missing"
                                else None,
                                identification_type=row["Identification type " + run],
                            ),
                        )
            if dict(counts) != dict(manifest.coverage["tables"][name]["measurement_statuses"]):
                raise ValueError("MS quantity/missingness reconciliation failed")

    return QuantitativeDataset.create(
        output,
        manifest,
        samples=samples,
        assays=assays,
        pairings=pairs,
        mappings=mappings(),
        measurements=measurements(),
    )
