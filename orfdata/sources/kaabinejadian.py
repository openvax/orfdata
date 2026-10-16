"""DLA class I ligands from dog tumours and DLA-88 monoallelic human cells (Kaabinejadian et al. 2026).

PRIDE PXD074485 deposits PEAKS Studio 11.5 peptide tables for eight
immunoprecipitations. Each dog tumour was searched against four protein
databases; each table is one search of the same spectra, so the same
peptide in two of a sample's tables is one observation, not two.
"""

import csv
import io
import math
import re
import zipfile
from collections import Counter

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample, TargetMapping
from ..reference import ScopedTarget

COLUMNS = ("Peptide", "-10LgP", "Mass", "Tag Length", "CAA (%)", "Length", "ppm", "m/z", "z", "RT",
           "Delta RT", "MS2 Correlation", "Scan", "Source File", "Area", "Intensity", "#Spec", "#Spec",
           "#Feature", "#Feature", "PTM", "AScore", "Positional Confidence", "tag(>=0.0%)", "Accession",
           "Found By")
FOUND_BY = {
    "DB Search": "PEAKS database search: the peptide occurs in each listed protein",
    "Homolog": "PEAKS SPIDER homology search: the peptide differs from the listed proteins by substitutions",
    "DeepNovo": "PEAKS DeepNovo de novo sequence: the listed proteins need not contain it",
}
# Biognosys iRT kit peptides (Escher et al. 2012), spiked into every run. PEAKS also reports
# fragments of them, and searches whose database lacks the kit report them as other proteins.
IRT_PEPTIDES = ("LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "GAGSSEPVTGLDAK",
                "TPVISGGPYEYR", "VEATFGVDESNAK", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                "LFLQFGAQGSPFLK")
SDRF_SAMPLE = ("source name", "characteristics[organism]", "characteristics[organism part]",
               "characteristics[disease]", "characteristics[cell line]", "characteristics[age]",
               "characteristics[sex]", "characteristics[individual]", "characteristics[mhc typing]")
# Breeds and DLA typing of the dogs, from the article's Table 1; the deposit's SDRF omits breed.
BREEDS = {"Lola": "Black Russian Terrier", "163828A": "English Cocker Spaniel",
          "Lily": "English Cocker Spaniel", "Bogey": "American Cocker Spaniel"}


def ligands_manifest():
    """Return the pinned PXD074485 release, without downloading anything."""
    from ..datasets import DatasetRegistry

    return DatasetRegistry().info("kaabinejadian-dog-dla-ligands")


def build_ligands(*, cache_root=None, source_paths=None, force=False, download_options=None):
    """Explicitly rebuild the release from its pinned PRIDE files."""
    from ..datasets import build_dataset

    manifest = ligands_manifest()
    return build_dataset(manifest, lambda inputs, output: import_ligands(inputs, output, manifest),
                         cache_root=cache_root, source_paths=source_paths, force=force,
                         download_options=download_options)


def peptide_sequence(peptide):
    """The unmodified sequence of a PEAKS peptide such as ``M(+15.99)PEPTIDE``."""
    sequence = re.sub(r"\([+-][0-9.]+\)", "", peptide)
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
        raise ValueError(f"Unexpected PEAKS peptide {peptide!r}")
    return sequence


def irt_standard(sequence, accessions=()):
    """Whether a peptide is a Biognosys iRT standard or a fragment of one.

    A peptide that PEAKS attributes only to the iRT kit (a homolog of a kit
    peptide, say) also counts, whatever its sequence.
    """
    return (any(sequence in peptide for peptide in IRT_PEPTIDES)
            or (bool(accessions) and all(value.startswith("Biognosys|") for value in accessions)))


def _sdrf(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    missing = set(SDRF_SAMPLE) | {"assay name", "comment[data file]", "characteristics[antibody enrichment]"}
    if missing - set(header) or any(len(row) != len(header) for row in rows[1:]):
        raise ValueError("Unexpected PXD074485 SDRF columns")
    named = [dict((name, row[i]) for i, name in enumerate(header) if name in missing) for row in rows[1:]]
    return {row["comment[data file]"]: row for row in named}


def _tables(path, expected):
    """Read every results table, checking the archive holds exactly the declared tables."""
    with zipfile.ZipFile(path) as archive:
        members = {name for name in archive.namelist() if not name.endswith("/")}
        if members != {f"Results/{name}.csv" for name in expected}:
            raise ValueError("PXD074485 results archive members changed")
        for name in sorted(expected):
            with archive.open(f"Results/{name}.csv") as raw:
                reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8-sig", newline=""))
                header = next(reader)
                yield name, header, list(enumerate(reader, 2))


def _columns(header):
    """Map each per-sample column to its role; PEAKS names them after the run label."""
    roles = []
    for column in header:
        role = column.split(" ", 1)[0] if column.startswith(("Area ", "Intensity ", "#Spec ", "#Feature ")) else column
        roles.append(role)
    if tuple(roles) != COLUMNS:
        raise ValueError("Unexpected PEAKS peptide table columns")
    labels = {column.split(" ", 1)[1] for column in header if " " in column and column.split(" ", 1)[0]
              in ("Area", "Intensity", "#Spec", "#Feature") and column not in ("#Spec", "#Feature")}
    if len(labels) != 1:
        raise ValueError("A PEAKS table must describe exactly one sample")
    label, = labels
    return label


def _value(text, line):
    if text == "":
        return None
    value = float(text)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid PEAKS quantity on line {line}")
    return value


def import_ligands(inputs, output, manifest):
    """Import every reported peptide row with its quantities, accessions and identification type."""
    if manifest.processing.get("importer") != "kaabinejadian-dla-v1":
        raise ValueError("Not a Kaabinejadian et al. DLA ligand manifest")
    assets = {asset.asset_id: asset for asset in manifest.sources}
    references = {reference.taxon_id: reference for reference in manifest.references}
    if (set(references) != {9615, 9606} or len(manifest.references) != 2
            or any(r.is_genomic or r.source_sha256 != assets["Results.zip"].sha256 for r in references.values())):
        raise ValueError("DLA ligands need dog and human source references pinned to the results archive")
    tables = manifest.processing["tables"]
    sdrf = _sdrf(inputs["sdrf.tsv"])
    samples, assays, sample_rows = {}, [], {}
    for name, table in sorted(tables.items()):
        row = sdrf.get(table["data_file"])
        if row is None or row["characteristics[antibody enrichment]"] != table["antibody"]:
            raise ValueError(f"{name} disagrees with the deposit's SDRF")
        sample_id = row["source name"]
        if sample_id != table["sample"]:
            raise ValueError(f"{name} names a different SDRF source")
        facts = {key: row[key] for key in SDRF_SAMPLE}
        if sample_rows.setdefault(sample_id, facts) != facts:
            raise ValueError(f"SDRF sample facts for {sample_id} disagree between assays")
        reference = references[9606 if facts["characteristics[organism]"] == "Homo sapiens" else 9615]
        if reference.taxon_id == 9615 and facts["characteristics[organism]"] != "Canis lupus familiaris":
            raise ValueError("Unexpected SDRF organism")
        typing = tuple(facts["characteristics[mhc typing]"].split(";"))
        provenance = dict(source_asset="sdrf.tsv", data_file=table["data_file"])
        if sample_id not in samples:
            dog = reference.taxon_id == 9615
            samples[sample_id] = Sample(
                sample_id, reference, biological_unit_id=facts["characteristics[individual]"],
                tissue=facts["characteristics[organism part]"], condition=facts["characteristics[disease]"],
                cell_type=None if dog else facts["characteristics[cell line]"],
                attributes=dict(material="tumour tissue" if dog else "DLA-88-transduced HCT 116 cell line",
                                sex=facts["characteristics[sex]"], age=facts["characteristics[age]"],
                                breed=BREEDS[sample_id] if dog else None,
                                breed_source="article Table 1" if dog else None,
                                mhc_typing=typing,
                                mhc_typing_source="deposit SDRF; article Table 1 agrees" if dog
                                else "deposit SDRF; the article's transduced allele",
                                source_attributes=facts),
                provenance=provenance)
        for database in table["databases"]:
            assays.append(Assay(
                f"{name}:{database}", sample_id, "mhc_ms", reference, run_id=row["assay name"],
                attributes=dict(antibody=table["antibody"], search_database=database,
                                search_database_description=manifest.coverage["search_databases"][database],
                                search_fasta="not deposited", identification="PEAKS Studio 11.5",
                                fdr="1% peptide FDR (PEAKS decoy-fusion), per the article",
                                raw_archive=table["data_file"],
                                instrument="SDRF: TripleTOF 5600; the article also reports a ZenoTOF 7600 "
                                           "and doesn't say which runs used it"),
                normalization=dict(method="none; PEAKS label-free area and intensity as exported"),
                calibration=dict(absolute_molecule_counts=False),
                provenance=dict(provenance, source_asset="Results.zip", table=f"Results/{name}-{database}.csv")))

    expected = {f"{name}-{database}" for name, table in tables.items() for database in table["databases"]}
    counts = manifest.coverage["source_counts"]
    parsed, accessions = [], {}
    for member, header, rows in _tables(inputs["Results.zip"], expected):
        name, database = member.rsplit("-", 1)
        table = tables[name]
        reference = samples[table["sample"]].reference
        label = _columns(header)
        found, irt = Counter(), 0
        for line, row in rows:
            if len(row) != len(header):
                raise ValueError(f"Short PEAKS row {member}:{line}")
            if row[25] not in FOUND_BY:
                raise ValueError(f"Unknown PEAKS identification type in {member}:{line}")
            found[row[25]] += 1
            listed = [value for value in row[24].split(":") if value]
            sequence = peptide_sequence(row[0])
            standard = irt_standard(sequence, listed)
            irt += standard
            accessions.setdefault((reference.taxon_id, database, row[25], sequence), set()).update(listed)
            run = row[13].rsplit("-Fxn", 1)[0]
            if run not in table["runs"]:
                raise ValueError(f"{member}:{line} cites an undeclared run {run!r}")
            parsed.append((member, f"{name}:{database}", table["sample"], reference, label, line, row,
                           sequence, standard, table["runs"][run]))
        if dict(found) != counts[member]["found_by"] or irt != counts[member]["irt_standard_rows"]:
            raise ValueError(f"PEAKS table {member} reconciliation failed: {dict(found)!r}, {irt} iRT rows")

    def measurements():
        for member, assay, sample_id, reference, label, line, row, sequence, standard, deposited in parsed:
            target = ScopedTarget(reference, "peptide", sequence, "amino-acid-sequence")
            provenance = dict(source_asset="Results.zip", table=member, line=line, peptide=row[0],
                              found_by=row[25], score_10lgp=float(row[1]), length=int(row[5]),
                              charge=int(row[8]), mz=row[7], ppm=row[6], rt=row[9], scan=row[12],
                              source_file=row[13], raw_file_deposited=deposited, ptm=row[20] or None,
                              features=int(row[19]), irt_standard=standard)
            spectra = int(row[17])
            if spectra < 1 or int(row[16]) != spectra:
                raise ValueError(f"Spectrum counts disagree in {member}:{line}")
            for column, quantity, value, unit in (
                    ("Area", "peptide_abundance", _value(row[14], line), "PEAKS area"),
                    ("Intensity", "peptide_abundance", _value(row[15], line), "PEAKS intensity"),
                    ("#Spec", "spectrum_count", spectra, "PSM")):
                yield Measurement(
                    f"{member}:{line}:{column}", target, sample_id, assay, quantity, value, unit,
                    manifest.dataset_id, manifest.version, status="missing" if value is None else "measured",
                    provenance=dict(provenance, column=f"{column} {label}"))

    def mappings():
        for (taxon, database, found_by, sequence), listed in sorted(accessions.items()):
            reference = references[taxon]
            standards = sorted(value for value in listed if value.startswith("Biognosys|"))
            proteins = sorted(value for value in listed if not value.startswith("Biognosys|"))
            accessions_only = sorted(value.split("|")[0] for value in proteins)
            if len(set(accessions_only)) != len(accessions_only):
                raise ValueError(f"Two {database} entries share an accession")
            # Accession containment is unverifiable: the search FASTAs were not deposited.
            evidence = dict(found_by=found_by, identification=FOUND_BY[found_by], irt_standards=standards)
            if irt_standard(sequence, sorted(listed)):
                reason = "Biognosys iRT retention-time standard"
            elif found_by == "DeepNovo" and proteins:
                reason = "De novo sequence; PEAKS's listed proteins need not contain it"
            elif not proteins:
                reason = "PEAKS reported no protein"
            else:
                reason = None
            if reason is not None:
                evidence["reported_proteins"] = accessions_only
                targets = ()
            else:
                targets = tuple(ScopedTarget(reference, "protein", value, f"PXD074485:{database}")
                                for value in accessions_only)
            status = "unique" if len(targets) == 1 else "ambiguous" if targets else "unmapped"
            yield TargetMapping(f"{database}:{found_by}:{sequence}", f"PXD074485:{database}:{found_by}",
                                sequence, reference, targets, status,
                                None if reason else f"PEAKS-reported proteins ({found_by})", reason, evidence,
                                dict(source_asset="Results.zip"))

    return QuantitativeDataset.create(output, manifest, samples=samples.values(), assays=assays,
                                      measurements=measurements(), mappings=mappings())
