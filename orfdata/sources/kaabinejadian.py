"""DLA class I ligands from dog tumours and DLA-88 monoallelic human cells (Kaabinejadian et al. 2026).

PRIDE PXD074485 deposits PEAKS Studio 11.5 peptide tables for eight
immunoprecipitations. Each dog tumour was searched against four protein
databases; each table is one search of the same spectra, so the same
peptide in two of a sample's tables is one observation, not two.

Every run was spiked with Biognosys iRT retention-time standards. Importer
``kaabinejadian-dla-v2`` flags a row as ``irt_standard`` by the first rule
that holds (:func:`irt_rules`):

``irt_sequence``
    its sequence, with I and L equal, is a kit peptide or a fragment of one;
``kit_only``
    PEAKS attributes it only to the kit;
``irt_precursor``
    it sits at the precursor of a row flagged by those two rules from the same
    immunoprecipitation: same charge, m/z and retention time within the
    manifest's tolerances. This applies to de novo and homology rows, and to
    database-search rows of a search without the kit (PEAKS couldn't choose
    the kit peptide there) whose sequence no search with the kit identified;
``de_novo_irt``
    it is a de novo sequence that PEAKS matches to the kit among other
    proteins, or that shares the manifest's minimum run of residues (I and L
    equal) with a kit peptide.

Importer ``kaabinejadian-dla-v1`` flagged only by the first two rules, with
exact sequences.
"""

import bisect
import csv
import io
import math
import re
import zipfile
from collections import Counter, defaultdict

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
    peptide, say) also counts, whatever its sequence. Sequences must match
    exactly, as importer ``kaabinejadian-dla-v1`` required; :func:`irt_rules`
    holds the current rules.
    """
    return (any(sequence in peptide for peptide in IRT_PEPTIDES)
            or (bool(accessions) and all(value.startswith("Biognosys|") for value in accessions)))


IRT_RULES = ("irt_sequence", "kit_only", "irt_precursor", "de_novo_irt")
IRT_REASONS = dict(irt_sequence="Biognosys iRT retention-time standard",
                   kit_only="Biognosys iRT retention-time standard",
                   irt_precursor="Observed at a Biognosys iRT standard's precursor: probably a misread iRT spectrum",
                   de_novo_irt="De novo read of a Biognosys iRT standard")


def _shared_run(sequence, peptides):
    """The longest run of residues a sequence shares with any of ``peptides``."""
    best = 0
    for size in range(1, len(sequence) + 1):
        if not any(sequence[i:i + size] in peptide for peptide in peptides
                   for i in range(len(sequence) - size + 1)):
            break
        best = size
    return best


def irt_rules(rows, kit_databases, precursor_ppm, precursor_rt_minutes, de_novo_min_shared_residues):
    """Each row's iRT rule (see the module docs), or None.

    ``rows`` are ``(table, database, found_by, sequence, accessions, charge,
    mz, rt)`` tuples, where ``table`` names one immunoprecipitation, whose
    searches share its spectra. ``kit_databases`` are the searches whose
    database included the kit.
    """
    kit = tuple(peptide.replace("I", "L") for peptide in IRT_PEPTIDES)
    rules = []
    for table, database, found_by, sequence, accessions, *_ in rows:
        if any(sequence.replace("I", "L") in peptide for peptide in kit):
            rules.append("irt_sequence")
        elif accessions and all(value.startswith("Biognosys|") for value in accessions):
            rules.append("kit_only")
        else:
            rules.append(None)
    anchors, kit_found = defaultdict(list), defaultdict(set)
    for rule, (table, database, found_by, sequence, _, charge, mz, rt) in zip(rules, rows):
        if rule:
            anchors[table, charge].append((mz, rt))
        elif found_by == "DB Search" and database in kit_databases:
            kit_found[table].add(sequence)
    for values in anchors.values():
        values.sort()

    def at_anchor(table, charge, mz, rt):
        candidates = anchors.get((table, charge), ())
        for anchor_mz, anchor_rt in candidates[bisect.bisect_left(candidates, (mz * (1 - precursor_ppm * 1e-6),)):]:
            if anchor_mz > mz * (1 + precursor_ppm * 1e-6):
                return False
            if abs(anchor_rt - rt) <= precursor_rt_minutes:
                return True
        return False

    for index, (table, database, found_by, sequence, accessions, charge, mz, rt) in enumerate(rows):
        if rules[index]:
            continue
        eligible = found_by != "DB Search" or (database not in kit_databases and sequence not in kit_found[table])
        if eligible and at_anchor(table, charge, mz, rt):
            rules[index] = "irt_precursor"
        elif found_by == "DeepNovo" and (
                any(value.startswith("Biognosys|") for value in accessions)
                or _shared_run(sequence.replace("I", "L"), kit) >= de_novo_min_shared_residues):
            rules[index] = "de_novo_irt"
    return rules


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
    importer = manifest.processing.get("importer")
    if importer not in ("kaabinejadian-dla-v1", "kaabinejadian-dla-v2"):
        raise ValueError("Not a Kaabinejadian et al. DLA ligand manifest")
    irt = None if importer == "kaabinejadian-dla-v1" else manifest.processing["parameters"]["irt"]
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
        for line, row in rows:
            if len(row) != len(header):
                raise ValueError(f"Short PEAKS row {member}:{line}")
            if row[25] not in FOUND_BY:
                raise ValueError(f"Unknown PEAKS identification type in {member}:{line}")
            listed = [value for value in row[24].split(":") if value]
            sequence = peptide_sequence(row[0])
            accessions.setdefault((reference.taxon_id, database, row[25], sequence), set()).update(listed)
            run = row[13].rsplit("-Fxn", 1)[0]
            if run not in table["runs"]:
                raise ValueError(f"{member}:{line} cites an undeclared run {run!r}")
            parsed.append((member, f"{name}:{database}", table["sample"], reference, label, line, row,
                           sequence, listed, table["runs"][run]))
    if irt is None:
        rules = ["irt_sequence" if irt_standard(sequence, listed) else None
                 for *_, sequence, listed, _ in parsed]
    else:
        rules = irt_rules([(member.rsplit("-", 1)[0], member.rsplit("-", 1)[1], row[25], sequence, listed,
                            int(row[8]), float(row[7]), float(row[9]))
                           for member, _, _, _, _, _, row, sequence, listed, _ in parsed],
                          irt["kit_databases"], irt["precursor_ppm"], irt["precursor_rt_minutes"],
                          irt["de_novo_min_shared_residues"])
    found, flagged, key_rules = defaultdict(Counter), defaultdict(Counter), {}
    for (member, _, _, reference, _, _, row, sequence, _, _), rule in zip(parsed, rules):
        found[member][row[25]] += 1
        if rule:
            flagged[member][rule] += 1
            key = (reference.taxon_id, member.rsplit("-", 1)[1], row[25], sequence)
            key_rules[key] = min(key_rules.get(key, rule), rule, key=IRT_RULES.index)
    for member in sorted(expected):
        expect = counts[member]
        if (dict(found[member]) != expect["found_by"]
                or sum(flagged[member].values()) != expect["irt_standard_rows"]
                or (irt is not None and dict(flagged[member]) != expect["irt_rules"])):
            raise ValueError(f"PEAKS table {member} reconciliation failed: {dict(found[member])!r}, "
                             f"{dict(flagged[member])!r} iRT rows")

    def measurements():
        for (member, assay, sample_id, reference, label, line, row, sequence, _, deposited), rule in zip(parsed, rules):
            target = ScopedTarget(reference, "peptide", sequence, "amino-acid-sequence")
            provenance = dict(source_asset="Results.zip", table=member, line=line, peptide=row[0],
                              found_by=row[25], score_10lgp=float(row[1]), length=int(row[5]),
                              charge=int(row[8]), mz=row[7], ppm=row[6], rt=row[9], scan=row[12],
                              source_file=row[13], raw_file_deposited=deposited, ptm=row[20] or None,
                              features=int(row[19]), irt_standard=rule is not None)
            if irt is not None:
                provenance["irt_rule"] = rule
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
            rule = key_rules.get((taxon, database, found_by, sequence))
            if irt is not None and rule:
                evidence["irt_rule"] = rule
            if irt is None and irt_standard(sequence, sorted(listed)):
                reason = "Biognosys iRT retention-time standard"
            elif irt is not None and rule:
                reason = IRT_REASONS[rule]
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
