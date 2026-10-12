"""Native TMT protein-group abundance from Oziolor et al. 2026 (three species)."""

import csv
import hashlib
import json
import math
import re
from collections import Counter

from ..dataset import QuantitativeDataset
from ..datasets import DatasetRegistry, build_dataset
from ..measurements import Assay, Measurement, Sample
from ..reference import ScopedTarget, species_taxon
from ..reference_mapping import ReferenceMapping

SAMPLE_COLUMNS = ("Run", "Channel", "BioReplicate", "Condition", "TechRepMixture", "Mixture", "sample.type")
VALUE_COLUMNS = ("IBAQ", "rIBAQ")
TABLE_COLUMNS = {*SAMPLE_COLUMNS, "ProteinGroup", "Protein", "Abundance", "sumIBAQ", *VALUE_COLUMNS}


def protein_sample_id(record):
    """Keep the protein-sample v1 identity already used by OpenVax consumers."""
    payload = {key: record[key] for key in SAMPLE_COLUMNS}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    return "openvax:protein-sample:v1:sha256:" + hashlib.sha256(encoded.encode("ascii")).hexdigest()


def protein_manifest(species):
    """Return the pinned source manifest, without downloading or opening data."""
    taxon = species_taxon(species)
    registry = DatasetRegistry()
    names = {record.dataset_id for record in registry.list(taxon_id=taxon)
             if taxon is not None and record.processing.get("importer") == "pfizer-protein-v1"}
    if len(names) != 1:
        raise ValueError("Pfizer protein summaries are registered for human, mouse and dog")
    return registry.info(names.pop())


def build_protein(species, *, cache_root=None, source_paths=None, force=False, download_options=None):
    """Explicitly build a verified source-native release using DataCache."""
    manifest = protein_manifest(species)
    return build_dataset(manifest, lambda inputs, output: import_protein(inputs, output, manifest),
                         cache_root=cache_root, source_paths=source_paths, force=force,
                         download_options=download_options)


def search_references(manifest):
    """Declare the exact deposited protein FASTA, independent of reported groups."""
    from ..search_reference import SearchReferenceSpec

    if manifest.processing.get("importer") != "pfizer-protein-v1":
        raise ValueError("Not a Pfizer protein source manifest")
    reference, = manifest.references
    asset = next(asset for asset in manifest.sources if asset.asset_id == "search.fasta")
    return (SearchReferenceSpec(manifest.dataset_id, manifest.version, reference, asset,
                                "uniprot", "uniprot-fasta", manifest.coverage["search_sequences"]),)


def _rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames) or set(reader.fieldnames) != TABLE_COLUMNS:
            raise ValueError("Unexpected Pfizer protein table columns")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None or not value for value in row.values()):
                raise ValueError(f"Incomplete protein source row {line}")
            if row["sample.type"] != "sample":
                raise ValueError("Summary must contain biological sample channels, not bridge measurements")
            yield line, row


def _value(text):
    if text == "NA":
        return None
    value = float(text)
    if not math.isfinite(value) or value < 0:
        raise ValueError("Native iBAQ must be finite and nonnegative, or source NA")
    return value


def _fasta(path, taxon):
    records = {}
    header, sequence = None, []

    def save():
        if header is None:
            return
        fields = header.split()[0].split("|")
        found = re.search(r" OX=(\d+)(?: |$)", header)
        if len(fields) != 3 or found is None or int(found.group(1)) != taxon:
            raise ValueError("Search FASTA identifier/taxon disagrees with the source manifest")
        identifier = fields[1]
        if identifier in records or not sequence:
            raise ValueError("Duplicate/empty search FASTA record")
        seq = "".join(sequence)
        if not re.fullmatch("[A-Z]+", seq):
            raise ValueError("Unexpected search FASTA sequence characters")
        gene = re.search(r" GN=(\S+)", header)
        records[identifier] = dict(header=header, gene_label=gene.group(1) if gene else None,
                                   sequence_sha256=hashlib.sha256(seq.encode("ascii")).hexdigest(),
                                   amino_acids=len(seq))

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError("Search FASTA sequence precedes its header")
            else:
                sequence.append(line)
    save()
    if not records:
        raise ValueError("Empty search FASTA")
    return records


def _audit(path):
    samples, groups, seen = {}, {}, set()
    counts = Counter()
    for _, row in _rows(path):
        sample_id = protein_sample_id(row)
        sample = {key: row[key] for key in SAMPLE_COLUMNS}
        samples[sample_id] = sample
        group, representative = row["ProteinGroup"], row["Protein"]
        if group in groups and groups[group] != representative:
            raise ValueError("Conflicting protein-group representatives")
        groups[group] = representative
        if (sample_id, group) in seen:
            raise ValueError("Duplicate protein-group/sample rows cannot be aggregated")
        seen.add((sample_id, group))
        counts["rows"] += 1
        for column in VALUE_COLUMNS:
            value = _value(row[column])
            counts[column + ":" + ("missing" if value is None else "zero" if value == 0 else "positive")] += 1
        # Keep native transformed fields as strings, including '-Inf'.
        if not math.isfinite(float(row["sumIBAQ"])) or float(row["sumIBAQ"]) <= 0:
            raise ValueError("Invalid source normalization sum")
    counts.update(groups=len(groups), assays=len(samples), absent_matrix_cells=len(groups) * len(samples) - counts["rows"])
    return samples, groups, dict(counts)


def import_protein(inputs, output, manifest):
    """Build from DataCache-verified inputs; preserve group ambiguity and units.

    Public callers should use build_protein or DatasetRegistry.build so input
    hashes are checked before this transformation runs. Memory holds metadata
    and source row identities, never a dense protein-by-sample matrix.
    """
    if manifest.processing.get("importer") != "pfizer-protein-v1" or len(manifest.references) != 1:
        raise ValueError("Not a single-species Pfizer protein manifest")
    reference, = manifest.references
    assets = {asset.asset_id: asset for asset in manifest.sources}
    if reference.sequence_sha256 != assets["search.fasta"].sha256:
        raise ValueError("Search reference and source FASTA fingerprints disagree")
    native_samples, groups, counts = _audit(inputs["protein.tsv"])
    expected = dict(manifest.coverage["source_counts"])
    if any(counts.get(key, 0) != value for key, value in expected.items()):
        raise ValueError(f"Source count reconciliation failed: {counts!r}")
    sequences = _fasta(inputs["search.fasta"], reference.taxon_id)
    if len(sequences) != manifest.coverage["search_sequences"]:
        raise ValueError("Search FASTA count reconciliation failed")
    accession = manifest.coverage["accession"]
    namespace = accession + ":ProteinGroup"
    provenance = dict(source_asset="protein.tsv", study=accession)

    def samples():
        for identifier, row in sorted(native_samples.items()):
            # Biological identity is the source pool label, not a run count.
            yield Sample(identifier, reference, biological_unit_id=accession + ":pool:" + row["BioReplicate"],
                         tissue=row["Condition"], condition="source tissue pool", pool_members=(),
                         attributes=dict(pooling="pooled", pool_membership="unresolved",
                                         donor_id=None, sex=None, biological_strain=None,
                                         health=None, developmental_stage=None, age=None,
                                         source_attributes=row,
                                         source_condition_includes_sex=reference.taxon_id == 10090),
                         provenance=provenance)

    def assays():
        for identifier, row in sorted(native_samples.items()):
            yield Assay(identifier, identifier, "digest_ms", reference, run_id=row["Run"],
                        technical_replicate=row["TechRepMixture"], channel=row["Channel"], plex=row["Mixture"],
                        bridge="131C", search_reference=reference.key,
                        normalization=dict(method="source MSstatsTMT bridge normalization and weighted iBAQ",
                                           IBAQ="source column; native linear iBAQ",
                                           rIBAQ="source column; relative iBAQ; retain original scale",
                                           between_species_calibration=False),
                        calibration=dict(method="none to copies/cell", absolute_molecule_counts=False),
                        attributes=dict(fractions="24 per TMT batch; source table summarizes fractions",
                                        instrument="unresolved: paper Eclipse versus deposit Fusion Lumos",
                                        source_attributes=row), provenance=provenance)

    def mappings():
        for group, representative in sorted(groups.items()):
            found = sequences.get(representative)
            yield ReferenceMapping(
                "representative:" + group, manifest.dataset_id, manifest.version,
                ScopedTarget(reference, "protein_group", group, namespace), reference,
                (ScopedTarget(reference, "protein", representative, "uniprot"),) if found else (),
                "unresolved" if found else "unmapped", "source-reported representative lookup",
                dict(group_membership="not resolved by summary table", sequence_assignment="not tested",
                     representative_search_record=found),
                reason="Representative does not resolve all group members" if found else "Representative absent from deposited FASTA",
                provenance=dict(source_asset="protein.tsv", search_asset="search.fasta", representative=representative))

    targets = {group: ScopedTarget(reference, "protein_group", group, namespace) for group in groups}

    def measurements():
        for line, row in _rows(inputs["protein.tsv"]):
            identifier = protein_sample_id(row)
            for column in VALUE_COLUMNS:
                value = _value(row[column])
                yield Measurement(
                    f"row:{line}:{column}", targets[row["ProteinGroup"]], identifier, identifier,
                    "protein_abundance", value, "iBAQ" if column == "IBAQ" else "rIBAQ",
                    manifest.dataset_id, manifest.version, status="missing" if value is None else "measured",
                    provenance=dict(source_asset="protein.tsv", line=line, column=column,
                                    native_value=row[column], representative=row["Protein"],
                                    native_Abundance=row["Abundance"], native_sumIBAQ=row["sumIBAQ"]))

    return QuantitativeDataset.create(output, manifest, samples=samples(), assays=assays(),
                                      measurements=measurements(), reference_mappings=mappings())
