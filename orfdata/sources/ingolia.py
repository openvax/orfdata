"""Native yeast RNA/Ribo windows; genomic identity and aliquot pairing stay explicit."""

import csv
import gzip
import math
import re
from collections import Counter

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample
from ..reference import ScopedTarget

COLUMNS = ("yorf", "norm", "dens", "count", "len", "total")
FORMULAS = {"dens": "count * 1000 / target_length_nt",
            "norm": "count * 1e9 / (target_length_nt * total_CDS_aligned_reads)"}
METADATA_FIELDS = (
    "geo_accession", "title", "source_name_ch1", "organism_ch1", "taxid_ch1",
    "characteristics_ch1", "treatment_protocol_ch1", "growth_protocol_ch1",
    "extract_protocol_ch1", "description", "data_processing", "relation", "instrument_model",
)


def ribo_manifest():
    """Return the pinned yeast rich/starvation source release."""
    from ..datasets import DatasetRegistry

    return DatasetRegistry().info("ingolia-yeast-rna-ribo")


def build_ribo(*, cache_root=None, source_paths=None, force=False, download_options=None):
    """Explicitly rebuild the reported RNA/Ribo measurements from pinned inputs."""
    from ..datasets import build_dataset

    manifest = ribo_manifest()
    return build_dataset(manifest, lambda inputs, output: import_ribo(inputs, output, manifest),
                         cache_root=cache_root, source_paths=source_paths, force=force,
                         download_options=download_options)


def _metadata(path, library):
    fields, accessions = {}, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            key, separator, value = line.rstrip("\r\n").partition(" = ")
            if not separator:
                continue
            if key == "^SAMPLE":
                accessions.append(value)
            elif key.startswith("!Sample_") and key[8:] in METADATA_FIELDS:
                fields.setdefault(key[8:], []).append(value)
    if (accessions != [library["sample_id"]]
            or fields.get("geo_accession") != accessions
            or fields.get("title") != [library["title"]]
            or fields.get("organism_ch1") != ["Saccharomyces cerevisiae"]
            or fields.get("taxid_ch1") != ["4932"]
            or not any("BY4741" in value for value in fields.get("characteristics_ch1", ()))):
        raise ValueError("Yeast source sample identity, title, organism or strain disagrees")
    for label, accession in (("BioSample", library["biosample"]), ("SRA", library["sra_experiment"])):
        candidates = [value for value in fields.get("relation", ()) if value.startswith(label + ": ")]
        if len(candidates) != 1 or re.split(r"[/=]", candidates[0])[-1] != accession:
            raise ValueError("Yeast source BioSample/SRA identity disagrees")
    return fields


def _feature_type(identifier):
    for suffix, kind in (("_3'utr", "three_prime_UTR"), ("_5'utr", "five_prime_UTR"), ("_intron", "intron")):
        if identifier.endswith(suffix):
            return kind
    return "unclassified_source_feature"


def import_ribo(inputs, output, manifest):
    """Preserve reported regions, values and normalization without invented ORFs."""
    if manifest.processing.get("importer") != "ingolia-yeast-ribo-v1" or len(manifest.references) != 1:
        raise ValueError("Not a single-reference Ingolia yeast manifest")
    reference, = manifest.references
    assets = {asset.asset_id: asset for asset in manifest.sources}
    if (reference.taxon_id != 4932 or reference.genome_build is not None
            or reference.annotation_name != "SGD" or reference.annotation_version != "2008-06-22"
            or reference.source_sha256 != assets["readMe.txt.gz"].sha256):
        raise ValueError("Yeast source reference must retain its unresolved historical genomic context")
    with gzip.open(inputs["readMe.txt.gz"], "rt", encoding="utf-8") as handle:
        readme = handle.read()
    if not all(text in readme for text in ("Saccharomyces Genome Database, 2008-06-22",
                                          "Nagalakshmi et al., Science 2008", "total CDS-aligned reads")):
        raise ValueError("Unexpected source region/reference documentation")
    libraries = manifest.processing["libraries"]
    if len(libraries) != manifest.coverage["assays"] or len({x["sample_id"] for x in libraries}) != len(libraries):
        raise ValueError("Duplicate or changed yeast assay inventory")
    samples, assays = [], []
    for library in libraries:
        sid = library["sample_id"]
        modality, condition, replicate = library["title"].split("-")
        if modality not in ("mRNA", "Footprints") or condition not in ("rich", "starved") or replicate not in ("1", "2"):
            raise ValueError("Unexpected source library title")
        kind = "rna_seq" if modality == "mRNA" else "ribo_seq"
        if library["assay_kind"] != kind:
            raise ValueError("Source library modality disagrees with manifest")
        fields = _metadata(inputs[library["metadata_asset"]], library)
        provenance = dict(source_asset=library["metadata_asset"], source_accession=sid)
        samples.append(Sample(sid, reference, attributes=dict(
            source_title=library["title"], source_condition=condition, source_replicate_label=replicate,
            biological_strain="BY4741", biosample=library["biosample"],
            sra_experiment=library["sra_experiment"], source_geo=fields,
            cross_assay_aliquot_correspondence="unresolved"), provenance=provenance))
        assays.append(Assay(sid, sid, kind, reference, search_reference=reference.key,
                            normalization=dict(FORMULAS, count="reported aligned read count",
                                               total_cds_reads=manifest.coverage["libraries"][sid]["total_cds_reads"],
                                               reported_decimal_places=3, source_windows="SGD2008-06-22; Nagalakshmi2008 UTRs; exact coordinates unresolved"),
                            calibration=dict(absolute_molecule_counts=False, absolute_translation_flux=False),
                            provenance=provenance))

    def observations():
        lengths, total_rows = {}, 0
        for library in libraries:
            sid, asset = library["sample_id"], library["quant_asset"]
            expected = manifest.coverage["libraries"][sid]
            seen, types, zeros = set(), Counter(), 0
            quantity = "rna_abundance" if library["assay_kind"] == "rna_seq" else "ribosome_occupancy"
            with gzip.open(inputs[asset], "rt", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if tuple(reader.fieldnames or ()) != COLUMNS:
                    raise ValueError("Unexpected yeast quantification columns")
                for line, row in enumerate(reader, 2):
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError("Incomplete yeast source row")
                    identifier = row["yorf"]
                    if not identifier or identifier in seen:
                        raise ValueError("Duplicate or missing yeast source feature")
                    seen.add(identifier)
                    if any(not re.fullmatch(r"[0-9]+", row[column]) for column in ("count", "len", "total")):
                        raise ValueError("Counts, target length and denominator must be nonnegative integers")
                    count, length, denominator = (int(row[column]) for column in ("count", "len", "total"))
                    if length < 1 or denominator != expected["total_cds_reads"] or denominator < 1:
                        raise ValueError("Invalid source length or normalization denominator")
                    if identifier in lengths and lengths[identifier] != length:
                        raise ValueError("Source feature lengths disagree between assays")
                    lengths[identifier] = length
                    reported = {column: float(row[column]) for column in ("dens", "norm")}
                    calculated = dict(dens=count * 1000 / length, norm=count * 1e9 / (length * denominator))
                    if any(not math.isfinite(value) or value < 0 or abs(value - calculated[column]) > 0.000501
                           for column, value in reported.items()):
                        raise ValueError("Reported source normalization disagrees with counts/length/denominator")
                    zeros += count == 0
                    kind = _feature_type(identifier)
                    types[kind] += 1
                    target = ScopedTarget(reference, "region", identifier, "Ingolia2009:feature")
                    count_id = f"{sid}:{line}:count"
                    provenance = dict(source_asset=asset, line=line, native_feature_type=kind,
                                      feature_length_nt=length, total_cds_reads=denominator)
                    for column, value, unit in (("count", count, "counts"), ("dens", reported["dens"], "reads/kb"),
                                                ("norm", reported["norm"], "RPKM")):
                        derived = column != "count"
                        yield Measurement(
                            f"{sid}:{line}:{column}", target, sid, sid, quantity, value, unit,
                            manifest.dataset_id, manifest.version, origin="derived" if derived else "observed",
                            input_measurements=(count_id,) if derived else (),
                            derivation=dict(method="source-reported normalized read density",
                                            formula=FORMULAS[column],
                                            reported_decimal_places=3, source_value_preserved=True) if derived else {},
                            provenance=dict(provenance, column=column, native_value=row[column]))
                    total_rows += 1
            if (len(seen) != manifest.coverage["features"] or seen != set(lengths)
                    or zeros != expected["zero_counts"] or dict(types) != dict(manifest.coverage["feature_types"])):
                raise ValueError("Source feature/zero/region inventory reconciliation failed")
        if total_rows != manifest.coverage["source_rows"] or total_rows * 3 != manifest.coverage["observations"]:
            raise ValueError("Source observation count reconciliation failed")

    return QuantitativeDataset.create(output, manifest, samples=samples, assays=assays, measurements=observations())
