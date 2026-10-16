"""Peptides and proteins from the Canine PeptideAtlas (Kusebauch et al. 2026), build 573.

PeptideAtlas reprocessed 138 dog runs (Labrador tissues, plasma and urine)
with the Trans-Proteomic Pipeline at a PSM FDR of 0.001 per experiment. Its
TSV export names every peptide, its PSM count in each sample, every protein
it maps to (with ROS_Cfam_1.0 coordinates for Ensembl proteins) and the
build's protein calls. PeptideAtlas states no licence for the export, so
orfdata builds this dataset from the original file and never redistributes it.
"""

import csv
import gzip
import hashlib
import io
import re
import zipfile
from collections import Counter, defaultdict

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample, TargetMapping
from ..reference import ScopedTarget
from ..reference_mapping import ReferenceMapping

TABLES = ("atlas_build", "atlas_search_batch", "sample", "sample_category", "instrument", "proteases",
          "peptide", "peptide_instance", "peptide_instance_search_batch", "peptide_mapping", "biosequence",
          "protein_identification", "protein_presence_level")
NULL = "\\N"
LOCATION = ("matched_biosequence_id", "start_in_biosequence", "end_in_biosequence", "chromosome", "strand",
            "start_in_chromosome", "end_in_chromosome")


def build_manifest():
    """Return the pinned build-573 release, without downloading anything."""
    from ..datasets import DatasetRegistry

    return DatasetRegistry().info("peptideatlas-dog-peptides")


def build_peptides(*, cache_root=None, source_paths=None, force=False, download_options=None):
    """Download the pinned PeptideAtlas export and Ensembl proteins, then build locally."""
    from ..datasets import build_dataset

    manifest = build_manifest()
    return build_dataset(manifest, lambda inputs, output: import_build(inputs, output, manifest),
                         cache_root=cache_root, source_paths=source_paths, force=force,
                         download_options=download_options)


def _namespace(name):
    if name.startswith("DECOY_"):
        return "PeptideAtlas:decoy"
    if name.startswith("CONTAM_"):
        return "PeptideAtlas:contaminant"
    if re.fullmatch(r"ENSCAFP[0-9]+\.[0-9]+", name):
        return "ensembl"
    if re.fullmatch(r"[A-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{3}[0-9]|A0A[A-Z0-9]{7}", name):
        return "uniprot"
    raise ValueError(f"Unexpected PeptideAtlas protein {name!r}")


def _rows(archive, table):
    with archive.open(f"{table}.tsv") as raw:
        reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8", newline=""), delimiter="\t",
                            quoting=csv.QUOTE_NONE)
        header = next(reader)
        if header[-1] == "":  # The export ends each line with a tab.
            header = header[:-1]
        for row in reader:
            if row and row[-1] == "" and len(row) == len(header) + 1:
                row = row[:-1]
            if len(row) != len(header):
                raise ValueError(f"Malformed PeptideAtlas {table} row")
            yield dict(zip(header, row))


def _fasta(path):
    found, name, sequence = {}, None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    found[name] = "".join(sequence)
                name, sequence = line[1:].split()[0], []
            else:
                sequence.append(line.strip())
    if name is not None:
        found[name] = "".join(sequence)
    return found


def _sha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _int(text):
    return None if text == NULL else int(text)


def _float(text):
    return None if text == NULL else float(text)


def import_build(inputs, output, manifest):
    """Import per-sample peptide PSM counts, peptide-to-protein locations and protein calls."""
    if manifest.processing.get("importer") != "peptideatlas-build-v1":
        raise ValueError("Not a PeptideAtlas build manifest")
    assets = {asset.asset_id: asset for asset in manifest.sources}
    source, ensembl = manifest.references
    if (source.source_sha256 != assets["atlas_build.tsv.zip"].sha256 or source.annotation_name == "ensembl"
            or (ensembl.annotation_name, ensembl.annotation_version) != ("ensembl", "110")
            or {source.taxon_id, ensembl.taxon_id} != {9615}):
        raise ValueError("PeptideAtlas needs its export-pinned source reference, then dog Ensembl 110")
    expected = manifest.coverage["source_counts"]
    build_id = manifest.processing["atlas_build_id"]
    with zipfile.ZipFile(inputs["atlas_build.tsv.zip"]) as archive:
        if set(TABLES) - {name[:-4] for name in archive.namelist()}:
            raise ValueError("PeptideAtlas export lacks a required table")
        build, = (row for row in _rows(archive, "atlas_build") if row["atlas_build_id"] == build_id)
        batches = {row["atlas_search_batch_id"]: row for row in _rows(archive, "atlas_search_batch")}
        native_samples = {row["sample_id"]: row for row in _rows(archive, "sample")}
        categories = {row["id"]: row["name"] for row in _rows(archive, "sample_category")}
        instruments = {row["instrument_id"]: row["instrument_name"] for row in _rows(archive, "instrument")}
        proteases = {row["id"]: row["name"] for row in _rows(archive, "proteases")}
        levels = {row["protein_presence_level_id"]: row["level_name"]
                  for row in _rows(archive, "protein_presence_level")}
        sequences = {row["peptide_id"]: (row["peptide_accession"], row["peptide_sequence"])
                     for row in _rows(archive, "peptide")}
        instances = {row["peptide_instance_id"]: row for row in _rows(archive, "peptide_instance")
                     if row["atlas_build_id"] == build_id}
        proteins = {row["biosequence_id"]: (row["biosequence_name"], row["biosequence_seq"])
                    for row in _rows(archive, "biosequence")}
        locations = defaultdict(list)
        for row in _rows(archive, "peptide_mapping"):
            locations[row["peptide_instance_id"]].append(tuple(row[key] for key in LOCATION))
        calls = [row for row in _rows(archive, "protein_identification") if row["atlas_build_id"] == build_id]
        counts = [(row["peptide_instance_id"], row["atlas_search_batch_id"], int(row["n_observations"]))
                  for row in _rows(archive, "peptide_instance_search_batch")]

    # The export stores the FDR as a 32-bit float (0.00100000004749745).
    if abs(float(build["protpro_PSM_FDR_per_expt"]) - manifest.coverage["psm_fdr_per_experiment"]) > 1e-9:
        raise ValueError("PeptideAtlas build FDR changed")
    if set(batches) != set(manifest.processing["search_batches"]):
        raise ValueError("PeptideAtlas search batches changed")
    observed = Counter(dict(peptides=len(instances), protein_calls=len(calls), search_batch_rows=len(counts),
                            peptide_mappings=sum(len(v) for v in locations.values()),
                            psms=sum(count for _, _, count in counts),
                            searched_spectra=sum(int(row["n_searched_spectra"]) for row in batches.values()),
                            runs=sum(int(row["n_runs"]) for row in batches.values())))
    observed.update(levels[row["presence_level_id"]] for row in calls)
    if dict(observed) != dict(expected):
        raise ValueError(f"PeptideAtlas reconciliation failed: {dict(observed)!r}")

    peptide_target = {}
    for identifier, row in instances.items():
        accession, sequence = sequences[row["peptide_id"]]
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
            raise ValueError(f"Unexpected PeptideAtlas peptide {accession}")
        peptide_target[identifier] = (accession, ScopedTarget(source, "peptide", sequence, "amino-acid-sequence"))

    def protein(biosequence_id):
        name, _ = proteins[biosequence_id]
        return ScopedTarget(source, "protein", name, _namespace(name))

    provenance = dict(source_asset="atlas_build.tsv.zip", atlas_build_id=build_id)
    sample_of = {}
    samples, assays = [], []
    for batch_id, batch in sorted(batches.items()):
        native = native_samples[batch["sample_id"]]
        sample_id = native["sample_accession"]
        sample_of[batch_id] = sample_id
        # Tags name the tissue, then any preparation: "LargeIntestine", "Plasma_depleted_timsTOF".
        tissue, *preparation = native["sample_tag"].split("_")
        tissue = re.sub(r"(?<!^)(?=[A-Z])", " ", tissue).lower()
        fluid = tissue in ("plasma", "urine")
        samples.append(Sample(
            sample_id, source, tissue=tissue,
            attributes=dict(sample_title=native["sample_title"], preparation=" ".join(preparation) or None,
                            sample_category=categories.get(native["sample_category_id"]),
                            sample_tag=native["sample_tag"], peptide_source_type=native["peptide_source_type"],
                            breed=manifest.coverage["breed"],
                            breed_source="the article" if fluid else "the deposit's SDRF and the article",
                            repository_identifiers=native["repository_identifiers"]),
            provenance=dict(provenance, table="sample.tsv", sample_id=batch["sample_id"])))
        assays.append(Assay(
            f"search-batch-{batch_id}", sample_id, "digest_ms", source, run_id=None,
            attributes=dict(runs=int(batch["n_runs"]), searched_spectra=int(batch["n_searched_spectra"]),
                            instrument=instruments.get(native["instrument_model_id"]),
                            protease=proteases.get(native["protease_id"]),
                            data_location=batch["data_location"],
                            psm_fdr_per_experiment=float(build["protpro_PSM_FDR_per_expt"])),
            normalization=dict(method="none; PSM counts"), calibration=dict(absolute_molecule_counts=False),
            provenance=dict(provenance, table="atlas_search_batch.tsv", atlas_search_batch_id=batch_id)))
    build_sample = f"build-{build_id}"
    samples.append(Sample(build_sample, source, condition="all samples of the build, combined by PeptideAtlas",
                          pool_members=tuple(sorted(sample_of.values())),
                          attributes=dict(pooling="combined analysis, not a physical pool"),
                          provenance=dict(provenance, table="atlas_build.tsv")))
    assays.append(Assay(build_sample, build_sample, "digest_ms", source,
                        attributes=dict(protein_inference="ProteinProphet presence levels"),
                        calibration=dict(absolute_molecule_counts=False),
                        provenance=dict(provenance, table="protein_identification.tsv")))

    def measurements():
        for instance, batch, count in counts:
            accession, target = peptide_target[instance]
            yield Measurement(
                f"search-batch-{batch}:{accession}", target, sample_of[batch], f"search-batch-{batch}",
                "spectrum_count", count, "PSM", manifest.dataset_id, manifest.version,
                provenance=dict(source_asset="atlas_build.tsv.zip", peptide_accession=accession))
        names = {key: value[0] for key, value in proteins.items()}
        for row in calls:
            yield Measurement(
                f"protein:{row['protein_identification_id']}", protein(row["biosequence_id"]),
                build_sample, build_sample, "spectrum_count", int(row["n_observations"]), "PSM",
                manifest.dataset_id, manifest.version,
                provenance=dict(provenance, table="protein_identification.tsv",
                                presence_level=levels[row["presence_level_id"]],
                                distinct_peptides=_int(row["n_distinct_peptides"]),
                                probability=_float(row["probability"]),
                                protein_group=_int(row["protein_group_number"]),
                                represented_by=names.get(row["represented_by_biosequence_id"]),
                                subsumed_by=names.get(row["subsumed_by_biosequence_id"]),
                                counting="PSMs of every observed peptide that maps to the protein, shared or not"))

    def mappings():
        for identifier, row in sorted(instances.items(), key=lambda item: peptide_target[item[0]][0]):
            accession, peptide = peptide_target[identifier]
            sequence = peptide.identifier
            targets, places, decoys = {}, [], 0
            for biosequence, start, end, chromosome, strand, first, last in sorted(
                    locations[identifier], key=lambda place: (proteins[place[0]][0], int(place[1]))):
                name, protein_sequence = proteins[biosequence]
                start, end = int(start), int(end)
                found = protein_sequence[start - 1:end]
                # PeptideAtlas maps isoleucine and leucine as equal, since MS can't tell them apart.
                if found.replace("I", "L") != sequence.replace("I", "L"):
                    raise ValueError(f"{accession} isn't at its stated position in {name}")
                if name.startswith("DECOY_"):
                    decoys += 1
                    continue
                targets.setdefault(name, protein(biosequence))
                places.append([name, start, end, None if chromosome == NULL else chromosome,
                               None if strand in ("?", NULL) else strand, _int(first), _int(last),
                               found == sequence])
            found = tuple(targets.values())
            status = "unique" if len(found) == 1 else "ambiguous" if found else "unmapped"
            yield TargetMapping(
                accession, "amino-acid-sequence", sequence, source, found, status,
                None if status == "unmapped" else "PeptideAtlas peptide mapping, checked against each sequence",
                "Maps only to decoy sequences" if status == "unmapped" else None,
                evidence=dict(locations=places, decoy_locations=decoys, best_probability=_float(row["best_probability"]),
                              psms=int(row["n_observations"]), samples=int(row["n_samples"]),
                              genome_locations=_int(row["n_genome_locations"]),
                              exon_spanning=row["is_exon_spanning"] == "Y",
                              empirical_proteotypic_score=_float(row["empirical_proteotypic_score"])),
                provenance=dict(source_asset="atlas_build.tsv.zip"))

    ensembl_proteins = _fasta(inputs["ensembl.pep.fa.gz"])
    by_sequence = defaultdict(list)
    for name, sequence in ensembl_proteins.items():
        by_sequence[sequence].append(name)

    def reference_mappings():
        called = {row["biosequence_id"] for row in calls}
        for biosequence_id in sorted(called, key=lambda key: proteins[key][0]):
            name, sequence = proteins[biosequence_id]
            namespace = _namespace(name)
            if namespace not in ("uniprot", "ensembl"):
                continue
            matches = sorted(by_sequence.get(sequence, ()))
            if namespace == "ensembl":
                matches = [match for match in matches if match == name]
            targets = tuple(ScopedTarget(ensembl, "protein", match.split(".")[0], "ensembl") for match in matches)
            status = "unique" if len(targets) == 1 else "ambiguous" if targets else "unmapped"
            yield ReferenceMapping(
                f"protein:{name}", manifest.dataset_id, manifest.version, protein(biosequence_id), ensembl,
                targets, status, "identical full-length protein sequence",
                dict(sequence=dict(source_sha256=_sha(sequence), target_sha256=_sha(sequence),
                                   all_candidates_checked=True) if targets else "no identical Ensembl 110 protein",
                     versioned_ensembl_ids=matches),
                reason=None if targets else "No Ensembl 110 protein has this sequence",
                provenance=dict(provenance, ensembl_asset="ensembl.pep.fa.gz"))

    return QuantitativeDataset.create(output, manifest, samples=samples, assays=assays,
                                      measurements=measurements(), mappings=mappings(),
                                      reference_mappings=reference_mappings())
