"""Native atlas RNA counts, without ortholog projection or inferred protein pairs."""

import csv
import gzip
import re
from collections import Counter

from ..dataset import QuantitativeDataset
from ..measurements import Assay, Measurement, Sample
from ..reference import ScopedTarget, species_taxon

COUNT_COLUMNS = ("sid", "species", "tnm", "egid", "rcnt")
SAMPLE_COLUMNS = ("sid", "species", "strain", "batch", "tnm", "grp", "type", "gene.count",
                  "mapped.reads", "mapped.percent", "run", "bases", "biosample")
SOFT_FIELDS = ("geo_accession", "title", "organism_ch1", "taxid_ch1", "characteristics_ch1",
               "description", "relation", "data_processing")


def rna_manifest(species):
    """Return the pinned native-count release; human source coverage is absent."""
    from ..datasets import DatasetRegistry

    names = {9615: "dog", 10090: "mouse"}
    taxon = species_taxon(species)
    if taxon not in names:
        raise ValueError("Native Pfizer RNA releases currently cover dog and mouse; human counts are not deposited")
    return DatasetRegistry().info(f"pfizer-atlas-{names[taxon]}-rna")


def build_rna(species, *, cache_root=None, source_paths=None, force=False, download_options=None):
    """Explicitly rebuild native counts from verified source files."""
    from ..datasets import build_dataset

    manifest = rna_manifest(species)
    return build_dataset(manifest, lambda inputs, output: import_rna(inputs, output, manifest),
                         cache_root=cache_root, source_paths=source_paths, force=force,
                         download_options=download_options)


def _geo_samples(inputs, assets):
    samples = {}
    for asset in assets:
        current = None
        with gzip.open(inputs[asset], "rt", encoding="utf-8") as handle:
            for line, text in enumerate(handle, 1):
                key, separator, value = text.rstrip("\r\n").partition(" = ")
                if not separator:
                    continue
                if key == "^SAMPLE":
                    if value in samples:
                        raise ValueError("Duplicate GEO sample accession across inputs")
                    current = dict(accession=value, source_asset=asset, source_line=line, fields={})
                    samples[value] = current
                elif key.startswith("^"):
                    current = None
                elif current is not None and key.startswith("!Sample_") and key[8:] in SOFT_FIELDS:
                    current["fields"].setdefault(key[8:], []).append(value)
    by_sid = {}
    for accession, sample in samples.items():
        fields = sample["fields"]
        if fields.get("geo_accession") != [accession] or len(fields.get("title", ())) != 1:
            raise ValueError("GEO sample accession/title is missing or inconsistent")
        sid, separator, _ = fields["title"][0].partition(",")
        if not separator or not re.fullmatch(r"[0-9]+", sid) or sid in by_sid:
            raise ValueError("Duplicate/invalid native sample ID in GEO metadata")
        by_sid[sid] = sample
    return by_sid


def _sample_table(path):
    from openpyxl import load_workbook
    book = load_workbook(path, read_only=True, data_only=True)
    records = {}
    try:
        rows = book["TableS5"].values
        if next(rows) != SAMPLE_COLUMNS:
            raise ValueError("Unexpected atlas sample-table columns")
        for line, row in enumerate(rows, 2):
            record = dict(zip(SAMPLE_COLUMNS, row))
            if type(record["sid"]) is not int or record["sid"] < 1 or record["type"] != "tx":
                raise ValueError("Invalid native RNA sample-table identity")
            sid = str(record["sid"])
            if sid in records:
                raise ValueError("Duplicate sample-table ID")
            records[sid] = dict(record, source_line=line)
    finally:
        book.close()
    return records


def _genes(path, reference, expected):
    genes, header = set(), {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#!"):
                key, _, value = line[2:].rstrip().partition(" ")
                header[key] = value
            elif not line.startswith("#"):
                columns = line.rstrip().split("\t")
                if len(columns) != 9:
                    raise ValueError("Invalid source GTF record")
                if columns[2] == "gene":
                    attributes = dict(re.findall(r'(\w+) "([^\"]*)";', columns[8]))
                    gene = attributes.get("gene_id")
                    if not gene or gene in genes:
                        raise ValueError("Missing/duplicate GTF gene ID")
                    genes.add(gene)
    if (header.get("genome-build") != reference.genome_build
            or header.get("genome-build-accession") != reference.assembly_accession):
        raise ValueError("GTF assembly disagrees with native RNA reference")
    if len(genes) != expected:
        raise ValueError("GTF gene-count reconciliation failed")
    return genes


def import_rna(inputs, output, manifest):
    """Stream one species from pinned native counts; retain original observation scope."""
    if manifest.processing.get("importer") != "pfizer-rna-v1" or len(manifest.references) != 1:
        raise ValueError("Not a single-species native Pfizer RNA manifest")
    reference, = manifest.references
    processing, coverage = manifest.processing, manifest.coverage
    species = processing["native_species"]
    assets = {source.asset_id: source for source in manifest.sources}
    annotation = assets["annotation.gtf.gz"]
    if (reference.taxon_id != species_taxon(species) or reference.annotation_name != "ensembl"
            or reference.annotation_version != "110" or reference.source_sha256 != annotation.sha256
            or annotation.metadata.get("ensembl_release") != "110"):
        raise ValueError("Native RNA species/reference and annotation pins disagree")
    genes = _genes(inputs["annotation.gtf.gz"], reference, coverage["reference_genes"])
    geo = _geo_samples(inputs, processing["geo_metadata_assets"])
    table = _sample_table(inputs["samples.xlsx"])
    if len(geo) != coverage["metadata_samples"] or set(geo) != set(table):
        raise ValueError("GEO and sample-table identities do not reconcile")
    selected = {sid: row for sid, row in table.items() if row["species"].split("_")[0] == species}
    if len(selected) != coverage["samples"]:
        raise ValueError("Native species sample count changed")
    samples, assays = [], []
    for sid, row in selected.items():
        sample = geo[sid]
        fields = sample["fields"]
        if fields.get("taxid_ch1") != [str(reference.taxon_id)]:
            raise ValueError("GEO sample taxon disagrees with native count species")
        organism, = fields["organism_ch1"]
        if species_taxon(organism) != reference.taxon_id:
            raise ValueError("GEO sample organism disagrees with native count species")
        biosamples = [value.rsplit("/", 1)[-1] for value in fields.get("relation", ())
                      if value.startswith("BioSample: ")]
        if biosamples != [row["biosample"]]:
            raise ValueError("GEO and sample table BioSample identities disagree")
        provenance = dict(source_asset="samples.xlsx", sheet="TableS5", line=row["source_line"],
                          geo_asset=sample["source_asset"], geo_line=sample["source_line"])
        samples.append(Sample(sid, reference, biological_unit_id=row["biosample"], tissue=row["tnm"],
                              attributes=dict(native_sid=sid, gsm=sample["accession"], biosample=row["biosample"],
                                              biological_strain=row["strain"], sex=None, donor_id=None,
                                              source_sample_table=row, source_geo=fields,
                                              cross_assay_sample_correspondence="unresolved"), provenance=provenance))
        assays.append(Assay(sid, sid, "rna_seq", reference, run_id=row["run"], search_reference=reference.key,
                            normalization=dict(method="source integer gene count (rcnt)",
                                               pipeline="modified nf-core RNAseq 3.14.0; Nextflow 23.10.0; Salmon 1.10.1",
                                               gene_aggregation_and_rounding="not specified in deposited count table",
                                               ortholog_projection=False, length_normalized=False),
                            calibration=dict(absolute_molecule_counts=False), provenance=provenance))

    def observations():
        seen_samples, seen_genes, all_genes, current = set(), set(), set(), None
        source_rows = Counter()
        measured, zeros = 0, 0
        with gzip.open(inputs["counts.csv.gz"], "rt", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != COUNT_COLUMNS:
                raise ValueError("Unexpected native RNA count columns")
            for line, row in enumerate(reader, 2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError("Incomplete native RNA row")
                source_rows[row["species"]] += 1
                if row["species"] != species:
                    continue
                sid, gene = row["sid"], row["egid"]
                if sid not in selected or row["tnm"] != selected[sid]["tnm"]:
                    raise ValueError("Native RNA sample/tissue disagrees with sample metadata")
                if sid != current:
                    if sid in seen_samples:
                        raise ValueError("Native RNA sample blocks must not repeat")
                    if current is not None and len(seen_genes) != coverage["genes"]:
                        raise ValueError("Native RNA per-sample gene count changed")
                    seen_samples.add(sid)
                    current, seen_genes = sid, set()
                if gene not in genes or gene in seen_genes:
                    raise ValueError("Unmapped/duplicate gene in native RNA sample")
                seen_genes.add(gene)
                all_genes.add(gene)
                if not re.fullmatch(r"[0-9]+", row["rcnt"]):
                    raise ValueError("Native RNA counts must be nonnegative integers")
                value = int(row["rcnt"])
                measured += 1
                zeros += value == 0
                yield Measurement(f"row:{line}", ScopedTarget(reference, "gene", gene, "ensembl"),
                                  sid, sid, "rna_abundance", value, "counts", manifest.dataset_id, manifest.version,
                                  provenance=dict(source_asset="counts.csv.gz", line=line, native_value=row["rcnt"]))
        if (dict(source_rows) != dict(coverage["source_species_rows"]) or seen_samples != set(selected)
                or len(all_genes) != coverage["genes"] or len(seen_genes) != coverage["genes"]
                or measured != coverage["rows"] or zeros != coverage["zeros"]):
            raise ValueError("Native RNA source count reconciliation failed")

    return QuantitativeDataset.create(output, manifest, samples=samples, assays=assays, measurements=observations())
