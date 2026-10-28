"""Healthy human tissue measurements with source target and calibration scope."""

import gzip
import math
import re
import zipfile

from ..reference import ScopedTarget
from ..search_reference import SearchReference, SearchReferenceSpec

SERIES = {
    "protein_groups": ("protein.zip", "Table_EV1.xlsx", "B. Protein groups", 4),
    "protein_genes": ("protein.zip", "Table_EV1.xlsx", "C. Genes", 2),
    "rna_transcripts": ("rna.zip", "Table_EV2.xlsx", "A. Transcripts", 4),
    "rna_genes": ("rna.zip", "Table_EV2.xlsx", "B. Genes", 2),
    "protein_copies": ("copies.zip", "Table_EV5.xlsx", "A. Protein copies", 3),
    "rna_copies": ("copies.zip", "Table_EV5.xlsx", "B. Transcript copies", 3),
}
ID_HEADERS = {
    "protein_groups": ("Protein IDs", "Protein ID (representatives)", "Gene ID (Ensembl)", "Gene name"),
    "protein_genes": ("Gene ID", "Gene name"),
    "rna_transcripts": ("Transcript ID", "Protein ID", "Gene ID", "Gene name"),
    "rna_genes": ("Gene ID", "Gene name"),
    "protein_copies": ("Protein ID", "Gene ID", "Gene name"),
    "rna_copies": ("Transcript ID", "Gene ID", "Gene name"),
}
SAMPLE_HEADERS = ("Tissue name", "Synonyms", "Tissue ID", "Tissue status", "RNA-Seq Data ID",
                  "Gender of donor", "Age of donor", "MS Experiment ID", "Proteases",
                  "Fragmentation techniques", "Mass spectrometer", "Search engine", "Used for…")


def search_references(manifest):
    """Declare the complete native search set, including alternate-reference loci."""
    if manifest.processing.get("importer") != "wang-human-tissue-v1":
        raise ValueError("Not a Wang healthy-tissue source manifest")
    reference, = manifest.references
    assets = {a.asset_id: a for a in manifest.sources}
    return (SearchReferenceSpec(manifest.dataset_id, manifest.version, reference,
                                assets["search.fasta.gz"], "ensembl-protein-versioned",
                                "ensembl-pep-fasta-gzip", manifest.coverage["search_sequences"]),)


def _sheet(path, member, sheet):
    from openpyxl import load_workbook
    with zipfile.ZipFile(path) as archive, archive.open(member) as stream:
        book = load_workbook(stream, read_only=True, data_only=True)
        try:
            yield from book[sheet].values
        finally:
            book.close()


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError("Abundance cell must be an explicit number")
    if isinstance(value, str) and not re.fullmatch(r"[+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", value):
        raise ValueError("Unknown native abundance sentinel")
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError("Native abundance must be finite and nonnegative")
    if number == 0 and isinstance(value, str) and any(c in "123456789" for c in value.split("e")[0].split("E")[0]):
        raise ValueError("Native numeric text underflows; do not replace it with zero")
    return number


def _rows(inputs, manifest, series):
    asset, member, sheet, count = SERIES[series]
    rows = _sheet(inputs[asset], member, sheet)
    header = next(rows)
    n = len(manifest.processing["tissues"])
    if header[:count] != ID_HEADERS[series]:
        raise ValueError("Unexpected native abundance identifier columns")
    tissue_columns = header[count:count + n]
    aliases = manifest.processing["tissue_aliases"]
    tissues = tuple(aliases.get(tissue, tissue) for tissue in tissue_columns)
    if len(set(tissues)) != n or set(tissues) != set(manifest.processing["tissues"]):
        raise ValueError("Native tissue columns do not reconcile")
    seen, found, empty = set(), 0, 0
    for line, row in enumerate(rows, 2):
        if all(value is None for value in row):
            empty += 1
            continue
        if row[0] is None or row[0] in seen:
            raise ValueError("Missing/duplicate native abundance target")
        seen.add(row[0])
        found += 1
        yield line, row[:count], tuple(zip(tissues, tissue_columns, row[count:count + n]))
    coverage = manifest.coverage["series"][series]
    if found != coverage["rows"] or empty != coverage["empty_rows"]:
        raise ValueError("Native abundance table row coverage changed")


def _samples(inputs, manifest):
    rows = _sheet(inputs["protein.zip"], "Table_EV1.xlsx", "A. Sample information")
    if next(rows) != SAMPLE_HEADERS:
        raise ValueError("Unexpected healthy-tissue sample columns")
    found, matched = {}, {}
    for line, row in enumerate(rows, 2):
        if "Quantitative proteome" not in str(row[-1]):
            continue
        tissue = row[0]
        if tissue in found:
            raise ValueError("Duplicate quantitative tissue specimen")
        found[tissue] = row
        if tissue in manifest.processing["tissues"]:
            if row[4] in (None, "-"):
                raise ValueError("Matched tissue lacks an RNA sample identifier")
            matched[tissue] = dict(zip(SAMPLE_HEADERS, row), source_line=line)
    if (len(found) != manifest.coverage["quantitative_metadata_rows"]
            or set(matched) != set(manifest.processing["tissues"])):
        raise ValueError("Quantitative tissue/sample inventory disagrees with coverage")
    return matched


def _annotation(inputs, manifest):
    reference, = manifest.references
    assets = {a.asset_id: a for a in manifest.sources}
    if (reference.taxon_id != 9606 or reference.annotation_name != "ensembl"
            or reference.annotation_version != "83" or reference.source_sha256 != assets["annotation.gtf.gz"].sha256):
        raise ValueError("Healthy human tissue requires its declared Ensembl 83 reference")
    genes, transcripts, proteins, header = set(), {}, {}, {}
    with gzip.open(inputs["annotation.gtf.gz"], "rt", encoding="utf-8") as handle:
        for text in handle:
            if text.startswith("#!"):
                key, _, value = text[2:].rstrip().partition(" ")
                header[key] = value
            elif not text.startswith("#"):
                columns = text.rstrip().split("\t")
                if len(columns) != 9:
                    raise ValueError("Malformed native Ensembl GTF record")
                attrs = dict(re.findall(r'(\w+) "([^\"]*)";', columns[8]))
                genes.add(attrs["gene_id"])
                if "transcript_id" in attrs:
                    transcripts[attrs["transcript_id"]] = attrs["gene_id"]
                if "protein_id" in attrs:
                    proteins[attrs["protein_id"]] = (attrs["transcript_id"], attrs["gene_id"])
    accession = header.get("genome-build-accession", "").removeprefix("NCBI:")
    if (header.get("genome-version") != reference.genome_build or accession != reference.assembly_accession
            or header.get("genome-build") != manifest.processing["source_assembly_header"]):
        raise ValueError("Native GTF assembly does not match the declared reference")
    actual = dict(genes=len(genes), transcripts=len(transcripts), proteins=len(proteins))
    if actual != dict(manifest.coverage["annotation_counts"]):
        raise ValueError("Historical annotation inventory changed")
    spec, = search_references(manifest)
    search = SearchReference(inputs["search.fasta.gz"], spec)
    by_protein, search_genes = {}, set()
    for entry in search.entries():
        identifier = entry.fasta_id.rsplit(".", 1)[0]
        attrs = dict(field.split(":", 1) for field in entry.header.split()[1:] if ":" in field)
        gene = attrs["gene"].rsplit(".", 1)[0]
        transcript = attrs["transcript"].rsplit(".", 1)[0]
        if identifier in by_protein:
            raise ValueError("Multiple native versions prevent unversioned protein assignment")
        by_protein[identifier] = dict(entry=entry, gene=gene, transcript=transcript)
        search_genes.add(gene)
    return reference, genes, transcripts, proteins, by_protein, search_genes


def _target(reference, kind, identifier):
    return ScopedTarget(reference, kind, identifier,
                        "Wang2019:EV1-protein-group" if kind == "protein_group" else "ensembl")


def import_tissue(inputs, output, manifest):
    """Import the authors' six abundance tables without attributing groups to ORFs."""
    from ..dataset import QuantitativeDataset
    from ..measurements import Assay, Measurement, Sample
    from ..pairing import AssayPair, TargetRelation
    from ..reference_mapping import ReferenceMapping

    reference, genes, transcripts, proteins, search, search_genes = _annotation(inputs, manifest)
    metadata = _samples(inputs, manifest)
    groups = tuple(_rows(inputs, manifest, "protein_groups"))
    gene_rows = tuple(_rows(inputs, manifest, "protein_genes"))
    group_by_copy_id = {}
    for line, identifiers, values in groups:
        native, representatives, gene_ids, label = identifiers
        if representatives in group_by_copy_id:
            raise ValueError("Protein-copy IDs cannot resolve to one native group row")
        members = representatives.split(";")
        if any(member not in search for member in members):
            raise ValueError("Protein representative is absent from the native Ensembl search")
        if any(gene not in search_genes for gene in gene_ids.split(";")):
            raise ValueError("Protein group gene is absent from the native Ensembl search")
        group_by_copy_id[representatives] = (line, identifiers, dict((t, v) for t, _, v in values))
    if any(ids[0] not in search_genes for _, ids, _ in gene_rows):
        raise ValueError("Protein gene is absent from the native Ensembl search")

    sample_ids = {tissue: record["Tissue ID"] for tissue, record in metadata.items()}
    samples, assays, pairings = [], [], []
    series = ("rna_transcripts", "rna_genes", "protein_groups", "protein_genes")
    normalization = {
        "rna_transcripts": "source Cufflinks 2.1.1 transcript FPKM; TopHat 2.0.8b; native Ensembl 83",
        "rna_genes": "source gene aggregation and median normalization; reported linear FPKM; methods state a 1-FPKM cutoff, but its normalized per-tissue value is not supplied",
        "protein_groups": "source MaxQuant 1.5.3.30 group iBAQ; matching between runs enabled",
        "protein_genes": "source gene aggregation and log10 median centering; table reports linear normalized iBAQ",
    }
    for tissue, record in metadata.items():
        sid = sample_ids[tissue]
        provenance = dict(source_asset="protein.zip", member="Table_EV1.xlsx",
                          sheet="A. Sample information", row=record["source_line"])
        samples.append(Sample(sid, reference, biological_unit_id=sid, tissue=tissue,
                              condition=record["Tissue status"],
                              attributes=dict(source_sample=record, donor_id=None,
                                              specimen_level="adjacent cryosections of one source tissue specimen",
                                              sex=None if record["Gender of donor"] == "-" else record["Gender of donor"],
                                              age_years=None if record["Age of donor"] == "-" else record["Age of donor"]),
                              provenance=provenance))
        for name in series:
            rna = name.startswith("rna_")
            calibration = dict(absolute_standard_measurement=False)
            if name == "rna_transcripts":
                calibration.update(copies_method="source RNA copies from FPKM and estimated cellular RNA",
                                   assumed_mrna_fraction_of_total_rna_mass=0.02,
                                   cellular_rna_estimation="total ribosomal-protein intensity; source proteomic-ruler method",
                                   source_detection_threshold_fpkm=1)
            elif name == "protein_groups":
                calibration.update(copies_method="source proteomic ruler from intensity values",
                                   calibration_not_supplied="per-sample calibration constants and uncertainty")
            assays.append(Assay(f"{name}:{sid}", sid, "rna_seq" if rna else "digest_ms", reference,
                                run_id=record["RNA-Seq Data ID"] if rna else record["MS Experiment ID"],
                                search_reference=reference.key,
                                normalization=dict(method=normalization[name]), calibration=calibration,
                                attributes=dict(source_view=name, additional_biological_replicate=False),
                                provenance=provenance))
        for left in series[:2]:
            for right in series[2:]:
                pairings.append(AssayPair(f"{left}:{right}:{sid}", reference, manifest.dataset_id, manifest.version,
                                           f"{left}:{sid}", f"{right}:{sid}", "adjacent_donor_tissue",
                                           dict(provenance, methods="adjacent cryosections of the same histologically healthy specimen",
                                                comparison_views="alternative summaries of the same acquisitions; not extra biological replicates"),
                                           "concurrent"))

    def relations():
        for line, identifiers, _ in gene_rows:
            target = _target(reference, "gene", identifiers[0])
            yield TargetRelation(f"gene:{line}", manifest.dataset_id, manifest.version, target, target,
                                 "same_target", "one_to_one",
                                 dict(source_asset="protein.zip", sheet="C. Genes", row=line,
                                      scope="source gene-level protein summary; not an individual protein or ORF"))
        for line, identifiers, _ in groups:
            target = _target(reference, "protein_group", f"row:{line}")
            parents = {("gene", gene) for gene in identifiers[2].split(";")}
            parents.update(("transcript", search[pid]["transcript"]) for pid in identifiers[1].split(";"))
            for kind, identifier in sorted(parents):
                yield TargetRelation(f"group:{line}:{kind}:{identifier}", manifest.dataset_id, manifest.version,
                                     _target(reference, kind, identifier), target,
                                     "gene_protein_group_context", "many_to_many",
                                     dict(source_asset="protein.zip", sheet="B. Protein groups", row=line,
                                          method="source representative/gene fields and pinned FASTA transcript headers",
                                          scope="candidate group attribution, not individual protein quantification"))

    def mappings():
        for line, identifiers, _ in groups:
            yield ReferenceMapping(f"group:{line}", manifest.dataset_id, manifest.version,
                                   _target(reference, "protein_group", f"row:{line}"), reference,
                                   tuple(_target(reference, "protein", pid) for pid in identifiers[1].split(";")),
                                   "unresolved", "source representative IDs checked against native Ensembl 83 FASTA",
                                   dict(native_protein_ids=identifiers[0], native_representatives=identifiers[1],
                                        native_gene_ids=identifiers[2], native_gene_name=identifiers[3],
                                        representative_search_records=[dict(fasta_id=search[pid]["entry"].fasta_id,
                                                                           sequence_sha256=search[pid]["entry"].sequence_sha256)
                                                                       for pid in identifiers[1].split(";")],
                                        unique_peptide_attribution="not supplied in this table"),
                                   reason="Protein-group representatives are candidates, not independent protein measurements",
                                   provenance=dict(source_asset="protein.zip", member="Table_EV1.xlsx",
                                                   sheet="B. Protein groups", row=line))

    def measurements():
        transcript_rows = {}
        # Inputs precede source-derived copy estimates, preserving acyclic lineage.
        order = ("rna_transcripts", "rna_genes", "protein_groups", "protein_genes", "protein_copies", "rna_copies")
        for name in order:
            zeros = numeric_text = count = 0
            rows = groups if name == "protein_groups" else gene_rows if name == "protein_genes" else _rows(inputs, manifest, name)
            asset, member, sheet, _ = SERIES[name]
            for line, identifiers, values in rows:
                assay_series = name
                inputs_by_tissue = {}
                associated_by_tissue = {}
                derivation = None
                if name == "rna_transcripts":
                    transcript, protein, gene, _ = identifiers
                    if transcripts.get(transcript) != gene or proteins.get(protein) != (transcript, gene):
                        raise ValueError("Native RNA transcript/protein/gene relationship disagrees with Ensembl 83")
                    target = _target(reference, "transcript", transcript)
                    transcript_rows[transcript] = line
                elif name == "rna_genes":
                    if identifiers[0] not in genes:
                        raise ValueError("Native RNA gene is absent from the historical GTF")
                    target = _target(reference, "gene", identifiers[0])
                elif name == "protein_groups":
                    target = _target(reference, "protein_group", f"row:{line}")
                elif name == "protein_genes":
                    target = _target(reference, "gene", identifiers[0])
                elif name == "protein_copies":
                    if identifiers[0] not in group_by_copy_id:
                        raise ValueError("Protein copies lack their corresponding native group")
                    source_line, source_ids, source_values = group_by_copy_id[identifiers[0]]
                    if tuple(identifiers[:2]) != tuple(source_ids[1:3]):
                        raise ValueError("Protein-copy identifiers disagree with the corresponding native group")
                    target = _target(reference, "protein_group", f"row:{source_line}")
                    assay_series = "protein_groups"
                    associated_by_tissue = {t: f"protein_groups:{source_line}:{sample_ids[t]}" for t, _, _ in values}
                    if any((_number(v) == 0) != (_number(source_values[t]) == 0) for t, _, v in values):
                        raise ValueError("Protein-copy detection disagrees with native group iBAQ")
                    derivation = dict(method="author-reported proteomic-ruler estimate from intensity values",
                                      assumptions="source calibration; not isotope-standard absolute quantification",
                                      uncertainty="not supplied", calibration_constants="not supplied",
                                      unavailable_inputs="raw intensity inputs are not present in these processed tables; iBAQ is an associated estimate, not the ruler input")
                else:
                    transcript, gene, _ = identifiers
                    if transcripts.get(transcript) != gene or transcript not in transcript_rows:
                        raise ValueError("Transcript copies lack their native transcript/gene input")
                    target = _target(reference, "transcript", transcript)
                    assay_series = "rna_transcripts"
                    inputs_by_tissue = {t: f"rna_transcripts:{transcript_rows[transcript]}:{sample_ids[t]}" for t, _, _ in values}
                    derivation = dict(method="author-reported transcript copies from FPKM and estimated cellular RNA",
                                      assumed_mrna_fraction_of_total_rna_mass=0.02,
                                      cellular_rna_estimation="ribosomal protein intensities/proteomic ruler",
                                      source_detection_threshold_fpkm=1, uncertainty="not supplied")
                rna = name.startswith("rna_")
                unit = "copies/cell" if name.endswith("copies") else "FPKM" if rna else "iBAQ"
                for tissue, column, native in values:
                    value = _number(native)
                    zeros += value == 0
                    numeric_text += isinstance(native, str)
                    count += 1
                    status = "measured"
                    if value == 0 and name != "rna_transcripts":
                        status = "below_detection" if rna else "not_measured"
                    sid = sample_ids[tissue]
                    kwargs = {}
                    if derivation is not None:
                        kwargs.update(origin="derived", derivation=derivation,
                                      input_measurements=(inputs_by_tissue[tissue],) if inputs_by_tissue else ())
                    provenance = dict(source_asset=asset, member=member, sheet=sheet, row=line,
                                      column=column, native_value=native,
                                      zero_semantics="source censored/unquantified value" if status != "measured" else "native quantitative estimate")
                    if associated_by_tissue:
                        provenance.update(associated_iBAQ_measurement=associated_by_tissue[tissue],
                                          native_protein_ids=identifiers[0], native_gene_ids=identifiers[1],
                                          native_gene_label=identifiers[2])
                    yield Measurement(f"{name}:{line}:{sid}", target, sid, f"{assay_series}:{sid}",
                                      "rna_abundance" if rna else "protein_abundance",
                                      value if status == "measured" else None, unit,
                                      manifest.dataset_id, manifest.version, status=status,
                                      provenance=provenance,
                                      **kwargs)
            expected = manifest.coverage["series"][name]
            if (count != expected["rows"] * len(metadata) or zeros != expected["zero_cells"]
                    or numeric_text != expected["numeric_text_cells"]):
                raise ValueError("Native abundance cell coverage changed")

    return QuantitativeDataset.create(output, manifest, samples=samples, assays=assays,
                                      pairings=pairings, relations=relations(), reference_mappings=mappings(),
                                      measurements=measurements())
