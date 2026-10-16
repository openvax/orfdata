"""Ribo-seq ORFs of Yang et al. 2024 (RibORF 2.0), placed through Ensembl's transcripts or genome.

Yang, Li, Stroup, Wang and Ji, "Widespread stable noncanonical peptides
identified by integrated analyses of ribosome profiling and ORF features", Nat
Commun 15:1932 (2024), list the translated ORFs RibORF 2.0 identified in the
Ribo-seq data of their Supplementary Data 1 as Supplementary Data 3, one sheet
per species: human (hg38, GENCODE V39), mouse (mm10, GENCODE vM25),
zebrafish (GRCz11), worm (WBcel235) and yeast (sacCer3). Each row names its
transcript and its transcript coordinates: a 1-based start and an exclusive
end, the stop codon included.

The source gives no protein sequences, so each ORF is placed through the same
transcript in a pinned Ensembl annotation (:func:`import_orfs`): its
transcript coordinates map through that transcript's exons, and its protein is
translated from that transcript's sequence. For mouse this also moves the ORFs
from mm10 to GRCm39 without a chain file with usage restrictions. A row is
left out, with the reason, when the transcript is missing, has no sequence or
changed length since the source's annotation, or when its sequence there
doesn't start with the stated codon and end at its first stop.

Yang et al.'s yeast transcripts carry UTRs and include unstable transcripts
(CUTs, SUTs) that Ensembl's yeast annotation lacks, so yeast ORFs are placed
from each row's genomic span instead (:func:`import_genome_orfs`), translated
from the pinned genome.

Annotated CDSs (``canonical``), likely sequencing errors and stop-codon
readthrough are not imported, nor are proteins shorter than
``processing["min_protein_length"]``: Chang et al.'s shortest is 5 residues,
and the source has thousands of 1-4 residue ones, bare start-stops among them.
Types Chang et al. name differently (ouORF, odORF, noncoding) take Chang's
names, so a query by type finds both sources' ORFs. IDs are ``riborf_`` and the
start and stop, never a name Chang et al. use. Presence in the source's set is
aggregate Ribo-seq support: it names no study and counts no libraries.
"""

import gzip
import hashlib
import json
from collections import Counter

from ..catalog import ORFCatalog, ORFEvidence, ORFRecord
from ..evidence_summary import EvidenceSummary
from ..model import ORFModel
from ..reference import ScopedTarget
from ..sequence import genetic_code
from .ncorf_mammals import _blocks, _transcripts

PRESENCE = ("In the source's genome-wide set of RibORF 2.0-identified translated ORFs (Supplementary Data 3), "
            "from the mouse Ribo-seq data of Supplementary Data 1 (224 samples, 26 GEO series). Presence, not a "
            "count; no study is named per ORF.")
STOPS = {"TAA", "TAG", "TGA"}


def _rows(path, sheet):
    try:
        from openpyxl import load_workbook
    except ImportError as error:
        raise ImportError("Install orfdata[build-data] to build this catalogue from its spreadsheet") from error

    with open(path, "rb") as handle:  # Cached files have no extension for openpyxl to go by.
        workbook = load_workbook(handle, read_only=True)
        rows = [row for row in workbook[sheet].iter_rows(values_only=True) if row and row[3]]
        workbook.close()
    return rows


def _sequences(paths, wanted):
    """Wanted sequences from FASTA files, named by their ID with or without its version."""
    found, name = {}, None
    for path in paths:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    # Names with dots of their own, such as C. elegans 2L52.1a.1, are matched whole first.
                    name = name if name in wanted else name.rsplit(".", 1)[0]
                    name = name if name in wanted else None
                    if name is not None:
                        found[name] = []
                elif name is not None:
                    found[name].append(line.strip())
    return {k: "".join(v) for k, v in found.items()}


def _versioned(stable, entry):
    return f"{stable}.{entry['version']}" if entry["version"] else stable


def _parsed(processing, sheet_rows):
    """The rows to import, as (transcript, length, start, end, type, start codon, row), and how many were short."""
    skip_types, names = set(processing["skip_types"]), processing["biotypes"]
    parsed, short = [], 0
    for row in sheet_rows:
        transcript, _, coordinates, orf_type, start_codon = row[3].split("|")
        if orf_type in skip_types:
            continue
        length, start, end = (int(x) for x in coordinates.split(":"))
        if (end - start) // 3 - 1 < processing["min_protein_length"]:
            short += 1
            continue
        parsed.append((transcript, length, start, end, names.get(orf_type, orf_type), start_codon, row))
    return parsed, short


def import_orfs(inputs, output, manifest):
    """Write the catalogue ``manifest`` describes from the source table and the pinned Ensembl files."""
    processing = manifest.processing
    reference = manifest.references[0]
    # Ensembl IDs carry a version (ENST00000000233.10); C. elegans names have dots of their own (2L52.1a.1).
    versioned = processing.get("versioned_ids", True)
    source_annotation = processing.get("source_annotation", "GENCODE vM25")
    parsed, short = _parsed(processing, _rows(inputs["table"], processing["sheet"]))

    def stable_id(source_transcript):
        name = source_transcript.split(":")[0]
        return name.split(".")[0] if versioned else name

    wanted = {stable_id(p[0]) for p in parsed}
    transcripts = _transcripts(inputs["gtf"], wanted)
    sequences = _sequences([inputs["cdna"], inputs["ncrna"]], wanted)
    code = genetic_code(1)

    paths, unresolved = {}, []
    for source_transcript, length, start, end, orf_type, start_codon, _ in parsed:
        stable = stable_id(source_transcript)

        def skip(reason):
            unresolved.append(dict(transcript=source_transcript.split(":")[0], start=start, end=end, reason=reason))

        entry, sequence = transcripts.get(stable), sequences.get(stable)
        if entry is None:
            skip(f"Transcript absent from Ensembl {reference.annotation_version}")
            continue
        if sequence is None:  # TEC transcripts, for one, are in the GTF but in neither FASTA.
            skip(f"Transcript has no sequence in Ensembl {reference.annotation_version}'s cDNA or ncRNA FASTA")
            continue
        if len(sequence) != length or sum(b - a + 1 for a, b in entry["exons"]) != length:
            skip(f"Transcript changed length since {source_annotation} ({length} nt, now {len(sequence)})")
            continue
        coding = sequence[start - 1 : end - 1]  # 1-based start, exclusive end.
        if set(coding) - set("ACGT"):
            skip("The transcript's sequence has ambiguous bases in the ORF")
            continue
        if (len(coding) % 3 or coding[:3] != start_codon or coding[-3:] not in STOPS
                or any(coding[i : i + 3] in STOPS for i in range(0, len(coding) - 3, 3))):
            skip("The transcript's sequence doesn't start with the stated codon and end at its first stop")
            continue
        declared = None if start_codon in code.start_codons else start_codon
        protein = code.translate(coding, start_codon=declared)[:-1]
        blocks = _blocks(entry, start, end - 1)
        key = (entry["contig"], entry["strand"], blocks)
        found = paths.setdefault(key, dict(protein=protein, start_codon=declared, rows=[]))
        found["rows"].append((stable, orf_type))

    models = {key: (_versioned(stable, transcripts[stable]), transcripts[stable]["gene_id"],
                    transcripts[stable]["gene_name"], orf_type, found["protein"], found["start_codon"])
              for key, found in paths.items() for stable, orf_type in [min(found["rows"])]}
    counts = dict(source_rows=len(parsed) + short, short_rows=short, resolved_rows=len(parsed) - len(unresolved),
                  unresolved_rows=len(unresolved))
    _write(models, unresolved, counts, output, manifest)


def import_genome_orfs(inputs, output, manifest):
    """Write the catalogue from each row's genomic span, translated from the pinned genome.

    For a source whose transcripts the pinned annotation lacks (Yang et al.'s
    yeast transcripts carry UTRs, and some are unstable CUTs and SUTs). The
    rows are BED lines whose span is the ORF, stop codon included. An ORF on
    several exons is left out: the source gives its block sizes but not where
    each block starts. Each record names the Ensembl transcript of the
    source's transcript where there is one (``processing["transcript_suffix"]``
    joins a systematic name to it), else the source's own transcript name.
    """
    processing = manifest.processing
    parsed, short = _parsed(processing, _rows(inputs["table"], processing["sheet"]))
    suffix = processing.get("transcript_suffix", "")
    names = {p[0].split(":")[0] for p in parsed}
    transcripts = _transcripts(inputs["gtf"], {name + suffix for name in names})
    genome = _genome(inputs["genome"], {str(p[6][0]) for p in parsed})
    code = genetic_code(1)

    paths, unresolved = {}, []
    for source_transcript, _, start, end, orf_type, start_codon, row in parsed:
        name = source_transcript.split(":")[0]
        contig, first, last, strand, block_count = str(row[0]), int(row[6]) + 1, int(row[7]), row[5], int(row[9])

        def skip(reason):
            unresolved.append(dict(transcript=name, start=start, end=end, reason=reason))

        if block_count != 1:
            skip("Spliced ORF: the source gives its block sizes but not where each block starts")
            continue
        if contig not in genome:
            skip("Contig absent from the pinned genome")
            continue
        coding = genome[contig][first - 1 : last]
        if strand == "-":
            coding = _reverse_complement(coding)
        if set(coding) - set("ACGT"):
            skip("The genome's sequence has ambiguous bases in the ORF")
            continue
        if (len(coding) % 3 or coding[:3] != start_codon or coding[-3:] not in STOPS
                or any(coding[i : i + 3] in STOPS for i in range(0, len(coding) - 3, 3))):
            skip("The genome's sequence doesn't start with the stated codon and end at its first stop")
            continue
        declared = None if start_codon in code.start_codons else start_codon
        protein = code.translate(coding, start_codon=declared)[:-1]
        key = (contig, strand, ((first, last),))
        found = paths.setdefault(key, dict(protein=protein, start_codon=declared, rows=[]))
        found["rows"].append((name, orf_type))

    models = {}
    for key, found in paths.items():
        name, orf_type = min(found["rows"])
        entry = transcripts.get(name + suffix)
        models[key] = ((name + suffix, entry["gene_id"], entry["gene_name"]) if entry is not None
                       else (name, name, "")) + (orf_type, found["protein"], found["start_codon"])
    counts = dict(source_rows=len(parsed) + short, short_rows=short, resolved_rows=len(parsed) - len(unresolved),
                  unresolved_rows=len(unresolved))
    _write(models, unresolved, counts, output, manifest)


def _genome(path, contigs):
    """The wanted contigs of a gzipped FASTA, named by their first word."""
    found, name = {}, None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                name = name if name in contigs else None
                if name is not None:
                    found[name] = []
            elif name is not None:
                found[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in found.items()}


def _reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def _write(models, unresolved, counts, output, manifest):
    """Records named by start and stop, checked against the manifest's coverage, as a gzipped catalogue.

    ``models`` maps each (contig, strand, blocks) to its transcript, gene ID and
    name, type, protein and declared start codon.
    """
    processing = manifest.processing
    reference = manifest.references[0]
    table = next(s for s in manifest.sources if s.asset_id == "table")
    presence = processing.get("presence", PRESENCE)

    def locus(contig, strand, blocks):
        first, last = (blocks[0][0], blocks[-1][1]) if strand == "+" else (blocks[-1][1], blocks[0][0])
        return f"riborf_{contig}_{first}_{last}_{strand}"

    # Named by start and stop; paths that share both also name their transcript, and answer to the shared name.
    shared = Counter(locus(*key) for key in models)
    records = []
    for (contig, strand, blocks), (transcript, gene_id, gene_name, orf_type, protein, start_codon) in sorted(
            models.items()):
        orf_id = locus(contig, strand, blocks)
        if shared[orf_id] > 1:
            orf_id = f"{orf_id}@{transcript}"
        model = ORFModel(orf_id, transcript, gene_id, gene_name, contig, strand, blocks, protein, orf_type,
                         reference.genome_build, reference.taxon_id, reference.species, source=processing["source"],
                         reference_context=reference, start_codon=start_codon)
        summary = EvidenceSummary(processing["source"], processing["release"],
                                  ScopedTarget(reference, "orf", orf_id, processing["source"]),
                                  hashlib.sha256(protein.encode("ascii")).hexdigest(), "ribosome",
                                  "RibORF 2.0 translated ORF", 1, presence, "orf", table.url, table.sha256)
        records.append(ORFRecord.from_model(model, ORFEvidence(summaries=(summary,))))

    for key, value in counts.items():
        if manifest.coverage.get(key) != value:
            raise ValueError(f"Source {key} is {value}, but the manifest declares {manifest.coverage.get(key)}")
    metadata = dict(
        genome_build=reference.genome_build, taxon_id=reference.taxon_id, species=reference.species,
        coordinates="1-based-inclusive", dataset_version=processing["release"],
        sources=[dict(id=processing["source"], release=processing["release"], url=table.url,
                      scope=processing["scope"])],
        unresolved=unresolved, unresolved_reasons=dict(Counter(u["reason"].split(" (")[0] for u in unresolved)),
        catalog_manifest=manifest.to_dict())
    catalog = ORFCatalog(records, metadata)
    raw = json.dumps(catalog.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    with open(output, "wb") as handle, gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
        compressed.write(raw)
