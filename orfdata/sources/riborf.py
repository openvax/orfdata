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
translated from that transcript's sequence. Mouse uses the separate native
GENCODE M25 reconstruction and exact-sequence placement in ``riborf_native``;
the unvalidated historical mouse rebuilds are rejected here. A row is
left out, with the reason, when the transcript is missing, has no sequence or
changed length since the source's annotation, or when its sequence there
doesn't start with the stated codon and end at its first stop.

Yang et al.'s yeast transcripts carry UTRs and include unstable transcripts
(CUTs, SUTs) that Ensembl's yeast annotation lacks, so yeast ORFs are placed
from each row's genomic span instead (:func:`import_genome_orfs`), translated
from the pinned genome.

Annotated CDSs (``canonical``), likely sequencing errors and stop-codon
readthrough are not imported. The current catalogue retains valid short source
ORFs; release 1 applied ``processing["min_protein_length"]`` during import.
Consumers can apply a length threshold explicitly when querying or building a
proteome.
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
    from openpyxl import load_workbook

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
    """Rows to import and counts before reference validation.

    RibORF catalogue release 1 used ``min_protein_length`` as an import filter.
    Release 2 reports rows below ``report_protein_lengths_below`` but retains
    them, so consumer-specific peptide length policy does not erase ORF models.
    """
    skip_types, names = set(processing["skip_types"]), processing["biotypes"]
    minimum = processing.get("min_protein_length")
    report_below = processing.get("report_protein_lengths_below", minimum)
    parsed, source_rows, short = [], 0, 0
    for row in sheet_rows:
        transcript, _, coordinates, orf_type, start_codon = row[3].split("|")
        if orf_type in skip_types:
            continue
        source_rows += 1
        length, start, end = (int(x) for x in coordinates.split(":"))
        protein_length = (end - start) // 3 - 1
        if report_below is not None and protein_length < report_below:
            short += 1
        if minimum is not None and protein_length < minimum:
            continue
        parsed.append((transcript, length, start, end, names.get(orf_type, orf_type), start_codon, row))
    return parsed, source_rows, short


def validate_transcript_manifest(manifest):
    """Reject superseded mouse builds before acquiring their source files."""
    if manifest.references[0].taxon_id == 10090:
        raise ValueError("Mouse RibORF generations 1/2 lack native sequence validation; build generation 3")


def import_orfs(inputs, output, manifest):
    """Write the catalogue ``manifest`` describes from the source table and the pinned Ensembl files."""
    validate_transcript_manifest(manifest)
    processing = manifest.processing
    reference = manifest.references[0]
    # Ensembl IDs carry a version (ENST00000000233.10); C. elegans names have dots of their own (2L52.1a.1).
    versioned = processing.get("versioned_ids", True)
    source_annotation = processing.get("source_annotation", "GENCODE vM25")
    parsed, source_rows, short = _parsed(processing, _rows(inputs["table"], processing["sheet"]))

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
    counts = dict(source_rows=source_rows, short_rows=short, resolved_rows=len(parsed) - len(unresolved),
                  unresolved_rows=len(unresolved))
    _write(models, unresolved, counts, output, manifest)


def import_genome_orfs(inputs, output, manifest):
    """Write the catalogue from each row's genomic span, translated from the pinned genome.

    For a source whose transcripts the pinned annotation lacks (Yang et al.'s
    yeast transcripts carry UTRs, and some are unstable CUTs and SUTs). The
    rows are BED lines whose span is the ORF, stop codon included, with block
    sizes but no block starts. An ORF on two exons is rebuilt from its first
    block's size at its start and its last block's at its end, when the gap
    between them is a GT-AG intron; one on more exons, or across a gap that
    isn't an intron (such as Ty1's programmed frameshift), is left out. A
    record names the Ensembl transcript of the source's transcript
    (``processing["transcript_suffix"]`` joins a systematic name to it) when
    that transcript carries the ORF's coding path, else the source's own
    transcript name, with the Ensembl gene's ID and name where there is one:
    Ensembl's yeast transcripts span only the CDS, so most uORFs and dORFs
    lie outside them.
    """
    processing = manifest.processing
    parsed, source_rows, short = _parsed(processing, _rows(inputs["table"], processing["sheet"]))
    suffix = processing.get("transcript_suffix", "")
    names = {p[0].split(":")[0] for p in parsed}
    # Ensembl names most yeast transcripts with the suffix, transposable elements without it.
    transcripts = _transcripts(inputs["gtf"], {name + suffix for name in names} | names)
    genome = _genome(inputs["genome"], {str(p[6][0]) for p in parsed})
    code = genetic_code(1)

    paths, unresolved = {}, []
    for source_transcript, _, start, end, orf_type, start_codon, row in parsed:
        name = source_transcript.split(":")[0]
        contig, first, last, strand, block_count = str(row[0]), int(row[6]) + 1, int(row[7]), row[5], int(row[9])

        def skip(reason):
            unresolved.append(dict(transcript=name, start=start, end=end, reason=reason))

        if contig not in genome:
            skip("Contig absent from the pinned genome")
            continue
        sequence = genome[contig]
        if block_count == 1:
            blocks = ((first, last),)
        elif block_count == 2 and (int(row[1]), int(row[2])) == (int(row[6]), int(row[7])):
            head, tail = (int(x) for x in str(row[10]).rstrip(",").split(","))
            blocks = ((first, first + head - 1), (last - tail + 1, last))
            gap = sequence[blocks[0][1] : blocks[1][0] - 1]
            gap = _reverse_complement(gap) if strand == "-" else gap
            if len(gap) < 2 or gap[:2] not in ("GT", "GC") or gap[-2:] != "AG":
                skip(f"The source's two blocks are joined across {len(gap)} nt that aren't an intron "
                     "(a programmed frameshift?)")
                continue
        else:
            skip("Spliced ORF on more than two exons: the source gives block sizes but not where they start")
            continue
        coding = "".join(sequence[a - 1 : b] for a, b in blocks)
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
        key = (contig, strand, blocks)
        found = paths.setdefault(key, dict(protein=protein, start_codon=declared, rows=[]))
        found["rows"].append((name, orf_type))

    models = {}
    for (contig, strand, blocks), found in paths.items():
        name, orf_type = min(found["rows"])
        ensembl = next((t for t in (name + suffix, name) if t in transcripts), None)
        entry = transcripts.get(ensembl)
        if entry is None:
            named = (name, name, "")
        elif _carries(entry, strand, blocks):
            named = (ensembl, entry["gene_id"], entry["gene_name"])
        else:  # On the source's longer model of the gene's transcript, outside Ensembl's.
            named = (name, entry["gene_id"], entry["gene_name"])
        models[contig, strand, blocks] = named + (orf_type, found["protein"], found["start_codon"])
    counts = dict(source_rows=source_rows, short_rows=short, resolved_rows=len(parsed) - len(unresolved),
                  unresolved_rows=len(unresolved))
    _write(models, unresolved, counts, output, manifest)


def _carries(entry, strand, blocks):
    """Whether a transcript's exons carry these coding blocks: each in one exon, joined only at splice sites."""
    if entry["strand"] != strand:
        return False
    exons = sorted(entry["exons"])
    found = []
    for a, b in blocks:
        index = next((i for i, (x, y) in enumerate(exons) if x <= a and b <= y), None)
        if index is None:
            return False
        found.append(index)
    return all(j == i + 1 and b == exons[i][1] and c == exons[j][0]
               for (_, b), (c, _), i, j in zip(blocks, blocks[1:], found, found[1:]))


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
    report_below = processing.get("report_protein_lengths_below")
    if "short_records" in manifest.coverage:
        counts["short_records"] = sum(len(value[-2]) < report_below for value in models.values())

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
