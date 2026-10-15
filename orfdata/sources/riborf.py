"""Ribo-seq ORFs of Yang et al. 2024 (RibORF 2.0), placed on GRCm39 through Ensembl's transcripts.

Yang, Li, Stroup, Wang and Ji, "Widespread stable noncanonical peptides
identified by integrated analyses of ribosome profiling and ORF features", Nat
Commun 15:1932 (2024), list the translated ORFs RibORF 2.0 identified in the
mouse Ribo-seq data of their Supplementary Data 1 (224 samples in 26 GEO
series) as Supplementary Data 3, in mm10 / GENCODE vM25 coordinates. Each row names its transcript and its
transcript coordinates: a 1-based start and an exclusive end, the stop codon
included.

The source gives no protein sequences, and lifting mm10 would need a chain
file with usage restrictions, so each ORF is placed on GRCm39 through the same
transcript in the pinned Ensembl annotation: its transcript coordinates map
through that transcript's exons, and its protein is translated from that
transcript's sequence. A row is left out, with the reason, when the transcript
is missing, has no sequence or changed length since GENCODE vM25, or when its
sequence there doesn't start with the stated codon and end at its first stop.

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
    found, name = {}, None
    for path in paths:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split()[0].split(".")[0]
                    name = name if name in wanted else None
                    if name is not None:
                        found[name] = []
                elif name is not None:
                    found[name].append(line.strip())
    return {k: "".join(v) for k, v in found.items()}


def import_orfs(inputs, output, manifest):
    """Write the catalogue ``manifest`` describes from the source table and the pinned Ensembl files."""
    processing = manifest.processing
    reference = manifest.references[0]
    table = next(s for s in manifest.sources if s.asset_id == "table")
    skip_types, names = set(processing["skip_types"]), processing["biotypes"]
    parsed, short = [], 0
    for row in _rows(inputs["table"], processing["sheet"]):
        transcript, _, coordinates, orf_type, start_codon = row[3].split("|")
        if orf_type in skip_types:
            continue
        length, start, end = (int(x) for x in coordinates.split(":"))
        if (end - start) // 3 - 1 < processing["min_protein_length"]:
            short += 1
            continue
        parsed.append((transcript, length, start, end, names.get(orf_type, orf_type), start_codon))
    wanted = {p[0].split(":")[0].split(".")[0] for p in parsed}
    transcripts = _transcripts(inputs["gtf"], wanted)
    sequences = _sequences([inputs["cdna"], inputs["ncrna"]], wanted)
    code = genetic_code(1)

    paths, unresolved = {}, []
    for source_transcript, length, start, end, orf_type, start_codon in parsed:
        stable = source_transcript.split(":")[0].split(".")[0]

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
            skip(f"Transcript changed length since GENCODE vM25 ({length} nt, now {len(sequence)})")
            continue
        coding = sequence[start - 1 : end - 1]  # 1-based start, exclusive end.
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

    def locus(contig, strand, blocks):
        first, last = (blocks[0][0], blocks[-1][1]) if strand == "+" else (blocks[-1][1], blocks[0][0])
        return f"riborf_{contig}_{first}_{last}_{strand}"

    # Named by start and stop; paths that share both also name their transcript, and answer to the shared name.
    shared = Counter(locus(*key) for key in paths)
    records = []
    for (contig, strand, blocks), found in sorted(paths.items()):
        stable, orf_type = min(found["rows"])
        entry = transcripts[stable]
        orf_id = locus(contig, strand, blocks)
        if shared[orf_id] > 1:
            orf_id = f"{orf_id}@{stable}.{entry['version']}"
        model = ORFModel(orf_id, f"{stable}.{entry['version']}", entry["gene_id"], entry["gene_name"], contig,
                         strand, blocks, found["protein"], orf_type, reference.genome_build, reference.taxon_id,
                         reference.species, source=processing["source"], reference_context=reference,
                         start_codon=found["start_codon"])
        summary = EvidenceSummary(processing["source"], processing["release"],
                                  ScopedTarget(reference, "orf", orf_id, processing["source"]),
                                  hashlib.sha256(found["protein"].encode("ascii")).hexdigest(), "ribosome",
                                  "RibORF 2.0 translated ORF", 1, PRESENCE, "orf", table.url, table.sha256)
        records.append(ORFRecord.from_model(model, ORFEvidence(summaries=(summary,))))

    counts = dict(source_rows=len(parsed) + short, short_rows=short, resolved_rows=len(parsed) - len(unresolved),
                  unresolved_rows=len(unresolved))
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
