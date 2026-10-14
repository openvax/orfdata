"""OpenProt alternative ORFs as an ORF catalogue, from its pinned TSV, BED and FASTA downloads.

Each TSV row assigns one OpenProt protein to one transcript; the BED file has
one row per TSV row, in the same order. Coding blocks come from the BED row,
whose 0-based half-open blocks include the stop codon. The TSV's genomic
start and stop are one lower than the BED's 1-based ones at both ends,
although its README calls them 1-based (issue #55). Every row must agree
with that, with its transcript span and with its protein's length, or the
import fails rather than guessing.

Rows on transcripts of the pinned annotation (Ensembl) become records; the
rest, such as RefSeq assignments, are listed in ``metadata["unresolved"]``
with the reason, since no annotation is pinned for them. One protein at one
spliced coding path is one record, named by its OpenProt accession; when an
accession has several paths, such as copies on chrY, each record's ID adds
``@`` and its first transcript. The MS score and translation-event (TE) score
become aggregate :class:`EvidenceSummary` counts for the protein sequence:
they are not peptides or named studies, and their counting rules are
OpenProt's.
"""

import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict

from ..catalog import ORFCatalog, ORFEvidence, ORFRecord
from ..evidence_summary import EvidenceSummary
from ..model import ORFModel
from ..reference import ScopedTarget

MS_DEFINITION = "Sum of study-specific unique-peptide counts; not global peptide uniqueness."
TE_DEFINITION = "Count of studies reporting one or more significant translation detections."
STRANDS = {"1": "+", "-1": "-"}


def _member(path, name):
    with zipfile.ZipFile(path) as archive:
        if name not in archive.namelist():
            raise ValueError(f"{path} lacks {name}")
        return archive.read(name).decode("utf-8")


def _fasta(text):
    """Sequences by every accession a header names: its own and its PA= aliases."""
    entries = []
    for line in text.splitlines():
        if line.startswith(">"):
            fields = line[1:].split()
            names = [fields[0].split("|")[0]] + [a for f in fields if f.startswith("PA=") for a in f[3:].split(",")]
            entries.append((names, []))
        elif entries:
            entries[-1][1].append(line.strip())
    proteins = {}
    for names, lines in entries:
        protein = "".join(lines)
        for name in names:
            if proteins.setdefault(name, protein) != protein:
                raise ValueError(f"OpenProt FASTA gives {name} two sequences")
    return proteins


def _blocks(row):
    """1-based inclusive coding blocks from a BED12 row."""
    start = int(row[1])
    sizes = [int(x) for x in row[10].strip(",").split(",")]
    offsets = [int(x) for x in row[11].strip(",").split(",")]
    if len(sizes) != int(row[9]) or len(offsets) != len(sizes) or start + offsets[-1] + sizes[-1] != int(row[2]):
        raise ValueError(f"Malformed BED blocks for {row[3]}")
    return tuple((start + o + 1, start + o + s) for o, s in zip(offsets, sizes))


def _summary(reference, accession, protein, method, metric, value, definition, url, sha256):
    return EvidenceSummary("OpenProt", reference.annotation_version,
                           ScopedTarget(reference, "protein", accession, "openprot"),
                           hashlib.sha256(protein.encode("ascii")).hexdigest(), method, metric, int(value),
                           definition, "protein_sequence", url, sha256)


def import_orfs(inputs, output, manifest):
    """Write the catalogue ``manifest`` describes from its OpenProt TSV, BED and FASTA archives."""
    processing = manifest.processing
    members = processing["members"]
    rows = list(csv.DictReader(io.StringIO(_member(inputs["tsv"], members["tsv"])), delimiter="\t"))
    bed = [line.split("\t") for line in _member(inputs["bed"], members["bed"]).splitlines() if line.strip()]
    proteins = _fasta(_member(inputs["fasta"], members["fasta"]))
    record_reference = next(r for r in manifest.references if r.annotation_name == processing["annotation_name"])
    summary_reference = next(r for r in manifest.references if r.annotation_name == "OpenProt")
    tsv = next(s for s in manifest.sources if s.asset_id == "tsv")
    if len(bed) != len(rows):
        raise ValueError("OpenProt TSV and BED row counts disagree")

    paths, unresolved, scores = defaultdict(dict), [], {}
    for line, (row, region) in enumerate(zip(rows, bed), start=2):
        accession, transcript = row["protein accession"], row["transcript accession"]
        if row["chr"] in ("MT", "M"):  # Needs table 2 and its starts, e.g. human IP_306399 begins ATA (#55).
            raise ValueError(f"OpenProt TSV line {line} ({accession}) is mitochondrial; the importer handles "
                             "nuclear ORFs only")
        length = int(row["protein length (a.a.)"])
        blocks = _blocks(region)
        protein = proteins.get(accession)
        if (region[3] != accession or region[0] != "chr" + row["chr"] or region[5] != STRANDS.get(row["strand"])
                or (blocks[0][0], blocks[-1][1]) != (int(row["start genomic coordinates"]) + 1,
                                                     int(row["stop genomic coordinates"]) + 1)
                or sum(b - a + 1 for a, b in blocks) != 3 * (length + 1)
                or int(row["stop transcript coordinates"]) - int(row["start transcript coordinates"]) != 3 * (length + 1)
                or protein is None or len(protein) != length):
            raise ValueError(f"OpenProt TSV line {line} ({accession} on {transcript}) disagrees with its "
                             "BED row, protein or the verified coordinate convention")
        score = (row["MS score"], row["TE Translation Event score"])
        if scores.setdefault(accession, score) != score:
            raise ValueError(f"OpenProt reports two sets of scores for {accession}")
        if not transcript.startswith(processing["transcript_prefix"]):
            unresolved.append(dict(protein=accession, transcript=transcript, reason=processing["unresolved_reason"]))
            continue
        key = (row["chr"], region[5], blocks)
        found = paths[accession].setdefault(key, dict(protein=protein, genes={}, types=set()))
        found["genes"][transcript] = row["gene symbol"].split(".")[0]  # Overlapping genes can share a path.
        found["types"].add(row["protein type"])

    records = []
    for accession, found in sorted(paths.items()):
        ms, te = scores[accession]
        for (contig, strand, blocks), path in sorted(found.items(), key=lambda item: min(item[1]["genes"])):
            transcript = min(path["genes"])
            if len(path["types"]) != 1:
                raise ValueError(f"OpenProt gives {accession} several protein types")
            protein = path["protein"]
            model = ORFModel(accession if len(found) == 1 else f"{accession}@{transcript}", transcript,
                             path["genes"][transcript], "", contig, strand, blocks, protein, path["types"].pop(),
                             record_reference.genome_build, record_reference.taxon_id, record_reference.species,
                             source="openprot", reference_context=record_reference)
            summaries = tuple(_summary(summary_reference, accession, protein, *values, tsv.url, tsv.sha256) for values in (
                ("ms", "MS score", ms, MS_DEFINITION), ("ribosome", "TE Translation Event score", te, TE_DEFINITION)))
            records.append(ORFRecord.from_model(model, ORFEvidence(summaries=summaries)))

    expected = dict(source_rows=len(rows), resolved_rows=len(rows) - len(unresolved), unresolved_rows=len(unresolved),
                    accessions=len(scores), resolved_accessions=len(paths))
    for key, value in expected.items():
        if manifest.coverage.get(key) != value:
            raise ValueError(f"OpenProt {key} is {value}, but the manifest declares {manifest.coverage.get(key)}")
    metadata = dict(
        genome_build=record_reference.genome_build, taxon_id=record_reference.taxon_id,
        species=record_reference.species, coordinates="1-based-inclusive",
        dataset_version=f"openprot-{summary_reference.annotation_version}",
        sources=[dict(id="openprot", release=summary_reference.annotation_version, url=tsv.url,
                      scope=processing["scope"])],
        unresolved=unresolved, unresolved_reasons=dict(Counter(u["reason"] for u in unresolved)),
        catalog_manifest=manifest.to_dict())
    catalog = ORFCatalog(records, metadata)
    raw = json.dumps(catalog.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    with open(output, "wb") as handle, gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
        compressed.write(raw)
