"""Pinned NCBI genome-download feature tables and reference proteins."""

import csv
import gzip
import hashlib
import io
import re
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from .proteome import Occurrence, ProteomeEntry, _sha256
from .reference import ReferenceContext

_ACCESSION = re.compile(r"[A-Z]{2}_[0-9]+\.[1-9][0-9]*")


def _source(path, expected):
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("RefSeq sources require lowercase SHA-256 byte checksums")
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"RefSeq source checksum mismatch: {path}")
    return (gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw).decode("ascii")


def _proteins(text):
    result, accession, chunks = {}, None, []
    for line in text.splitlines() + [">"]:
        if line.startswith(">"):
            if accession is not None:
                sequence = "".join(chunks).removesuffix("*")
                if not sequence or not re.fullmatch(r"[A-Z]+", sequence):
                    raise ValueError(f"Invalid RefSeq protein sequence: {accession}")
                if accession in result:
                    raise ValueError(f"Duplicate RefSeq protein accession: {accession}")
                result[accession] = sequence
            accession, chunks = (line[1:].split() or [None])[0], []
            if accession is not None and not _ACCESSION.fullmatch(accession):
                raise ValueError("RefSeq proteins require versioned accessions")
        elif line.strip():
            if accession is None:
                raise ValueError("RefSeq FASTA sequence precedes its header")
            chunks.append(line.strip())
    if not result:
        raise ValueError("RefSeq protein inventory is empty")
    return result


def refseq_proteome(cls, protein_fasta, feature_table, *, reference, expected_sha256):
    if (not isinstance(reference, ReferenceContext) or reference.annotation_name != "NCBI RefSeq"
            or not reference.annotation_version or not reference.assembly_accession):
        raise ValueError("RefSeq requires explicit taxon, assembly accession and NCBI RefSeq annotation")
    if not isinstance(expected_sha256, dict) or set(expected_sha256) != {"protein", "feature_table"}:
        raise ValueError("RefSeq expected_sha256 requires protein and feature_table")
    # Read exactly the verified bytes. Subsequent file changes cannot affect this build.
    proteins = _proteins(_source(protein_fasta, expected_sha256["protein"]))
    text = _source(feature_table, expected_sha256["feature_table"])
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {"# feature", "class", "assembly", "product_accession", "related_accession",
                "GeneID", "symbol", "genomic_accession", "strand", "product_length", "start", "end"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise ValueError("Missing NCBI genome-download feature-table columns")
    if len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("Duplicate RefSeq feature-table columns")
    occurrences, transcripts = defaultdict(dict), {}
    untranslated = 0
    for row in reader:
        if row["assembly"] != reference.assembly_accession:
            raise ValueError("RefSeq feature-table assembly disagrees with reference")
        if row["# feature"] != "CDS":
            continue
        protein, transcript, gene = row["product_accession"], row["related_accession"], row["GeneID"]
        if row["class"] == "without_protein" and not protein and not transcript:
            untranslated += 1
            continue
        if (row["class"] != "with_protein" or not _ACCESSION.fullmatch(protein)
                or not _ACCESSION.fullmatch(transcript) or not re.fullmatch(r"[1-9][0-9]*", gene)
                or not row["genomic_accession"] or row["strand"] not in {"+", "-"}):
            raise ValueError("Incomplete RefSeq CDS accession/gene mapping")
        if protein not in proteins:
            raise ValueError(f"Mapped RefSeq protein absent from FASTA: {protein}")
        if row["product_length"] != str(len(proteins[protein])):
            raise ValueError(f"RefSeq feature-table protein length disagrees with FASTA: {protein}")
        assignment = (gene, protein)
        if transcript in transcripts and transcripts[transcript] != assignment:
            raise ValueError(f"Conflicting RefSeq transcript assignment: {transcript}")
        transcripts[transcript] = assignment
        if (not re.fullmatch(r"[1-9][0-9]*", row["start"])
                or not re.fullmatch(r"[1-9][0-9]*", row["end"])
                or int(row["start"]) > int(row["end"])):
            raise ValueError("Invalid RefSeq source genomic span")
        occurrence = Occurrence(transcript, gene, row["symbol"] or None, None,
                                row["genomic_accession"], row["strand"], None, None, (), True)
        key = (transcript, gene, row["genomic_accession"], row["strand"], int(row["start"]), int(row["end"]))
        if key in occurrences[protein] and occurrences[protein][key] != occurrence:
            raise ValueError(f"Conflicting RefSeq source occurrence: {transcript}")
        occurrences[protein][key] = occurrence
    missing = set(proteins) - set(occurrences)
    if missing:
        raise ValueError("RefSeq FASTA proteins lack feature-table mappings: " + ", ".join(sorted(missing)[:5]))
    source_hash, sequence_hash = _sha256(expected_sha256), _sha256(sorted(proteins.items()))
    if (reference.source_sha256 not in (None, source_hash)
            or reference.sequence_sha256 not in (None, sequence_hash)):
        raise ValueError("RefSeq source/sequence checksums disagree with reference")
    reference = replace(reference, source_sha256=source_hash, sequence_sha256=sequence_hash)
    source = f"NCBI RefSeq:{reference.annotation_version}"
    entries = [ProteomeEntry(accession, proteins[accession], "annotated", reference.taxon_id,
                source, occurrences=tuple(occurrences[accession][key] for key in sorted(occurrences[accession])),
                source_binding="feature-table-accession-join; spliced CDS coordinates unavailable")
               for accession in sorted(proteins)]
    layer = dict(kind="annotated", role="host", source=source, reference=reference.key,
                 added=len(entries), content_provenance=dict(expected_sha256),
                 identifier_namespaces={"gene": "NCBI Gene", "transcript": "NCBI RefSeq RNA",
                                        "protein": "NCBI RefSeq protein"},
                 cds_without_protein=untranslated, mapping="feature-table-accession-join")
    # Same accession can have multiple genomic placements. Keep their spans
    # indexed exactly like entry.occurrences without inventing spliced CDS blocks.
    layer['genomic_spans'] = {accession: [list(key[-2:]) for key in sorted(occurrences[accession])]
                             for accession in sorted(occurrences)}
    return cls(reference, entries, (layer,))
