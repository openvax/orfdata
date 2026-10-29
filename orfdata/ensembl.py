"""Pinned Ensembl GTF/RNA/peptide accession joins without annotation indexes."""

import gzip
import hashlib
import io
import re
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, replace

from .proteome import Occurrence, ProteinSequence, _sha256
from .reference import ReferenceContext


class _HashedReader(io.RawIOBase):
    """Hash the original bytes consumed by the decompressor, not decoded text."""

    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()

    def readable(self):
        return True

    def readinto(self, buffer):
        size = self.stream.readinto(buffer)
        self.digest.update(memoryview(buffer)[:size])
        return size


@contextmanager
def _source(path, expected):
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("Ensembl sources require lowercase SHA-256 byte checksums")
    with open(path, "rb") as stream:
        reader = _HashedReader(stream)
        with io.BufferedReader(reader) as buffered:
            compressed = buffered.peek(2).startswith(b"\x1f\x8b")
            decoded = gzip.GzipFile(fileobj=buffered) if compressed else buffered
            with io.TextIOWrapper(decoded, encoding="utf8") as text:
                yield text
                # All parsers consume the complete file. Hashing the same open
                # byte stream prevents a separate hash/parse file-swap race.
                if text.read(1):
                    raise ValueError("Ensembl source parser did not consume the complete file")
                if reader.digest.hexdigest() != expected:
                    raise ValueError(f"Ensembl source checksum mismatch: {path}")


def _versioned(attributes, kind):
    identifier, version = attributes.get(kind + "_id", ""), attributes.get(kind + "_version", "")
    letter = {"gene": "G", "transcript": "T", "protein": "P"}[kind]
    if not re.fullmatch(r"ENS[A-Z]*" + letter + r"[0-9]+", identifier) or not re.fullmatch(r"[1-9][0-9]*", version):
        raise ValueError(f"Ensembl GTF requires a versioned {kind} identity")
    return identifier + "." + version


@dataclass(frozen=True, slots=True)
class _Transcript:
    gene: str
    name: str | None
    biotype: str | None
    contig: str
    strand: str
    start: int
    end: int


def _annotation(path, expected, reference):
    transcripts, products = defaultdict(set), defaultdict(lambda: defaultdict(set))
    accession = None
    with _source(path, expected) as stream:
        for line in stream:
            if line.startswith("#!genome-build-accession "):
                value = line.split()[-1]
                if accession is not None and accession != value:
                    raise ValueError("Conflicting Ensembl assembly headers")
                accession = value
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed Ensembl GTF row")
            if fields[2] not in {"transcript", "CDS"}:
                continue
            if (not fields[0] or fields[6] not in {"+", "-"}
                    or any(not re.fullmatch(r"[1-9][0-9]*", fields[i]) for i in (3, 4))
                    or int(fields[3]) > int(fields[4])):
                raise ValueError("Invalid Ensembl transcript/CDS locus")
            attributes = {}
            for key, value in re.findall(r'(\w+) "([^"]*)"', fields[8]):
                if key != "tag" and key in attributes and attributes[key] != value:
                    raise ValueError(f"Conflicting Ensembl GTF attribute: {key}")
                attributes[key] = value
            tid, gene = _versioned(attributes, "transcript"), _versioned(attributes, "gene")
            start, end = int(fields[3]), int(fields[4])
            if fields[2] == "transcript":
                transcripts[tid].add(_Transcript(gene, attributes.get("gene_name"),
                    attributes.get("transcript_biotype"), fields[0], fields[6], start, end))
            else:
                pid = _versioned(attributes, "protein")
                products[pid][(tid, gene, fields[0], fields[6])].add((start, end))
    if accession != reference.assembly_accession:
        raise ValueError("Ensembl GTF assembly accession disagrees with reference")
    if not transcripts or not products:
        raise ValueError("Ensembl GTF lacks a translated transcript inventory")
    return transcripts, products


def _fasta(paths, hashes, *, protein=False):
    records, duplicate_count = {}, 0

    def retain(header, chunks):
        nonlocal duplicate_count
        if header is None:
            return
        parts = header.split()
        identifier = parts[0]
        kind = "P" if protein else "T"
        if not re.fullmatch(r"ENS[A-Z]*" + kind + r"[0-9]+\.[1-9][0-9]*", identifier):
            raise ValueError("Ensembl FASTA requires exact versioned accessions")
        attributes = {}
        for part in parts[1:]:
            if ":" not in part:
                continue
            key, value = part.split(":", 1)
            if key in {"gene", "transcript"} and key in attributes and attributes[key] != value:
                raise ValueError("Conflicting Ensembl FASTA header identity")
            attributes[key] = value
        sequence = "".join(chunks)
        if not sequence or not re.fullmatch(r"[A-Z*]+" if protein else r"[A-Z]+", sequence):
            raise ValueError(f"Invalid Ensembl sequence: {identifier}")
        gene, transcript = attributes.get("gene"), attributes.get("transcript")
        if not gene or (protein and not transcript):
            raise ValueError("Ensembl FASTA lacks gene/transcript header identities")
        if not protein and transcript not in (None, identifier):
            raise ValueError("Ensembl RNA header transcript disagrees with its accession")
        row = (gene, transcript, len(sequence), hashlib.sha256(sequence.encode()).hexdigest(),
               sequence if protein else None)
        if identifier in records:
            if records[identifier] != row:
                raise ValueError(f"Conflicting Ensembl FASTA identity: {identifier}")
            duplicate_count += 1
        records[identifier] = row

    for path, expected in zip(paths, hashes):
        with _source(path, expected) as stream:
            header, chunks = None, []
            for line in stream:
                if line.startswith(">"):
                    retain(header, chunks)
                    header, chunks = line[1:].strip(), []
                    if not header:
                        raise ValueError("Empty Ensembl FASTA header")
                elif line.strip():
                    if header is None:
                        raise ValueError("Ensembl FASTA sequence precedes its header")
                    chunks.append(line.strip())
            retain(header, chunks)
    return records, duplicate_count


def ensembl_proteome(cls, protein_fasta, gtf, rna_fasta, ncrna_fasta, *, reference, expected_sha256):
    if (not isinstance(reference, ReferenceContext) or reference.annotation_name != "Ensembl"
            or not reference.annotation_version or not reference.assembly_accession):
        raise ValueError("Ensembl requires explicit taxon, assembly accession and Ensembl annotation")
    roles = {"protein", "gtf", "rna", "ncrna"}
    if not isinstance(expected_sha256, dict) or set(expected_sha256) != roles:
        raise ValueError("Ensembl expected_sha256 requires protein, gtf, rna and ncrna")
    transcripts, products = _annotation(gtf, expected_sha256["gtf"], reference)
    proteins, protein_duplicates = _fasta([protein_fasta], [expected_sha256["protein"]], protein=True)
    rna, rna_duplicates = _fasta([rna_fasta, ncrna_fasta],
                                [expected_sha256["rna"], expected_sha256["ncrna"]])
    if not proteins or set(proteins) != set(products):
        raise ValueError("Ensembl protein FASTA and GTF protein inventories disagree")
    for tid, row in rna.items():
        if tid not in transcripts or {t.gene for t in transcripts[tid]} != {row[0]}:
            raise ValueError(f"Ensembl RNA/GTF gene identity conflict: {tid}")
    source_hash = _sha256(expected_sha256)
    sequence_hash = _sha256(sorted((pid, row[4]) for pid, row in proteins.items()))
    if (reference.source_sha256 not in (None, source_hash)
            or reference.sequence_sha256 not in (None, sequence_hash)):
        raise ValueError("Ensembl source/sequence checksums disagree with reference")
    reference = replace(reference, source_sha256=source_hash, sequence_sha256=sequence_hash)
    entries, spans = [], {}
    source = f"Ensembl:{reference.annotation_version}"
    for pid, (gene, tid, _, _, sequence) in sorted(proteins.items()):
        if {(p[0], p[1]) for p in products[pid]} != {(tid, gene)}:
            raise ValueError(f"Ensembl protein/GTF identity conflict: {pid}")
        if tid not in rna or rna[tid][0] != gene:
            raise ValueError(f"Ensembl protein lacks its exact RNA identity: {pid}")
        placements = []
        for (_, _, contig, strand), blocks in products[pid].items():
            loci = [t for t in transcripts[tid] if (t.gene, t.contig, t.strand) == (gene, contig, strand)]
            if any(not any(t.start <= a <= b <= t.end for t in loci) for a, b in blocks):
                raise ValueError(f"Ensembl CDS has no containing transcript locus: {pid}")
            placements.extend(t for t in loci if any(t.start <= a <= b <= t.end for a, b in blocks))
        placements = sorted(set(placements), key=lambda t: (t.contig, t.strand, t.start, t.end))
        if not placements:
            raise ValueError(f"Ensembl protein lacks transcript occurrences: {pid}")
        identities = [(t.gene, t.contig, t.strand, t.start, t.end) for t in placements]
        if len(set(identities)) != len(identities):
            raise ValueError(f"Conflicting Ensembl transcript metadata: {tid}")
        occurrences = tuple(Occurrence(tid, t.gene.rsplit(".", 1)[0], t.name, t.biotype,
            t.contig, t.strand, None, None, (), True) for t in placements)
        entries.append(ProteinSequence(pid, sequence, "annotated", reference.taxon_id, source,
            occurrences=occurrences, source_transcript_id=tid,
            source_binding="gtf-rna-peptide-accession-join; spliced CDS coordinates not reconstructed"))
        spans[pid] = [[t.start, t.end] for t in placements]
    layer = dict(kind="annotated", role="host", source=source, reference=reference.key,
        added=len(entries), content_provenance=dict(expected_sha256),
        identifier_namespaces={"gene": "Ensembl gene", "transcript": "Ensembl transcript", "protein": "Ensembl protein"},
        mapping="gtf-rna-peptide-accession-join", transcript_spans=spans,
        annotated_transcripts=len(transcripts), rna_records=len(rna),
        untranslated_transcripts_without_rna=sorted(set(transcripts) - set(rna)),
        identical_protein_fasta_duplicates=protein_duplicates, identical_rna_fasta_duplicates=rna_duplicates)
    return cls(reference, entries, (layer,))
