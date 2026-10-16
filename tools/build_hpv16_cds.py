"""Rebuild the small HPV16 CDS catalogue from its versioned RefSeq GenBank record.

Requires Biopython only for building, not for loading the released catalogue.
The GDC placement is an exact circular rotation, verified against the DNA MD5
of the pinned GDC reference, not an alignment-based annotation transfer.
"""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path

from Bio import SeqIO

ACCESSION = "NC_001526.4"
NATIVE_MD5 = "08d2b3849496d12cc45fd4df713d5c2f"
GDC_MD5 = "f2463721818c29998ef5ab917e6df74a"
ROTATION = 864
SOURCE = ("https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"
          "id=NC_001526.4&db=nuccore&report=genbank&retmode=text")
HOST_SOURCE = ("PaVE HPV16REF /host=Human (Homo sapiens), checked 2026-10-16: "
               "https://pave.niaid.nih.gov/locus_viewer?feature_id=HPV16_E7&seq_id=HPV16REF")
TAXONOMY_SOURCE = "https://ictv.global/report/chapter/papillomaviridae/taxonomy/papillomaviridae"


def build(genbank):
    record = SeqIO.read(io.StringIO(genbank.decode()), "genbank")
    native = str(record.seq).upper()
    gdc = native[-ROTATION:] + native[:-ROTATION]
    if record.id != ACCESSION or hashlib.md5(native.encode()).hexdigest() != NATIVE_MD5:
        raise ValueError("The source is not the pinned NC_001526.4 DNA")
    if hashlib.md5(gdc.encode()).hexdigest() != GDC_MD5:
        raise ValueError("The circular rotation does not match the GDC DNA")
    definitions = []
    native_blocks, gdc_blocks = {}, {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        q = feature.qualifiers
        if feature.location.strand != 1 or int(q.get("codon_start", [1])[0]) != 1:
            raise ValueError("This source builder requires complete plus-strand CDSs")
        protein_id = q["protein_id"][0]
        blocks = [[int(part.start) + 1, int(part.end)] for part in feature.location.parts]
        moved = [[(a - 1 + ROTATION) % len(native) + 1,
                  (b - 1 + ROTATION) % len(native) + 1] for a, b in blocks]
        if any(a > b for a, b in moved):
            raise ValueError("A rotated coding block crosses the origin; qualify it explicitly")
        coding = feature.extract(record.seq)
        protein = q["translation"][0]
        table = int(q.get("transl_table", [1])[0])
        if str(coding.translate(table=table, cds=True)) != protein:
            raise ValueError(f"Source translation mismatch: {protein_id}")
        definitions.append(dict(protein_id=protein_id, gene=q.get("gene", [None])[0],
                                product=q["product"][0], protein_sequence=protein,
                                translation_table=table, strand="+", phase=0,
                                start_complete=True, stop_complete=True))
        native_blocks[protein_id] = blocks
        gdc_blocks[protein_id] = moved
    if len(definitions) != 9:
        raise ValueError("Expected the nine annotated CDSs in this source version")
    source = dict(accession=ACCESSION, url=SOURCE, sha256=hashlib.sha256(genbank).hexdigest())
    references = []
    for name, sequence, blocks, offset in (
            (ACCESSION, native, native_blocks, 0),
            ("GRCh38.d1.vd1:HPV16", gdc, gdc_blocks, ROTATION)):
        references.append(dict(id=name, taxon_id=333760, species="human_papillomavirus_16",
            length=len(sequence), sequence_md5=hashlib.md5(sequence.encode()).hexdigest(),
            sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            source=source, native_to_reference_rotation=offset, coding_blocks=blocks,
            reference_sequence=sequence, host_taxon_ids=[9606], host_source=HOST_SOURCE,
            virus=dict(group="HPV", type="HPV16", ictv_species="Alphapapillomavirus 9",
                       taxonomy_source=TAXONOMY_SOURCE)))
    return dict(schema="orfdata.viral_cds.v2", source=source,
                scope="Nine source-annotated HPV16 CDSs; not all HPV ORFs or complete transcript models",
                references=references, cds=definitions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genbank", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raw = args.genbank.read_bytes()
    if args.genbank.suffix == ".gz":
        raw = gzip.decompress(raw)
    payload = json.dumps(build(raw), sort_keys=True, separators=(",", ":")).encode() + b"\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress(payload, mtime=0))
    print(hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
