"""Rebuild the rCRS mitochondrion annotation that the t2t-chm13 genome carries as contig NC_012920.1.

Ensembl's rapid-release annotation of T2T-CHM13v2.0 has no mitochondrial genes,
and the T2T Consortium's ``chm13v2.0_maskedY_rCRS`` analysis set, the one its
README recommends for short-read variant calling, replaces CHM13's own
mitochondrion (CP068254.1) with the rCRS (NC_012920.1). The rCRS is GRCh38's
MT, so Ensembl's GRCh38 annotation of MT describes that file's chrM exactly.

This copies Ensembl GRCh38 release 115's MT records (37 genes: 13 protein
coding, 22 tRNA and 2 rRNA) and names the contig ``NC_012920.1``, so they can
never be read as CHM13's own mitochondrion, which Ensembl would call ``MT``.
Before writing, it checks:

- every source file against its pinned SHA-256;
- Ensembl's MT DNA against the rCRS sequence;
- every transcript's spliced exons against its cDNA or ncRNA record;
- every CDS, translated with the vertebrate mitochondrial code (table 2), against
  its protein record.

Run it from a checkout; without ``--source`` it downloads the five Ensembl
files (about 270 MB):

    python tools/build_rcrs_mitochondrion.py --source gtf=Homo_sapiens.GRCh38.115.gtf.gz ...
"""

import argparse
import gzip
import hashlib
import re
import tempfile
import urllib.request
from pathlib import Path

from orfdata.sequence import genetic_code
from orfdata.viral_annotation import extract, translate

RELEASE = 115
ENSEMBL = f"https://ftp.ensembl.org/pub/release-{RELEASE}"
# SHA-256 of each file as published; each also matches the BSD sum in Ensembl's CHECKSUMS.
SOURCES = {
    "gtf": (f"{ENSEMBL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{RELEASE}.gtf.gz",
            "2f8e31578c3aa2f35646927c4a3b3b0dcf0321e57c0ebd3ecc81afcbc836d1a8"),
    "cdna": (f"{ENSEMBL}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
             "93af088ed22a46daa9645fc6292fb120c1e4ffe791d24abbcbd35062666c03f1"),
    "ncrna": (f"{ENSEMBL}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
              "7825db09257902df82172928537c1da16cfe79edcca56425556f74f1507437a0"),
    "pep": (f"{ENSEMBL}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
            "7991ab37d09986eefdb5c77032b15370ab0e99918b78f235f3230c0302851c96"),
    "dna": (f"{ENSEMBL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
            "f446a3d571e05d3ded44bf1d6c630fdced5e584c6651cf31f62fd4e13cbd1574"),
}
CONTIG = "NC_012920.1"
# SHA-256 of NC_012920.1's sequence, upper case with no line breaks (NCBI EFetch, 2026-11-09).
RCRS_SHA256 = "f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66"
OUTPUT = Path(__file__).resolve().parent.parent / "orfdata/data/genomes"
PREFIX = f"rcrs-ensembl{RELEASE}"
TABLE = 2  # Vertebrate mitochondrial


def read_source(role, path, directory):
    url, sha256 = SOURCES[role]
    if path is None:
        path = Path(directory) / url.rsplit("/", 1)[-1]
        urllib.request.urlretrieve(url, path)
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha256:
        raise ValueError(f"{path} is not Ensembl's {url.rsplit('/', 1)[-1]} (SHA-256 mismatch)")
    return gzip.decompress(raw).decode()


def fasta_records(text, keep):
    """The records whose header ``keep`` accepts, as (header, sequence) pairs in file order."""
    found, header, lines = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None and keep(header):
                found.append((header, "".join(lines)))
            header, lines = line[1:], []
        else:
            lines.append(line.strip())
    if header is not None and keep(header):
        found.append((header, "".join(lines)))
    return found


def on_mt(header):
    return ":GRCh38:MT:" in header


def attributes(field):
    return dict(re.findall(r'(\S+) "([^"]*)"', field))


def build(texts):
    dna = "".join(sequence for _, sequence in fasta_records(texts["dna"], lambda h: True)).upper()
    if hashlib.sha256(dna.encode()).hexdigest() != RCRS_SHA256:
        raise ValueError("Ensembl's GRCh38 MT is not the rCRS (NC_012920.1)")
    rows = [line.split("\t") for line in texts["gtf"].splitlines()
            if not line.startswith("#") and line.split("\t", 1)[0] == "MT"]
    exons, cds, stops = {}, {}, {}
    for row in rows:
        transcript = attributes(row[8]).get("transcript_id")
        span = (int(row[3]), int(row[4]), row[6])
        if row[2] == "exon":
            exons.setdefault(transcript, []).append((int(attributes(row[8])["exon_number"]), span))
        elif row[2] == "CDS":
            cds.setdefault(transcript, []).append(span)
        elif row[2] == "stop_codon":
            stops.setdefault(transcript, []).append(span)
    transcripts = fasta_records(texts["cdna"], on_mt) + fasta_records(texts["ncrna"], on_mt)
    proteins = fasta_records(texts["pep"], on_mt)
    by_id = {header.split()[0].split(".")[0]: sequence for header, sequence in transcripts}
    if set(by_id) != set(exons) or len(exons) != 37:
        raise ValueError(f"MT transcripts disagree: {len(exons)} in the GTF, {len(by_id)} in the FASTAs")
    for transcript, numbered in exons.items():
        spliced = extract(dna, [span for _, span in sorted(numbered)])
        if spliced != by_id[transcript]:
            raise ValueError(f"{transcript}'s exons do not spell its cDNA")
    code = genetic_code(TABLE)
    peptides = {re.search(r"transcript:(\S+)", header)[1].split(".")[0]: sequence for header, sequence in proteins}
    if set(peptides) != set(cds) or len(cds) != 13:
        raise ValueError(f"MT proteins disagree: {len(cds)} CDSs in the GTF, {len(peptides)} proteins")
    for transcript, spans in cds.items():
        # Some MT stops are completed by polyadenylation, so the GTF's CDS is the coding sequence without them.
        coding = extract(dna, spans)
        protein = translate(coding, table=TABLE)
        if coding[:3] in code.start_codons:
            protein = "M" + protein[1:]
        if protein != peptides[transcript]:
            raise ValueError(f"{transcript}'s CDS does not translate to its protein with table {TABLE}")
    header = (f"#!orfdata: Ensembl GRCh38 release {RELEASE}'s MT annotation, its contig renamed {CONTIG} "
              f"(the rCRS); from {SOURCES['gtf'][0]}\n")
    gtf = header + "".join("\t".join([CONTIG, *row[1:]]) + "\n" for row in rows)
    return dict(gtf=gtf, transcripts=fasta_text(transcripts), proteins=fasta_text(proteins))


def fasta_text(records):
    return "".join(f">{header}\n" + "".join(sequence[i:i + 60] + "\n" for i in range(0, len(sequence), 60))
                   for header, sequence in records)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source", action="append", default=[], metavar="ROLE=PATH",
                        help=f"A local copy of a source file instead of downloading it; roles: {', '.join(SOURCES)}")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    paths = dict(item.split("=", 1) for item in args.source)
    with tempfile.TemporaryDirectory() as directory:
        texts = {role: read_source(role, paths.get(role), directory) for role in SOURCES}
    built = build(texts)
    names = dict(gtf=f"{PREFIX}.gtf.gz", transcripts=f"{PREFIX}.cdna.fa.gz", proteins=f"{PREFIX}.pep.fa.gz")
    args.output.mkdir(parents=True, exist_ok=True)
    for role, text in built.items():
        raw = gzip.compress(text.encode(), mtime=0)
        (args.output / names[role]).write_bytes(raw)
        print(f'"{role}": {{"resource": "data/genomes/{names[role]}", '
              f'"sha256": "{hashlib.sha256(raw).hexdigest()}", "size": {len(raw)}}}')


if __name__ == "__main__":
    main()
