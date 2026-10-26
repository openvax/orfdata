"""Build the dog search database: annotated proteins, retention-time standards and candidate
ncORF translations.

usage: python -I build_db.py --gtf GTF --cdna CDNA_FA --ncrna NCRNA_FA --pep PEP_FA
                             --peptideatlas-fasta PA_FASTA --pyensembl-cache DIR
                             --out-dir DIR [--min-length 10]

Canonical class (FASTA ids CANON000000...):
  Ensembl 110 ROS_Cfam_1.0 pep.all, then the PeptideAtlas build FASTA's targets (UniProt
  UP000805418, Ensembl and 499 contaminants; DECOY_ entries dropped), then one entry per
  retention-time standard peptide (standards.py). Identical sequences are merged.
ncORF class (ids NCORF000000..., ordered by contig and start): ATG-initiated ORFs of at least
--min-length amino acids (stop excluded) that end at a stop inside the transcript:
  uORF    starts in the 5' UTR and ends before the CDS start
  uoORF   starts in the 5' UTR out of frame and ends inside or after the CDS
  dORF    starts in the 3' UTR
  lncORF  any frame of a transcript of a lncRNA gene
  ncTxORF any frame of a retained_intron or processed_transcript transcript
  For each stop codon only the most upstream ATG is kept (nested in-frame ATGs are suffixes),
  and in-frame upstream ATGs with no stop before the CDS (N-terminal extensions) are skipped.
  Candidates are deduplicated by sequence and dropped if identical to or contained in a
  canonical sequence, or if they contain a standard peptide that no canonical protein contains
  (I and L treated as equal).

Writes OUT_DIR/dog_db.fasta, OUT_DIR/ncorf_candidates.tsv and OUT_DIR/build_db_summary.json.
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from standards import STANDARDS, contained_standards  # noqa: E402

CODON = {}
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _a in "TCAG":
    for _b in "TCAG":
        for _c in "TCAG":
            CODON[_a + _b + _c] = _AAS[_i]
            _i += 1
STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_canonical(pep, pa_fasta):
    """Canonical sequences -> ids, in insertion order (Ensembl, PeptideAtlas, standards), and
    the standards that no canonical protein contains."""
    canon = {}
    n = collections.Counter()
    for name, s in read_fasta(pep):
        canon.setdefault(s.replace("*", ""), []).append("ens|" + name.split()[0])
        n["ensembl"] += 1
    for name, s in read_fasta(pa_fasta):
        if name.startswith("DECOY_"):
            continue
        pid = name.split()[0]
        tag = "contam|" if pid.startswith("CONTAM_") else "pa|"
        n["contaminant" if tag == "contam|" else "peptideatlas"] += 1
        canon.setdefault(s, []).append(tag + pid)
    synthetic = add_standards(canon)
    n["standard"] = len(STANDARDS)
    return canon, n, synthetic


def add_standards(canon):
    """Add one canonical entry per standard peptide, so its spectra never go to an ncORF.
    Returns the standards that no sequence already in canon contains (I and L equal)."""
    natural = set(contained_standards("\n".join(canon)))
    for kit, peptide in STANDARDS:
        canon.setdefault(peptide, []).append(f"std|{kit}|{peptide}")
    return {p for _, p in STANDARDS if p not in natural}


def translate_to_stop(seq, p):
    """Translate from p; return (aa, exclusive end offset including the stop), or (None, None)."""
    aa = []
    for j in range(p, len(seq) - 2, 3):
        c = seq[j:j + 3]
        if c in STOPS:
            return "".join(aa), j + 3
        aa.append(CODON.get(c, "X"))
    return None, None


def spliced_to_genomic(t, off):
    """0-based spliced offset -> 1-based genomic position."""
    exons = sorted(t.exons, key=lambda e: e.start, reverse=(t.strand == "-"))
    for e in exons:
        n = e.end - e.start + 1
        if off < n:
            return e.start + off if t.strand == "+" else e.end - off
        off -= n
    return None


def scan_transcripts(genome, min_length):
    cands = collections.defaultdict(list)  # aa sequence -> origins
    stats = collections.Counter()
    mism = 0
    for t in genome.transcripts():
        try:
            seq = t.sequence
        except Exception:
            seq = None
        if not seq:
            stats["no_seq_" + t.biotype] += 1
            continue
        seq = seq.upper()
        gene_bt = t.gene.biotype
        regions = []  # (class, scan start, scan end, CDS start)
        if t.biotype in ("protein_coding", "nonsense_mediated_decay") and t.complete:
            try:
                cs = t.first_start_codon_spliced_offset
                ce = t.last_stop_codon_spliced_offset + 1  # Exclusive, includes the stop.
            except Exception:
                stats["bad_cds"] += 1
                continue
            if t.protein_sequence and t.protein_sequence.rstrip("*") != \
                    "".join(CODON.get(seq[j:j + 3], "X") for j in range(cs, ce - 3, 3)):
                mism += 1
            stats["coding_tx_" + t.biotype] += 1
            regions.append(("5UTR", 0, cs, cs))
            regions.append(("3UTR", ce, len(seq), cs))
        elif t.biotype == "lncRNA" or gene_bt == "lncRNA":
            stats["lnc_tx"] += 1
            regions.append(("lncORF", 0, len(seq), None))
        elif t.biotype in ("retained_intron", "processed_transcript"):
            stats["nctx_" + t.biotype] += 1
            regions.append(("ncTxORF", 0, len(seq), None))
        else:
            stats["skipped_" + t.biotype] += 1
            continue
        for cls, a, b, cs in regions:
            seen_stops = set()
            for p in range(a, b - 2):
                if seq[p:p + 3] != "ATG":
                    continue
                aa, end = translate_to_stop(seq, p)
                if aa is None:
                    stats[cls + "_nostop"] += 1
                    continue
                if end in seen_stops:  # Nested in-frame ATG.
                    continue
                seen_stops.add(end)
                if len(aa) < min_length:
                    continue
                kind = cls
                if cls == "5UTR":
                    if end <= cs:
                        kind = "uORF"
                    elif (cs - p) % 3 == 0:
                        stats["Nterm_extension_skipped"] += 1
                        continue
                    else:
                        kind = "uoORF"
                elif cls == "3UTR":
                    kind = "dORF"
                cands[aa].append(dict(kind=kind, tx=t.transcript_id, gene=t.gene_id, gene_name=t.gene_name,
                                      tx_biotype=t.biotype, contig=t.contig, strand=t.strand,
                                      start=spliced_to_genomic(t, p), stop_end=spliced_to_genomic(t, end - 1),
                                      tx_off=p))
    return cands, stats, mism


K = 10
_AA = "ACDEFGHIKLMNPQRSTVWYUOBZXJ"
_ENC = np.full(256, 31, dtype=np.int64)
for _k, _ch in enumerate(_AA):
    _ENC[ord(_ch)] = _k


def kmers(s):
    a = _ENC[np.frombuffer(s.encode(), dtype=np.uint8)]
    if len(a) < K:
        return np.zeros(0, dtype=np.int64)
    v = np.zeros(len(a) - K + 1, dtype=np.int64)
    for j in range(K):
        v = (v << 5) | a[j:len(a) - K + 1 + j]
    return v


def filter_candidates(cands, canon, synthetic):
    """Drop candidates identical to or contained in a canonical sequence, or containing a
    synthetic standard (one no canonical protein contains). A candidate that shares a standard's
    sequence with a dog protein stays: that peptide is filtered later like any shared one."""
    allk = np.unique(np.concatenate([kmers(s) for s in canon]))
    big = "\n".join(canon)
    removed = collections.Counter()
    keep = {}
    for aa, origins in cands.items():
        if aa in canon:
            removed["identical"] += 1
            continue
        km = kmers(aa)
        idx = np.searchsorted(allk, km)
        idx[idx >= len(allk)] = 0
        if np.all(allk[idx] == km) and big.find(aa) >= 0:
            removed["contained"] += 1
            continue
        if synthetic.intersection(contained_standards(aa)):
            removed["contains_standard"] += 1
            continue
        keep[aa] = origins
    return keep, removed


def write_outputs(out_dir, canon, keep):
    order = sorted(keep, key=lambda s: (keep[s][0]["contig"], keep[s][0]["start"] or 0))
    kinds = collections.Counter()
    nres = 0
    with open(os.path.join(out_dir, "dog_db.fasta"), "w") as fa, \
            open(os.path.join(out_dir, "ncorf_candidates.tsv"), "w") as tsv:
        tsv.write("orf_id\tlength\tkinds\tgenes\tgene_names\ttranscripts\tlocus\tsequence\n")
        for cid, s in enumerate(canon):
            fa.write(f">CANON{cid:06d} {' '.join(canon[s][:5])}\n{s}\n")
        for n, s in enumerate(order):
            o = keep[s]
            oid = f"NCORF{n:06d}"
            ks = sorted({x["kind"] for x in o})
            kinds.update(ks)
            nres += len(s)
            loci = sorted({f"{x['contig']}:{x['start']}-{x['stop_end']}:{x['strand']}" for x in o})
            fa.write(f">{oid} {','.join(ks)} {o[0]['gene_name'] or o[0]['gene']}\n{s}\n")
            tsv.write("\t".join([oid, str(len(s)), ",".join(ks), ",".join(sorted({x['gene'] for x in o})),
                                 ",".join(sorted({x['gene_name'] or '' for x in o})),
                                 ",".join(sorted({x['tx'] for x in o})), ";".join(loci), s]) + "\n")
    return kinds, nres, len(order)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="Ensembl 110 dog GTF (.gtf.gz)")
    ap.add_argument("--cdna", required=True, help="Ensembl 110 dog cdna.all.fa.gz")
    ap.add_argument("--ncrna", required=True, help="Ensembl 110 dog ncrna.fa.gz")
    ap.add_argument("--pep", required=True, help="Ensembl 110 dog pep.all.fa.gz")
    ap.add_argument("--peptideatlas-fasta", required=True,
                    help="UP000805418_2023-09-25_ENS-ROS_Cfam_1.0_CONTAM_targetdecoy.fasta from PXD052932")
    ap.add_argument("--pyensembl-cache", required=True, help="directory for PyEnsembl's index of the GTF")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-length", type=int, default=10, help="minimum ncORF length in aa, stop excluded")
    a = ap.parse_args()

    import pyensembl  # Imported here so the module's helpers work without it.

    os.makedirs(a.out_dir, exist_ok=True)
    canon, n, synthetic = load_canonical(a.pep, a.peptideatlas_fasta)
    log = {"inputs": {k: {"path": os.path.abspath(p), "sha256": sha256(p)} for k, p in (
        ("gtf", a.gtf), ("cdna", a.cdna), ("ncrna", a.ncrna), ("pep", a.pep),
        ("peptideatlas_fasta", a.peptideatlas_fasta))}, "min_length": a.min_length}
    print(f"canonical: {n['ensembl']} Ensembl pep, {n['peptideatlas']} PeptideAtlas/UniProt targets, "
          f"{n['contaminant']} contaminants, {n['standard']} standard peptides -> {len(canon)} unique sequences",
          file=sys.stderr)
    natural = sorted(p for _, p in STANDARDS if p not in synthetic)
    print(f"standards that canonical proteins already contain: {natural}", file=sys.stderr)
    log["standards_in_canonical_proteins"] = natural

    genome = pyensembl.Genome(reference_name="ROS_Cfam_1.0", annotation_name="ensembl", annotation_version=110,
                              gtf_path_or_url=a.gtf, transcript_fasta_paths_or_urls=[a.cdna, a.ncrna],
                              protein_fasta_paths_or_urls=[a.pep], cache_directory_path=a.pyensembl_cache)
    genome.index()
    cands, stats, mism = scan_transcripts(genome, a.min_length)
    print(f"transcript stats: {dict(stats)}; CDS translation mismatches vs pep: {mism}", file=sys.stderr)
    print(f"unique candidate ORF sequences before canonical filter: {len(cands)}", file=sys.stderr)
    keep, removed = filter_candidates(cands, canon, synthetic)
    print(f"removed: {dict(removed)}; kept {len(keep)}", file=sys.stderr)
    kinds, nres, n_orfs = write_outputs(a.out_dir, canon, keep)
    canon_res = sum(len(s) for s in canon)
    print(f"ncORF kinds (an ORF can carry several): {dict(kinds)}", file=sys.stderr)
    print(f"FASTA: {len(canon)} canonical ({canon_res} aa) + {n_orfs} ncORF ({nres} aa)", file=sys.stderr)
    log.update(canonical_inputs=dict(n), canonical_sequences=len(canon), canonical_residues=canon_res,
               transcript_stats=dict(stats), cds_translation_mismatches=mism,
               unique_candidates=len(cands), removed=dict(removed), ncorf_sequences=n_orfs,
               ncorf_residues=nres, ncorf_kinds=dict(kinds),
               dog_db_fasta_sha256=sha256(os.path.join(a.out_dir, "dog_db.fasta")))
    with open(os.path.join(a.out_dir, "build_db_summary.json"), "w") as f:
        json.dump(log, f, indent=1)


if __name__ == "__main__":
    main()
