"""EXPERIMENTAL. Aggregate per-run Sage rank-1 PSMs (parquet) into class-specific FDR, filtered
ncORF peptides and a machine-readable table of supported ORFs. Can run at any point of a
search; it uses whichever runs have finished.

usage: python -I aggregate.py CONFIG_JSON

CONFIG_JSON keys (README.md describes each one and every output):
  required: psm_dir, runs_tsv, out_dir, db_fasta, candidates_tsv, human_pep, mouse_pep,
            gtf, cdna, ncrna, pep, pyensembl_cache
  optional: manifest_tsv, extra_dog_fastas ([]), fdr (0.01), homology_min_score (100),
            homology_min_identity (0.4), isobaric_tolerance_da (0.02)
"""
import collections
import glob
import gzip
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdr import qvalues  # noqa: E402
from homology import blosum62_sequence  # noqa: E402
from standards import STANDARDS  # noqa: E402

with open(sys.argv[1]) as _f:
    cfg = json.load(_f)
FDR = cfg.get("fdr", 0.01)
HOMOLOGY_MIN_SCORE = cfg.get("homology_min_score", 100)
HOMOLOGY_MIN_IDENTITY = cfg.get("homology_min_identity", 0.4)
ISOBARIC_TOLERANCE = cfg.get("isobaric_tolerance_da", 0.02)
out = cfg["out_dir"]
os.makedirs(out, exist_ok=True)
T0 = time.time()
summary = []
STANDARDS_IL = {p.replace("I", "L") for _, p in STANDARDS}


def log(s):
    print(s, flush=True)
    summary.append(s)


STRIP = r"\[[^\]]*\]|[^A-Z]"


def tissue_of(run):
    if re.search(r"urine", run, re.I):
        return "urine"
    if re.search(r"plasma", run, re.I):
        return "plasma"
    m = re.match(r"\d{8}_(.+?)_WB", run) or re.match(r"\d{8}-Canine-(.+?)-\d\d-DDA", run)
    t = m.group(1).lower() if m else "unknown"
    return "adipose" if t.startswith("adipose") else t


# ---------------- load PSMs ----------------
runs = pd.read_csv(cfg["runs_tsv"], sep="\t", header=None,
                   names=["run", "group", "zip_url", "member", "mzml_bytes", "compressed_bytes"])
files = sorted(glob.glob(os.path.join(cfg["psm_dir"], "*.parquet")))
if not files:
    raise SystemExit(f"no per-run PSM parquet files found in {cfg['psm_dir']}")
NC_COLS = ["run", "filename", "scannr", "peptide", "proteins", "num_proteins", "label", "charge", "expmass", "calcmass",
           "precursor_ppm", "isotope_error", "hyperscore", "delta_next", "matched_peaks", "rt", "ion_mobility",
           "sage_discriminant_score", "spectrum_q", "peptide_q", "posterior_error"]
canon_parts, nc_parts = [], []
for fn in files:
    d = pd.read_parquet(fn, columns=NC_COLS)
    is_c = d.proteins.str.contains("CANON", regex=False)
    c = d.loc[is_c, ["run", "label", "sage_discriminant_score", "peptide"]].copy()
    c["stripped"] = c.peptide.str.replace(STRIP, "", regex=True).astype("string[pyarrow]")
    c = c.drop(columns="peptide")
    c["label"] = c.label.astype("int8")
    c["sage_discriminant_score"] = c.sage_discriminant_score.astype("float32")
    canon_parts.append(c)
    nc_parts.append(d.loc[~is_c].copy())
canon = pd.concat(canon_parts, ignore_index=True)
canon["run"] = canon.run.astype("category")
nc = pd.concat(nc_parts, ignore_index=True)
nc["stripped"] = nc.peptide.str.replace(STRIP, "", regex=True)
del canon_parts, nc_parts
log(f"runs aggregated: {len(files)} of {len(runs)} listed; rank-1 PSMs: canonical class {len(canon)}, ncORF class {len(nc)}")
if cfg.get("manifest_tsv") and os.path.exists(cfg["manifest_tsv"]):
    manifest = pd.read_csv(cfg["manifest_tsv"], sep="\t")
    missing = set(manifest.run) - {os.path.basename(f)[:-len(".parquet")] for f in files}
    log(f"run manifest: {len(manifest)} runs searched"
        + (f"; {len(missing)} of them have no PSM file in psm_dir: {sorted(missing)[:5]}" if missing else ""))

# ---------------- class-specific FDR ----------------
canon["decoy"] = canon.label == -1
nc["decoy"] = nc.label == -1
canon["q_class"] = qvalues(canon.sage_discriminant_score.values.astype(float), canon.decoy.values)
nc["q_class"] = qvalues(nc.sage_discriminant_score.values, nc.decoy.values)


def peptide_q(df):
    best = df.groupby(["stripped", "decoy"], observed=True).sage_discriminant_score.max().reset_index()
    best["q_pep"] = qvalues(best.sage_discriminant_score.values.astype(float), best.decoy.values)
    return best


cpep = peptide_q(canon)
npep = peptide_q(nc)
ct = canon[~canon.decoy]
log(f"canonical class: rank-1 targets {len(ct)}, decoys {int(canon.decoy.sum())}; PSMs at class q<={FDR}: "
    f"{int((ct.q_class <= FDR).sum())}; peptides at peptide-level class q<={FDR}: "
    f"{int(((cpep.q_pep <= FDR) & ~cpep.decoy).sum())}")
nt = nc[~nc.decoy]
log(f"ncORF class: rank-1 targets {len(nt)}, decoys {int(nc.decoy.sum())}; PSMs at class q<={FDR}: "
    f"{int((nt.q_class <= FDR).sum())}; target peptides at peptide-level class q<={FDR}: "
    f"{int(((npep.q_pep <= FDR) & ~npep.decoy).sum())} (decoy peptides above threshold: "
    f"{int((npep.decoy & (npep.sage_discriminant_score >= npep[(npep.q_pep <= FDR) & ~npep.decoy].sage_discriminant_score.min())).sum()) if ((npep.q_pep <= FDR) & ~npep.decoy).any() else 0})")
canon_pep_counts = {
    "canonical_psms_q01": int((ct.q_class <= FDR).sum()),
    "canonical_peptides_q01": int(((cpep.q_pep <= FDR) & ~cpep.decoy).sum()),
}
per_run_canon = ct[ct.q_class <= FDR].groupby("run", observed=True).size()

# ---------------- QC: retention-time standard PSMs (tabulated below) ----------------
ct_q = ct[ct.q_class <= FDR]
nt_q = nt[nt.q_class <= FDR]
std_psms = pd.concat([
    ct_q.loc[ct_q.stripped.str.replace("I", "L").isin(STANDARDS_IL), ["run", "stripped"]].assign(cls="canonical"),
    nt_q.loc[nt_q.stripped.str.replace("I", "L").isin(STANDARDS_IL), ["run", "stripped"]].assign(cls="ncorf"),
], ignore_index=True)
std_psms["run"] = std_psms.run.astype(str)
std_psms["stripped"] = std_psms.stripped.astype(str)
per_run_ncorf = nt_q.groupby("run").size()
del canon, ct, ct_q, nt_q

# ---------------- candidate ncORF peptides ----------------
passing_pep = set(npep[(npep.q_pep <= FDR) & ~npep.decoy].stripped)
nt = nt.merge(npep[~npep.decoy][["stripped", "q_pep"]], on="stripped", how="left")
sel = nt[nt.stripped.isin(passing_pep)].copy()
log(f"ncORF target peptides at peptide-level q<={FDR}: {len(passing_pep)} ({len(sel)} rank-1 PSMs, "
    f"{int((sel.q_class <= FDR).sum())} of them at PSM-level class q<={FDR})")

cand = pd.read_csv(cfg["candidates_tsv"], sep="\t").set_index("orf_id")


def read_fasta_seqs(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, cur, name = {}, [], None
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(cur)
                name, cur = line[1:].strip(), []
            else:
                cur.append(line.strip())
    if name is not None:
        seqs[name] = "".join(cur)
    return seqs


dbseqs = read_fasta_seqs(cfg["db_fasta"])
canon_seqs = [s for n, s in dbseqs.items() if n.startswith("CANON")]
dog_ref = list(canon_seqs)
extra_n = 0
for fa in cfg.get("extra_dog_fastas", []):
    ex = read_fasta_seqs(fa)
    extra_n += len(ex)
    dog_ref.extend(ex.values())
# Standards too, in case the database lacks them: a peptide near one is no ncORF evidence.
dog_ref.extend(p for _, p in STANDARDS)
dog_il = "\n".join(dog_ref).replace("I", "L")
log(f"dog reference for filtering: {len(canon_seqs)} canonical search sequences + {extra_n} extra dog sequences "
    f"+ {len(STANDARDS)} standard peptides")

# Standards whose sequence occurs in a dog protein can't be told from endogenous peptides (the
# PRTC kit's are heavy-labelled, which these Sage parameters don't search for), so only the
# others count as spike-in PSMs. A database built by build_db.py holds every standard, so their
# PSMs belong to the canonical class; ncORF-class standard PSMs mean the database lacks them.
natural_il = "\n".join(s for n, s in dbseqs.items() if n.startswith("CANON")
                       and not all(t.startswith("std|") for t in n.split()[1:])).replace("I", "L")
natural_standards = sorted(p for _, p in STANDARDS if p.replace("I", "L") in natural_il)
del natural_il
std_psms["natural"] = std_psms.stripped.str.replace("I", "L").isin({p.replace("I", "L") for p in natural_standards})
spike = std_psms[~std_psms.natural]
spike_canonical = spike[spike.cls == "canonical"]
qc = pd.DataFrame({
    "canonical_psms": per_run_canon,
    "ncorf_psms": per_run_ncorf,
    "standard_psms": spike_canonical.groupby("run").size(),
    "standard_peptides": spike_canonical.groupby("run").stripped.nunique(),
    "standard_psms_ncorf_class": spike[spike.cls == "ncorf"].groupby("run").size(),
    "natural_standard_psms": std_psms[std_psms.natural].groupby("run").size(),
}).reindex(sorted({os.path.basename(f)[:-len(".parquet")] for f in files})).fillna(0).astype(int)
qc.index.name = "run"
qc.insert(0, "tissue", [tissue_of(r) for r in qc.index])
qc.to_csv(os.path.join(out, "run_qc.tsv"), sep="\t")
log(f"retention-time standards at class q<={FDR}: {int(qc.standard_psms.sum())} canonical-class PSMs in "
    f"{int((qc.standard_psms > 0).sum())} of {len(qc)} runs ({spike_canonical.stripped.nunique()} peptides); "
    f"ncORF-class PSMs of standards: {int(qc.standard_psms_ncorf_class.sum())}; not counted because dog proteins "
    f"contain them: {natural_standards} ({int(qc.natural_standard_psms.sum())} PSMs)")

MASS = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=160.03065, L=113.08406,
            I=113.08406, N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891,
            F=147.06841, R=156.10111, Y=163.06333, W=186.07931, U=150.95364)
MODS = {"none": 0.0, "Met-ox": 15.99491, "deamid": 0.98402, "acetyl": 42.01057, "isotope": 1.00335}


def near_hits(pep, ref, maxd):
    """windows of ref (I->L) within Hamming maxd (>=1) of pep (I->L); pigeonhole on maxd+1 parts"""
    p = pep.replace("I", "L")
    L = len(p)
    cuts = [round(i * L / (maxd + 1)) for i in range(maxd + 2)]
    hits = {}
    for a, b in zip(cuts, cuts[1:]):
        for m in re.finditer("(?=" + re.escape(p[a:b]) + ")", ref):
            s = m.start() - a
            if s < 0 or s + L > len(ref):
                continue
            w = ref[s:s + L]
            if "\n" in w or w in hits:
                continue
            mm = [i for i in range(L) if w[i] != p[i]]
            if 1 <= len(mm) <= maxd:
                dm = sum(MASS.get(c, 0) for c in p) - sum(MASS.get(c, 0) for c in w)
                hits[w] = (mm, dm)
    return hits


def stop_keys(oid):
    keys = set()
    for loc in cand.at[oid, "locus"].split(";"):
        contig, rest = loc.split(":", 1)
        b = rest.rsplit(":", 1)[0].split("-")[1]
        keys.add((contig, b, loc[-1]))
    return keys


# reference strings for human / mouse
def ref_string(path):
    return "\n".join(read_fasta_seqs(path).values()).replace("I", "L")


human_il = ref_string(cfg["human_pep"])
mouse_il = ref_string(cfg["mouse_pep"])

rows = []
for pep, g in sel.groupby("stripped"):
    oids = sorted({p for ps in g.proteins for p in ps.split(";")})
    il = pep.replace("I", "L")
    in_dog = il in dog_il
    h2 = {} if in_dog else near_hits(pep, dog_il, 2)
    h1 = [w for w, (mm, dm) in h2.items() if len(mm) == 1]
    iso = [f"{w}(d={len(mm)},dm={dm:+.3f},{'/'.join(k for k, v in MODS.items() if abs(abs(dm) - v) < ISOBARIC_TOLERANCE)})"
           for w, (mm, dm) in h2.items()
           if len(mm) == 2 and any(abs(abs(dm) - v) < ISOBARIC_TOLERANCE for v in MODS.values())]
    common = set.intersection(*[stop_keys(o) for o in oids])
    hs_exact = il in human_il
    mm_exact = il in mouse_il
    hs1 = [] if hs_exact else [w for w, (mm, dm) in near_hits(pep, human_il, 1).items()]
    mm1 = [] if mm_exact else [w for w, (mm, dm) in near_hits(pep, mouse_il, 1).items()]
    gq = g[g.q_class <= FDR]
    per_run = gq.groupby("run").size().to_dict()
    best = g.sort_values("sage_discriminant_score", ascending=False).iloc[0]
    rows.append(dict(
        peptide=pep, q_peptide=float(g.q_pep.iloc[0]), best_q_psm=float(g.q_class.min()),
        n_psms_q01=len(gq), n_psms_rank1=len(g), n_runs=len(per_run),
        tissues=",".join(sorted({tissue_of(r) for r in per_run})),
        per_run_psms_q01=json.dumps(per_run, sort_keys=True),
        spectrum_ids=";".join(f"{f}:{s}" for f, s in zip(gq.filename, gq.scannr)),
        best_hyperscore=float(g.hyperscore.max()), best_score=float(g.sage_discriminant_score.max()),
        best_spectrum=f"{best.filename}:{best.scannr}", best_modified=best.peptide, best_charge=int(best.charge),
        best_precursor_ppm=round(float(best.precursor_ppm), 2), best_isotope_error=round(float(best.isotope_error), 2),
        best_delta_next=round(float(best.delta_next), 3), best_matched_peaks=int(best.matched_peaks),
        orf_ids=",".join(oids), orf_unique=bool(common),
        kinds=",".join(sorted({k for o in oids for k in cand.at[o, "kinds"].split(",")})),
        IL_in_dog_reference=in_dog, dog_1aa=";".join(h1[:5]), n_dog_1aa=len(h1),
        dog_near_isobaric_2aa=";".join(iso[:5]), n_dog_near_isobaric=len(iso),
        human_exact=hs_exact, human_1aa=";".join(hs1[:3]), mouse_exact=mm_exact, mouse_1aa=";".join(mm1[:3]),
    ))
P = pd.DataFrame(rows)
if len(P):
    P["fail_reason"] = np.select(
        [P.IL_in_dog_reference, P.n_dog_1aa > 0, P.n_dog_near_isobaric > 0, ~P.orf_unique],
        ["I/L-identical to dog reference", "1 aa from dog reference", "near-isobaric (2 aa) to dog reference",
         "not ORF-unique"], default="")
    P["passes"] = P.fail_reason == ""
    P["conserved_peptide"] = P.human_exact | P.mouse_exact | (P.human_1aa != "") | (P.mouse_1aa != "")
    P["classification"] = np.where(~P.passes, "filtered",
                                   np.where(P.conserved_peptide, "conserved_coding_probable_annotation_gap",
                                            "candidate_novel_ncORF"))
log("filter outcome for peptide-level ncORF peptides: " + json.dumps(P.fail_reason.replace("", "pass").value_counts().to_dict() if len(P) else {}))

# ---------------- ORF-level homology (secondary evidence) ----------------
from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

ok = P[P.passes] if len(P) else P
orf_ids = sorted({o for x in ok.orf_ids for o in x.split(",")}) if len(ok) else []
aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
K = 6
orf_hom = {}
for sp, path in (("human", cfg["human_pep"]), ("mouse", cfg["mouse_pep"])):
    seqs = read_fasta_seqs(path)
    names, sl = list(seqs), list(seqs.values())
    alignment_seqs = [blosum62_sequence(s) for s in sl]
    qk = collections.defaultdict(set)
    for o in orf_ids:
        s = blosum62_sequence(cand.at[o, "sequence"])
        for i in range(len(s) - K + 1):
            qk[s[i:i + K]].add(o)
    hits = collections.defaultdict(collections.Counter)
    for j, s in enumerate(alignment_seqs):
        seen = set()
        for i in range(len(s) - K + 1):
            k = s[i:i + K]
            if k in qk and k not in seen:
                seen.add(k)
                for o in qk[k]:
                    hits[o][j] += 1
    for o in orf_ids:
        q = blosum62_sequence(cand.at[o, "sequence"])
        bestaln = None
        for j, _ in hits[o].most_common(25):
            a = aligner.align(q, alignment_seqs[j])[0]
            if bestaln is None or a.score > bestaln[0]:
                bestaln = (a.score, j, a)
        if bestaln is None:
            orf_hom[(o, sp)] = dict(score=0.0, hit="", ident=0.0, aln_len=0)
            continue
        sc, j, a = bestaln
        pairs = [(qs + k, ts + k) for (qs, qe), (ts, te) in zip(*a.aligned) for k in range(qe - qs)]
        ident = sum(1 for qi, ti in pairs if q[qi] == alignment_seqs[j][ti]) / max(len(pairs), 1)
        gene = next((t.split(":", 1)[1] for t in names[j].split() if t.startswith("gene_symbol:")), "")
        orf_hom[(o, sp)] = dict(score=float(sc), hit=f"{names[j].split()[0]}({gene})", ident=round(ident, 3),
                                aln_len=len(pairs))


def orf_homologous(o):
    return any(orf_hom[(o, sp)]["score"] >= HOMOLOGY_MIN_SCORE and orf_hom[(o, sp)]["ident"] >= HOMOLOGY_MIN_IDENTITY
               for sp in ("human", "mouse"))


if len(P):
    P["orf_homologous_to_annotated_protein"] = [
        any(orf_homologous(o) for o in x.split(",")) if ps else None for x, ps in zip(P.orf_ids, P.passes)]
    P.to_csv(os.path.join(out, "ncorf_peptides_all_candidates.tsv"), sep="\t", index=False)
    P[P.passes].to_csv(os.path.join(out, "ncorf_peptides_passing.tsv"), sep="\t", index=False)

# ---------------- ORF table ----------------
import pyensembl  # noqa: E402

genome = pyensembl.Genome(reference_name="ROS_Cfam_1.0", annotation_name="ensembl", annotation_version=110,
                          gtf_path_or_url=cfg["gtf"], transcript_fasta_paths_or_urls=[cfg["cdna"], cfg["ncrna"]],
                          protein_fasta_paths_or_urls=[cfg["pep"]], cache_directory_path=cfg["pyensembl_cache"])
CODON = {}
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for x in "TCAG":
    for y in "TCAG":
        for z in "TCAG":
            CODON[x + y + z] = aas[i]
            i += 1


def placements(t, prot):
    seq = t.sequence.upper()
    res = []
    for p in range(0, len(seq) - 2):
        if seq[p:p + 3] != "ATG":
            continue
        aa = []
        for j in range(p, len(seq) - 2, 3):
            c = CODON.get(seq[j:j + 3], "X")
            if c == "*":
                if "".join(aa) == prot:
                    res.append((p, j + 3))
                break
            aa.append(c)
            if len(aa) > len(prot):
                break
    return res


def blocks_for(t, a, b):
    """spliced [a,b) -> ascending 1-based inclusive genomic blocks"""
    exons = sorted(t.exons, key=lambda e: e.start, reverse=(t.strand == "-"))
    blocks, off = [], 0
    for e in exons:
        n = e.end - e.start + 1
        lo, hi = max(a, off), min(b, off + n)
        if lo < hi:
            if t.strand == "+":
                blocks.append((e.start + lo - off, e.start + hi - 1 - off))
            else:
                blocks.append((e.end - (hi - 1 - off), e.end - (lo - off)))
        off += n
    return sorted(blocks)


orf_rows = []
for o in orf_ids:
    prot = cand.at[o, "sequence"]
    peps = ok[ok.orf_ids.str.split(",").map(lambda xs: o in xs)]
    placements_all = []
    for tid in cand.at[o, "transcripts"].split(","):
        t = genome.transcript_by_id(tid)
        for a, b in placements(t, prot):
            bl = blocks_for(t, a, b)
            assert sum(y - x + 1 for x, y in bl) == 3 * (len(prot) + 1)
            start = bl[0][0] if t.strand == "+" else bl[-1][1]
            placements_all.append(dict(transcript_id=tid, transcript_biotype=t.biotype, gene_id=t.gene_id,
                                       gene_name=t.gene_name, contig=t.contig, strand=t.strand,
                                       coding_blocks=[list(x) for x in bl], start_codon_position=start,
                                       transcript_offset=a))
    rep = placements_all[0]
    loci = {(p_["contig"], p_["strand"], tuple(map(tuple, p_["coding_blocks"]))) for p_ in placements_all}
    pep_evidence = []
    for _, r in peps.iterrows():
        ivals = [[m.start() + 1, m.start() + len(r.peptide)] for m in re.finditer("(?=" + r.peptide + ")", prot)]
        pep_evidence.append(dict(method="proteomics", sequence=r.peptide, protein_intervals=ivals, source="sage_reanalysis",
                                 accession="PXD052932", mapping_count=len(r.orf_ids.split(",")),
                                 spectrum_count=int(r.n_psms_q01), datasets=sorted(json.loads(r.per_run_psms_q01)),
                                 sequence_match="exact", attribution="source_orf", quality="reported",
                                 spectrum_ids=r.spectrum_ids.split(";") if r.spectrum_ids else [],
                                 q_peptide=r.q_peptide, best_q_psm=r.best_q_psm, best_hyperscore=r.best_hyperscore,
                                 classification=r.classification))
    hs, mm = orf_hom[(o, "human")], orf_hom[(o, "mouse")]
    runs_o = sorted({x for e in pep_evidence for x in e["datasets"]})
    cls_pep = {e["classification"] for e in pep_evidence}
    orf_rows.append(dict(
        orf_id=o, kind=cand.at[o, "kinds"], orf_length_aa=len(prot), protein_sequence=prot, start_codon="ATG",
        transcript_id=rep["transcript_id"], gene_id=rep["gene_id"], gene_name=rep["gene_name"] or "",
        transcript_biotype=rep["transcript_biotype"], contig=rep["contig"], strand=rep["strand"],
        coding_blocks=json.dumps(rep["coding_blocks"]), start_codon_position=rep["start_codon_position"],
        n_genomic_placements=len(loci), all_placements=json.dumps(placements_all),
        genome_build="ROS_Cfam_1.0", annotation="Ensembl 110", taxon_id=9615, species="canis_lupus_familiaris",
        translation_table=1,
        n_peptides=len(pep_evidence), n_psms_q01=sum(e["spectrum_count"] for e in pep_evidence), n_runs=len(runs_o),
        tissues=",".join(sorted({tissue_of(r) for r in runs_o})), peptides=",".join(e["sequence"] for e in pep_evidence),
        peptide_evidence=json.dumps(pep_evidence),
        classification=("conserved_coding_probable_annotation_gap"
                        if "conserved_coding_probable_annotation_gap" in cls_pep else "candidate_novel_ncORF"),
        orf_homologous_to_annotated_protein=orf_homologous(o),
        human_best_hit=hs["hit"], human_score=hs["score"], human_identity=hs["ident"], human_aln_len=hs["aln_len"],
        mouse_best_hit=mm["hit"], mouse_score=mm["score"], mouse_identity=mm["ident"], mouse_aln_len=mm["aln_len"],
    ))
orfs = pd.DataFrame(orf_rows)
if len(orfs):
    orfs = orfs.sort_values(["classification", "n_psms_q01"], ascending=[True, False])
    orfs.to_csv(os.path.join(out, "supported_orfs.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "supported_orfs.jsonl"), "w") as f:
        for r in orfs.to_dict("records"):
            for k in ("coding_blocks", "all_placements", "peptide_evidence"):
                r[k] = json.loads(r[k])
            f.write(json.dumps(r) + "\n")
    log(f"supported ORFs (ORF-unique passing peptides): {len(orfs)}; by classification: "
        f"{orfs.classification.value_counts().to_dict()}; candidate_novel ORFs that are nonetheless homologous to an "
        f"annotated human/mouse protein at ORF level: {int(((orfs.classification == 'candidate_novel_ncORF') & orfs.orf_homologous_to_annotated_protein).sum())}")
    log(f"ORFs by kind: {orfs.kind.value_counts().to_dict()}")
    nov = orfs[(orfs.classification == "candidate_novel_ncORF")]
    log(f"candidate novel ncORFs: {len(nov)}; with >=2 peptides: {int((nov.n_peptides >= 2).sum())}; "
        f"seen in >=2 runs: {int((nov.n_runs >= 2).sum())}; without ORF-level homology: "
        f"{int((~nov.orf_homologous_to_annotated_protein).sum())}")
with open(os.path.join(out, "aggregate_summary.json"), "w") as f:
    json.dump(dict(summary=summary, canonical=canon_pep_counts, per_run_canonical_psms_q01=per_run_canon.to_dict(),
                   runtime_s=round(time.time() - T0)), f, indent=1)
with open(os.path.join(out, "aggregate_summary.txt"), "w") as f:
    f.write("\n".join(summary) + "\n")
