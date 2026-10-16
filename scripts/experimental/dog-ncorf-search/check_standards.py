"""Self-test of the retention-time standard handling in standards.py and build_db.py, on tiny
made-up inputs. Needs numpy only; no network, no Ensembl files.

usage: python -I check_standards.py
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_db import filter_candidates, load_canonical, write_outputs  # noqa: E402
from standards import IRT, PRTC, STANDARDS, contained_standards, standard_of  # noqa: E402

assert len(IRT) == 11 and len(PRTC) == 15 and len(set(IRT + PRTC)) == 26
try:  # The orfdata package flags the same iRT peptides in the dog DLA ligand dataset.
    from orfdata.sources.kaabinejadian import IRT_PEPTIDES
    assert IRT == tuple(IRT_PEPTIDES), "IRT differs from orfdata.sources.kaabinejadian.IRT_PEPTIDES"
except ImportError:
    pass
assert standard_of("LGGNEQVTR") == ("iRT", "LGGNEQVTR")
assert standard_of("NGFLLDGFPR") == ("PRTC", "NGFILDGFPR")  # I and L are equal.
assert standard_of("LGGNEQVT") is None  # Fragments aren't standards.
assert contained_standards("MKLGGNEQVTRAAA") == ["LGGNEQVTR"]

# A dog protein that happens to contain the PRTC peptide IGDYAGIK.
protein = "MSTNPKPQRKTKRNTNRRPQDVKFPGGIGDYAGIKVLLAAQW"
with tempfile.TemporaryDirectory() as tmp:
    pep = os.path.join(tmp, "pep.fa")
    pa = os.path.join(tmp, "pa.fa")
    with open(pep, "w") as f:
        f.write(f">ENSCAFP1.1 pep\n{protein}*\n>ENSCAFP2.1 pep\nMAAAAAKWWWWWR*\n")
    with open(pa, "w") as f:
        f.write(f">P1 UniProt\n{protein}\n>DECOY_P1\n{protein[::-1]}\n>CONTAM_TRYP_PIG\nMKWVTFISLLLLFSSAYS\n")
    canon, n, synthetic = load_canonical(pep, pa)
    assert n == {"ensembl": 2, "peptideatlas": 1, "contaminant": 1, "standard": 26}, n
    seqs = list(canon)
    # Proteins first, in input order, merged by sequence; then one entry per standard.
    assert seqs[:3] == [protein, "MAAAAAKWWWWWR", "MKWVTFISLLLLFSSAYS"], seqs[:3]
    assert canon[protein] == ["ens|ENSCAFP1.1", "pa|P1"]
    assert seqs[3:] == [p for _, p in STANDARDS]
    assert canon["LGGNEQVTR"] == ["std|iRT|LGGNEQVTR"]
    assert synthetic == {p for _, p in STANDARDS} - {"IGDYAGIK"}

    cands = {
        "MPLGGNEQVTRSSSSA": [dict(contig="1", start=10)],      # Has a synthetic standard: dropped.
        "MPPIGDYAGIKSSSSA": [dict(contig="1", start=20)],      # Shares a dog peptide: kept.
        "MSTNPKPQRKTKRNT": [dict(contig="1", start=30)],       # Inside a dog protein: dropped.
        "MAAAAAKWWWWWR": [dict(contig="1", start=40)],         # Identical to one: dropped.
        "MQQQQQQQQQQQQ": [dict(contig="2", start=5)],          # Kept.
    }
    keep, removed = filter_candidates(cands, canon, synthetic)
    assert set(keep) == {"MPPIGDYAGIKSSSSA", "MQQQQQQQQQQQQ"}, set(keep)
    assert removed == {"contains_standard": 1, "contained": 1, "identical": 1}, removed

    for s in keep:
        keep[s] = [dict(o, kind="lncORF", tx="T", gene="G", gene_name="", stop_end=0, strand="+")
                   for o in keep[s]]
    write_outputs(tmp, canon, keep)
    with open(os.path.join(tmp, "dog_db.fasta")) as f:
        headers = [line.rstrip() for line in f if line.startswith(">")]
    assert headers[0] == ">CANON000000 ens|ENSCAFP1.1 pa|P1"
    assert headers[3] == ">CANON000003 std|iRT|LGGNEQVTR"
    assert headers[-2:] == [">NCORF000000 lncORF G", ">NCORF000001 lncORF G"]
print("standards checks passed")
