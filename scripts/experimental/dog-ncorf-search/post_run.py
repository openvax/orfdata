"""Post-process one Sage run: keep its rank-1 PSMs with all of Sage's columns (zstd parquet), and
append the run's receipt (download receipt plus PSM counts) to the run manifest. Writes the
manifest header when the manifest doesn't exist yet. Prints the appended row as JSON.

usage: python -I post_run.py RUN_ID SAGE_TSV RECEIPT_JSON SAGE_SECONDS PARQUET_OUT MANIFEST_TSV
"""
import json
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from standards import standard_of  # noqa: E402

STRIP = r"\[[^\]]*\]|[^A-Z]"


def main():
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    run, tsv, receipt, sage_s, parquet_out, manifest = sys.argv[1:7]
    with open(receipt) as f:
        rec = json.load(f)
    d = pd.read_csv(tsv, sep="\t")
    n_rows = len(d)
    d = d[d["rank"] == 1].copy()
    d.insert(0, "run", run)
    d.to_parquet(parquet_out + ".part", index=False, compression="zstd")
    os.rename(parquet_out + ".part", parquet_out)

    t = d[d.label == 1]
    canon = d.proteins.str.split(";").map(lambda ps: any(p.replace("rev_", "").startswith("CANON") for p in ps))
    # iRT peptides only: no dog protein contains one, while four PRTC peptides are also dog
    # peptides (aggregate.py's run_qc.tsv tells those apart).
    irt = t.peptide.str.replace(STRIP, "", regex=True).map(lambda p: (standard_of(p) or ("",))[0] == "iRT")
    row = dict(
        run=run, member=rec["member"], zip_url=rec["zip_url"], compressed_bytes=rec["compressed_bytes"],
        mzml_bytes=rec["mzml_bytes"], mzml_crc32=rec["mzml_crc32"], crc_ok=rec["crc_ok"],
        mzml_sha256=rec["mzml_sha256"], n_ms2=rec["n_ms2"], download_s=rec["download_s"],
        download_MBps=rec["MBps"], sage_s=float(sage_s), sage_rows=n_rows, rank1_psms=len(d),
        target_psms_q01=int((t.spectrum_q <= 0.01).sum()),
        target_peptides_q01=int(t[t.peptide_q <= 0.01].peptide.nunique()),
        target_proteins_q01=int(t[t.protein_q <= 0.01].proteins.nunique()),
        canonical_rank1_targets=int((canon & (d.label == 1)).sum()),
        canonical_rank1_decoys=int((canon & (d.label == -1)).sum()),
        ncorf_rank1_targets=int((~canon & (d.label == 1)).sum()),
        ncorf_rank1_decoys=int((~canon & (d.label == -1)).sum()),
        psm_parquet=os.path.basename(parquet_out),
        completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        irt_psms_q01=int((irt & (t.spectrum_q <= 0.01)).sum()),
    )
    new = not os.path.exists(manifest)
    with open(manifest, "a") as f:
        if new:
            f.write("\t".join(row) + "\n")
        f.write("\t".join(str(v) for v in row.values()) + "\n")
    print(json.dumps(row))


if __name__ == "__main__":
    main()
