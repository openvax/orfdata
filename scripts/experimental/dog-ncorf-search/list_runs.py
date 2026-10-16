"""List the mzML runs of PRIDE PXD052932 (Canine PeptideAtlas) as runs.tsv, the run list
driver.sh works through.

Reads each archive's ZIP central directory with HTTP range requests (a few hundred KB per
archive), or a saved listing given with --listing. Writes one row per .mzML member, with no
header: run, group, zip_url, member, mzml_bytes, compressed_bytes. Other members (the
Heart S-Trap runs' .pin files, directory entries) are reported on stderr and skipped.

usage: python -I list_runs.py [--base-url URL] [--group NAME ...] [--listing GROUP=PATH ...]
                              [--save-listings DIR] [--match REGEX] OUT_TSV
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remotezip import listing  # noqa: E402

BASE_URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD052932/"
# In the order driver.sh searches them: smallest runs first, so problems show up early.
GROUPS = (
    "mzML_QExactive_plasma_urine",
    "mzML_timsTOF_tissue_PAC_DDA",
    "mzML_timsTOF_tissue_STrap_DDA",
    "mzML_timsTOF_plasma_DDA",
)


def read_listing(path):
    """Rows of a saved listing: member, size, compressed size, method (tab-separated)."""
    with open(path) as f:
        return [(m, int(u), int(c), int(t)) for m, u, c, t in
                (line.rstrip("\n").split("\t") for line in f if line.strip())]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="runs.tsv to write")
    ap.add_argument("--base-url", default=BASE_URL, help="directory URL holding GROUP.zip (default: %(default)s)")
    ap.add_argument("--group", action="append", help="archive name without .zip; repeatable (default: all four)")
    ap.add_argument("--listing", action="append", default=[], metavar="GROUP=PATH",
                    help="use a saved listing (output of `remotezip.py list`) instead of the network")
    ap.add_argument("--save-listings", metavar="DIR", help="also save each archive's listing as DIR/GROUP.tsv")
    ap.add_argument("--match", metavar="REGEX", help="keep only runs whose name matches (re.search)")
    a = ap.parse_args()

    saved = dict(x.split("=", 1) for x in a.listing)
    groups = a.group or list(GROUPS)
    unknown = set(saved) - set(groups)
    if unknown:
        raise SystemExit(f"--listing for groups not being listed: {sorted(unknown)}")
    rows, skipped = [], []
    for group in groups:
        url = a.base_url.rstrip("/") + "/" + group + ".zip"
        members = read_listing(saved[group]) if group in saved else listing(url)
        if a.save_listings:
            os.makedirs(a.save_listings, exist_ok=True)
            with open(os.path.join(a.save_listings, group + ".tsv"), "w") as f:
                f.writelines("\t".join(map(str, m)) + "\n" for m in members)
        for member, size, csize, _ in members:
            if not member.endswith(".mzML"):
                skipped.append(member)
                continue
            run = os.path.basename(member)[: -len(".mzML")]
            if a.match and not re.search(a.match, run):
                continue
            rows.append((run, group, url, member, size, csize))
    runs = [r[0] for r in rows]
    if len(set(runs)) != len(runs):
        raise SystemExit("run names are not unique across archives")
    with open(a.out, "w") as f:
        f.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    for m in skipped:
        print(f"skipped (not mzML): {m}", file=sys.stderr)
    print(f"{len(rows)} runs, {sum(r[5] for r in rows) / 1e9:.1f} GB compressed, "
          f"{sum(r[4] for r in rows) / 1e9:.1f} GB mzML -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
