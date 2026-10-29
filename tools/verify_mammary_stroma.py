"""Independently compare every projected value with its original source cell."""

import argparse
import csv
import gzip
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from orfdata import QuantitativeDataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("source_directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = args.source_directory
    dataset = QuantitativeDataset(args.dataset)
    rows = {}
    for asset, filename, delimiter in [
        ("rna.csv.gz", "GSE135183_raw_counts.csv.gz", ","),
        ("proteinGroups.txt", "proteinGroups.txt", "\t"),
        ("peptides.txt", "peptides.txt", "\t"),
    ]:
        opener = gzip.open if filename.endswith(".gz") else open
        with opener(sources / filename, "rt", newline="") as handle:
            rows[asset] = {
                line: row for line, row in enumerate(csv.DictReader(handle, delimiter=delimiter), 2)
            }
    counts = Counter()
    checked = 0
    started = time.monotonic()
    for observation in dataset.measurements():
        p = observation.provenance
        raw = rows[p["source_asset"]][p["line"]][p["column"]]
        assert p["native_value"] == raw
        missing = observation.quantity != "rna_abundance" and (raw == "NaN" or float(raw) == 0)
        assert observation.status == ("missing" if missing else "measured")
        assert observation.value == (None if missing else float(raw))
        if observation.quantity != "rna_abundance":
            for name, value in p["source_qc"].items():
                assert value == rows[p["source_asset"]][p["line"]][name]
        counts[observation.quantity + ":" + observation.unit + ":" + observation.status] += 1
        checked += 1
        if checked % 200000 == 0:
            print("Checked source cells", checked, flush=True)
    assert checked == 993560
    for mapping in dataset.mappings():
        p = mapping.to_dict()["provenance"]
        assert p["native_record"] == rows[p["source_asset"]][p["line"]]
    pairs = list(dataset.pairings())
    assert len(pairs) == 30
    complete = [p for p in pairs if p.right_assay_id]
    missing = [p for p in pairs if not p.right_assay_id]
    assert len(complete) == 28 and len(missing) == 2
    assert len({p.evidence["rna_case"] for p in complete}) == 14
    assert {p.evidence["rna_case"] for p in missing} == {"127"}
    assert all(p.basis == "same_specimen" and p.reference != p.right_reference for p in pairs)
    assert not list(dataset.relations()) and not list(dataset.reference_mappings())
    with dataset.path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    receipt = dict(
        dataset_id=dataset.manifest.dataset_id,
        version=dataset.manifest.version,
        path=str(dataset.path),
        size=dataset.path.stat().st_size,
        sha256=digest,
        checked_source_cells=checked,
        quantities=dict(counts),
        counts=dataset.counts,
        specimens=14,
        complete_assay_pairs=28,
        missing_mates=2,
        native_ms_rows="All 6996 original MS records exactly preserved",
        reference_mapping="Unresolved; no gene-to-protein or current-reference projection",
        missingness="Source MS 0 and NaN retained with explicit missing status; RNA zeros measured",
        elapsed_seconds=time.monotonic() - started,
    )
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
