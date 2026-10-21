"""Build the metadata-only coverage ledger from GDC's published viral inventory.

The input is the original ``GRCh83.d1.vd1_virus_decoy.txt`` file. Its name
contains ``GRCh83`` even though the reference is ``GRCh38.d1.vd1``. No viral
sequences or coding models are inferred from this inventory.
"""

import argparse
import csv
import hashlib
import io
from pathlib import Path

PANEL = "GRCh38.d1.vd1"
SOURCE_URL = "https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt"
SOURCE_SHA256 = "f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56"
SOURCE_COLUMNS = ("Name", "Abbreviation", "RefSeq", "GenBank", "PaVE ID", "PaVE Notes")
OUTPUT_COLUMNS = ("index", "name", "abbreviation", "refseq", "genbank", "pave_id",
                  "pave_notes", "bundled_reference_id", "coding_models", "status")

# An inventory row becomes covered only after its GDC placement is released.
# This records the source CDS count, not a claim of exhaustive viral ORF coverage.
BUNDLED_CDS = {"HPV16": ("GRCh38.d1.vd1:HPV16", 9)}


def build(raw):
    """Return a deterministic TSV after checking the exact published input."""
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError("GDC viral inventory checksum mismatch")
    reader = csv.DictReader(io.StringIO("\n".join(raw.decode("utf-8").splitlines())),
                            delimiter="\t")
    if tuple(reader.fieldnames or ()) != SOURCE_COLUMNS:
        raise ValueError("Unexpected GDC viral inventory columns")
    rows = [row for row in reader if any(row.values())]
    if len(rows) != 200 or sum(row["Abbreviation"].startswith("HPV") for row in rows) != 189:
        raise ValueError("Expected 200 GDC viral rows, including 189 HPV rows")
    abbreviations = [row["Abbreviation"] for row in rows]
    if len(set(abbreviations)) != len(abbreviations):
        raise ValueError("Duplicate GDC viral abbreviation")
    if set(BUNDLED_CDS) - set(abbreviations):
        raise ValueError("A bundled reference is absent from the GDC inventory")

    output = io.StringIO()
    output.write("# schema: orfdata.viral_panel_coverage.v1\n")
    output.write(f"# panel: {PANEL}\n")
    output.write(f"# source: {SOURCE_URL}\n")
    output.write(f"# source_sha256: {SOURCE_SHA256}\n")
    writer = csv.DictWriter(output, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    for index, row in enumerate(rows, 1):
        abbreviation = row["Abbreviation"]
        reference_id, coding_models = BUNDLED_CDS.get(abbreviation, ("", ""))
        writer.writerow(dict(index=index, name=row["Name"], abbreviation=abbreviation,
                             refseq=row["RefSeq"], genbank=row["GenBank"],
                             pave_id=row["PaVE ID"], pave_notes=row["PaVE Notes"],
                             status="source_cds_bundled" if reference_id else "not_bundled",
                             bundled_reference_id=reference_id, coding_models=coding_models))
    return output.getvalue().encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(build(args.inventory.read_bytes()))
    print(hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
