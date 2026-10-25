"""Record what one source publishes for every entry in a collection's inventory.

A reference inventory can name an entry more than once — here, a GenBank
accession and a PaVE reference — and those are not always the same molecule.
GDC's decoy FASTA carries the PaVE one: PaVE's ``HPV16REF`` reproduces the GDC
HPV16 DNA MD5 pinned in ``tools/build_hpv16_cds.py`` exactly, while GenBank's
``K02718.1`` is a different, two-base-shorter sequence.

The ledger this writes pins each molecule's length and checksums and counts
its coding features, so a collection can tell offline which source it is
looking at, verify what it downloads, and measure one source's annotation
against another's. One ledger per source, all in one format:

    python tools/build_reference_ledger.py --source genbank \
        orfdata/data/viral/gdc-grch38-d1-vd1-genbank.tsv
    python tools/build_reference_ledger.py --source pave \
        orfdata/data/viral/gdc-grch38-d1-vd1-pave.tsv

Rebuild a ledger whenever the source is re-checked, or whenever ORFdata's own
readers change what they count as a coding feature.
"""

import argparse
import hashlib
import sys
from datetime import datetime, timezone

from orfdata.reference_ledger import PUBLISHED, read_ledger_bytes, write_ledger
from orfdata.viral_coverage import viral_panel_coverage
from orfdata.viral_registry import viral_collection


def build(collection_name="gdc-viruses", source_name=None, *, timeout=60):
    """Fetch every record one source publishes for an inventory, with its checksums."""
    collection = viral_collection(collection_name)
    source = collection.source(source_name)

    def report(entry, error):
        print(f"  {entry.abbreviation}: {error}", file=sys.stderr)

    ledger = source.build_ledger(
        viral_panel_coverage(collection=collection), timeout=timeout, report=report
    )
    for expected in collection.expected_for(source.name):
        expected.check(ledger)
    checked = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    notes = {f"expected_{expected.abbreviation.lower()}": expected.sequence_md5
             for expected in collection.expected_for(source.name)}
    content = write_ledger(ledger, panel=collection.panel, source=source.name,
                           source_url=source.url, checked=checked,
                           notes=notes)
    read_ledger_bytes(content, name=f"{collection.panel} {source.name} ledger")
    serving = int((ledger["state"] == PUBLISHED).sum())
    print(f"{source.name}: {serving} of {len(ledger)} entries published", file=sys.stderr)
    return content


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="Where to write the ledger")
    parser.add_argument("--collection", default="gdc-viruses")
    parser.add_argument("--source", help="Which of the collection's sources to record")
    parser.add_argument("--timeout", type=float, default=60)
    arguments = parser.parse_args()
    content = build(arguments.collection, arguments.source, timeout=arguments.timeout)
    with open(arguments.output, "wb") as handle:
        handle.write(content)
    print(f"{arguments.output}: {len(content)} bytes, "
          f"sha256 {hashlib.sha256(content).hexdigest()}", file=sys.stderr)


if __name__ == "__main__":
    main()
