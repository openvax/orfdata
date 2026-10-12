# Native protein search references

A search reference lists the sequences available to a particular MS search.
It is independent of reported abundance or experimental evidence that a protein
exists. ORFdata exposes the exact deposited FASTA, including unreported entries.

| Dataset | Search entries | FASTA bytes, downloaded and installed |
| --- | ---: | ---: |
| `pfizer-atlas-human-protein` | 20,600 | 13,831,073 |
| `pfizer-atlas-mouse-protein` | 21,989 | 14,325,720 |
| `pfizer-atlas-dog-protein` | 20,649 | 14,659,832 |

These are uncompressed, no-isoform UniProt search sets from data release
`2026-10-12.1`. The historical UniProt release and genomic assembly remain
unknown. ORFdata does not replace them with a current proteome.

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()
name, version = "pfizer-atlas-dog-protein", "2026-10-12.1"
spec, = registry.search_references(name, version)  # Metadata only; no I/O.
registry.download_reference(name, version)        # Explicit verified download.
search = registry.open_reference(name, version)   # Offline; raises if missing.

entry = search.get("A0A068LEV2")
print(entry.header, entry.sequence, entry.sequence_sha256)
assert search.membership(entry.target).status == "present"
for entry in search.entries():
    assert entry.target.reference == search.reference

absent = ScopedTarget(search.reference, "protein", "not-in-this-search", "uniprot")
assert search.membership(absent).status == "not_in_search"
```

Dog has 9,929 reported group representatives and 10,720 other search entries,
including `A0A068LEV2`. Those other entries have **no reported representative**;
that does not prove that their peptides were undetected. The summary table does
not resolve all protein-group members. Query group abundance separately through
[`dataset.measurements()`](pfizer-protein.md).

Membership accepts an explicitly scoped native protein target. A different
reference raises an error. Other namespaces, ORFs, genes and protein groups
return `unknown`: first supply an explicit, validated reference mapping. Search
membership never creates a zero, allocates a group's abundance, establishes
protein detectability, or supplies MHC presentation evidence. Exact accession
lookup does not perform sequence matching or infer coverage of a new ORF.

Each entry retains its full header, FASTA ID, residues, sequence hash, source
taxon, optional gene label and source line. `search.spec.to_dict()` records the
dataset version, reference identity, sequence-set SHA-256, source URL, size,
licence and attribution. The deposited records do not specify decoy/control
roles; those fields stay `None`. UniProt records retain CC-BY-4.0 attribution.

```sh
orfdata datasets references pfizer-atlas-dog-protein --version 2026-10-12.1
orfdata datasets download-reference pfizer-atlas-dog-protein --version 2026-10-12.1
orfdata datasets inspect-reference pfizer-atlas-dog-protein --version 2026-10-12.1
orfdata datasets query-reference pfizer-atlas-dog-protein --version 2026-10-12.1 --target A0A068LEV2
```

Pin the quantitative dataset version and retain the matching search spec when
exporting an analysis. Both use the same native `ReferenceContext`; select
`reference=spec.reference` and `asset_id=spec.source.asset_id` if a dataset has
multiple declared search sets (`--reference KEY --asset ID` in the CLI).
Reference downloads use a separate DataCache store and never rebuild or download
abundance. Local source overrides are checksum-verified; `force=True` creates a
new immutable bundle and preserves already-open readers. Ordinary discovery,
inspection and opening stay offline. Published quantitative hashes are unchanged.

The same reference model accepts any taxon, with no inferred cross-species joins.
Adapters must explicitly declare their native search sets and formats. The
initial parser supports deposited UniProt FASTAs. CoPY yeast's QconCAT FASTA is
synthetic calibration material, so it is not exposed as a native yeast proteome;
that historical search reference remains an explicit source gap.
