# Viral coding references

`viral_orfs()` returns existing `TranscriptORF` objects for source-annotated
viral CDSs. Supply the complete reference DNA that the RNA was aligned to:

```python
from orfdata import viral_orfs, viral_references

print([(r["id"], r["coding_models"]) for r in viral_references()])
# [('NC_001526.4', 9), ('GRCh38.d1.vd1:HPV16', 9)]

orfs = viral_orfs(reference_dna, contig="HPV16")
for orf in orfs:
    print(orf.orf_id, orf.gene_name, orf.coding_blocks, orf.protein_sequence)
```

The small bundled catalogue loads offline without PyEnsembl, Biopython or a
human genome. Its nine HPV16 CDSs come from RefSeq NC_001526.4. The loader
requires an exact full-DNA hash and separately qualifies native and GDC
circular rotations. A contig name or related strain is insufficient.
An unavailable reference raises `ValueError`, distinct from a valid empty query.
`viral_references()` exposes exact identities, coding coverage, source pins
and product names without loading other data.

Each model uses a `SourceRNA` containing only the available CDS, including
ordered splice blocks. Full transcript ends are not asserted. The existing
`TranscriptORF.at()` validates its source protein and coding boundaries.
These models have `orf_type="reference_CDS"`. Their internal kind `extra_orf`
describes a source-defined ORF outside a loaded PyEnsembl annotation, not a
claim that a RefSeq CDS is biologically noncanonical. No Ribo-seq or MS
support is inferred from the reference annotation.

The nine annotations are not an exhaustive HPV16 ORF catalogue. Other HPV
types/lineages, EBV/KSHV, CMV, hepatitis viruses including HDV, HTLV-1/HIV,
and polyomaviruses remain explicit coverage work in
[issues #91–98](https://github.com/openvax/orfdata/issues/91).

Isovar owns reconstructed sample RNA and read evidence. It can use these
coding paths to anchor a frame on observed sequence, preserving variants
and incomplete starts/stops. ORFdata does not assemble sample RNA, replace
missing bases, infer protein abundance or establish viral tumor causality.
RNA and coding-reference matches are not direct protein evidence.

Source attribution and rebuilding instructions are in
[`orfdata/data/viral/NOTICE.md`](https://github.com/openvax/orfdata/blob/main/orfdata/data/viral/NOTICE.md).
