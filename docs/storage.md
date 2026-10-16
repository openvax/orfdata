# Quantitative storage

ORFdata 0.8 writes compact SQLite storage version 4. `Measurement`,
`QuantitativeDataset`, the registry and CLI keep the same public records and
query semantics. Human, mouse, dog, yeast and arbitrary taxa retain separate
reference-scoped identities.

References, targets and assays are stored once and addressed by internal integer
keys. Measurements retain their values and metadata; repeated required fields
are reconstructed from indexed columns and the dataset manifest. Fields equal
to the fixed version-4 defaults are omitted from stored JSON and restored on read.
Zeros, nulls, censoring, uncertainty, nested provenance and derived inputs retain
their meaning. Foreign keys prevent cross-reference target/assay assignments.
Internal integer keys never appear as public biological identifiers.

ORFdata 0.32 writes storage version 5, which also stores each target mapping's
candidates as kind, identifier and namespace: they always share the mapping's
reference, which its row already names. The dog DLA ligand dataset, with 62,492
mappings, is 255 MB, against 389 MB in the version 4 layout. Version 4 files remain readable, and a
version 4 fixture written by the released ORFdata 0.30.0 tests that.

Published storage versions 1, 2, 3 and 4 remain readable. There is no automatic
migration or replacement of downloaded artifacts. A new source rebuild can use
a different storage layout and therefore a different file checksum while
returning the same records. Registry downloads always verify the exact published
artifact's hash; source builds retain their own DataCache receipts.

The implementation streams observations and bounds its target cache to 131,072
entries. SQLite index statistics are populated during creation so a narrow
target query can use the target index instead of scanning an entire reference.
Raw SQL table layouts are internal; use the Python/CLI query and export APIs.

## Actual-source benchmark

The first 100,000 rows of the native Pfizer/GEO RNA table contain four dog
libraries and 30,941 genes, including measured zeros. Source gzip SHA-256:
`a22780f92ffd0a028e30fe3392d73d01f787eb5a62b2043e8d26940191cbac70`.
This is a storage benchmark slice, not a newly released experimental dataset.

| Measurement | Storage v3 | Storage v4 |
|---|---:|---:|
| Installed size | 121.6 MB | 38.6 MB |
| Gzip size | 13.4 MB | 7.5 MB |
| Build time | 8.90 s | 4.28 s |
| One-gene query (four results) | 135 ms | 1.5 ms |
| Peak build RSS | 39.8 MB | 52.1 MB |

All 100,000 reconstructed observations were compared field for field. Timings
are local measurements on Python 3.12, not cross-machine guarantees; the bounded
target cache trades some memory for repeated-target lookup speed. Legacy tests
also use actual databases written by ORFdata 0.3.0, 0.4.0 and 0.5.0.
