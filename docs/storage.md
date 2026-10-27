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

## Native peptide vectors

Canine native-evidence builds use storage version 6 for the common measurement
API. Numeric vectors preserve their declared type, exact numeric bytes, independent
presence bits and null-vector state. Sample/assay axes and native source records
are stored once. Queries expand the requested cells into ordinary `Measurement`
objects, retaining missing values, source record IDs, native units and provenance.
Target and assay indexes select vectors before expansion. Counts report logical
observations rather than stored vectors. Mapped/paired queries and Parquet export
use the same projection. Opening and querying a release never accesses the remote
vendor file or requires the extraction dependencies.

Measured-only queries discard source-null vectors in SQL before expanding cells;
non-null vectors still receive per-cell presence and finite-value checks. Missing
queries and logical counts retain every null cell. On a 10,000-group spread across
the completed native table, the same 30,831 measured records (identical full-field
hash) took 30.2 seconds before this filter and 3.4 seconds after it. The query
expanded 1,758 non-null vectors instead of 30,000 total vectors. This is a bounded
query benchmark, not a full-dataset performance guarantee.

Version 6 declares its projection encoding and completion state explicitly. It
does not contain an empty flat `measurement` view: SQL access to that old internal
view fails. Generic dataset creation continues to write version 5, and versions
1–5 remain readable without migration.

The packaged canine importer also requires a final successful build receipt.
It pins normalized evidence, compiled lineage, searched FASTA and the existing
protein summary; validates the source selection and counts; builds membership
and candidate indexes; and projects the native vectors. `source_evidence()`
exposes the original records and assignments. Its `summary_group()` query connects
the baseline summary representative to native records and group links without
duplicating the baseline abundance series or asserting unique protein detection.

A bounded comparison of 100 actual native peptide groups (source row IDs 50,
100, …, 5000) reproduced all 13,200 cells from their three 44-channel abundance
vectors, field for field. Flat storage occupied 24,891,392 bytes; vector storage
occupied 6,557,696 bytes including exact selected source evidence and shared
sample/assay metadata. Peak benchmark process RSS was 127,729,664 bytes. This
is a storage slice, not a completed canine dataset release or a guarantee of
full-dataset size.

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
