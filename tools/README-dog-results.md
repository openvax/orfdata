# Canine vendor evidence extraction

This work implements ORFdata #23. The current code supports bounded schema
inspection, durable acquisition accounting, native numeric-vector decoding and
complete searched-reference peptide matching in bounded batches. The full
normalized dataset export, source-count reconciliation and release are pending.

Use a separate extraction environment with `.[canine-extraction]` installed.
The current audit uses APSW 3.54.1.0 (SQLite 3.54.1) and pyahocorasick 2.3.1.
Neither normal imports nor quantitative
dataset queries call this tool or acquire the vendor database.

From this checkout, an explicit schema audit is:

```sh
python -m tools.inspect_dog_results pdResultView /path/to/audit/pdResultView \
  --transfer-mib 64 --cache-mib 32 --scratch-mib 128 --seconds 600
```

Supported pinned assets are `pdResultView`, `msfView` and `pdResult` in public
MassIVE MSV000099744. The small views contain subset/sorting IDs; the main
database contains the biological evidence. The source size and MDTM are checked
before and after range acquisition. Each new chunk is checked against an
independent shifted overlapping request before entering the checksum cache.
Cached bytes are checked on every reuse. SQLite performs all page, index and
overflow traversal through a read-only APSW VFS; no second page parser is used.

Transfer and request limits are cumulative for a cache. Each range is durably
reserved before opening the stream; a crash leaves its full reservation charged
until actual receipt is known. Completed and failed attempts settle to their
received payload. Reopening a cache does not reset these counters. Migrating an
older cache requires every prior receipt explicitly with `--prior-receipt`;
duplicate receipts or a different source identity are rejected. Never supply
both a numbered receipt and its identical latest-receipt copy.

A process lock prevents concurrent use. Pending chunk writes recover only when
their saved length/hash match, and all reused chunks are verified again. The
`--chunk-kib` option can be reduced on resumption without splitting or discarding
existing verified chunks; 4 KiB reads avoid overfetch when index leaves are
scattered among large data records. Increasing the chunk size on that cache is
rejected because it could overlap its existing finer chunks.

Limits cover received payload, cached payload, overall scratch/output bytes,
request count, elapsed time and SQLite allocation. Scratch includes allocated
file blocks, metadata, atomic-write duplicates, and a reserved upper bound for
the transfer ledger and rollback journal. Receipts measure process peak RSS,
current and peak-observed disk use, and reserved capacity. They include tool
source fingerprints and runtime versions. Received payload includes retries and overlap verification, but excludes
TLS/control overhead and bytes sent ahead by the server. Range hashes and the
header/metadata pin do not constitute an authoritative whole-database checksum.
The elapsed-time limit applies to each invocation; payload/request counters
persist. The final export must use the same scratch-budget object for all output
writes and bound its SQLite output database/journal capacity too.

`orfdata.sources.pfizer_vendor.DistributionDecoder` reads column types, map IDs
and ordered boxes from the native database metadata. It keeps numeric presence
flags, null vectors and nonfinite values explicit. It does not equate native
reporter height, PSM abundance, normalized abundance or the separately released
iBAQ summary. Protein assignment must retain the reported members and the full
set of sequence candidates in the checksum-pinned searched FASTA separately.

`peptide_candidates` takes a verified `SearchReference` and a batch capped by
total peptide residues. It streams all accessions and repeated sites, using
zero-based, end-exclusive protein offsets. Exact residue matches and optional
isoleucine/leucine-equivalent matches have distinct labels. A sequence with no
match remains explicit; the target FASTA does not establish a complete decoy or
canine isoform background. Candidates never replace reported group memberships.

The offline tests include SQLite index/interior/overflow pages, failed/short and
wrong-offset reads, process interruption, source changes, corrupted cache,
restart/temporary-file quotas, shared peptide attribution, real 11-channel
PSM and 44-sample protein arrays, missing quantification and malformed metadata.
The native fixture's provenance is in `tests/data/pfizer/vendor/NOTICE.md`.

Source references:

- https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=062b109e9adf461ba619369dd0d85aa9
- https://doi.org/10.1038/s41597-026-07512-9
- https://rogerbinns.github.io/apsw/vfs.html
- https://www.sqlite.org/fileformat.html
- https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1297839115v1
- https://github.com/BenBruyneel/proteinDiscover/blob/main/R/transformRaw.R
