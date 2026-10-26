# Canine vendor evidence extraction

This work implements ORFdata #23. The current code supports bounded schema
inspection and tested native numeric-vector decoding. The full normalized
dataset export, source-count reconciliation and release are still pending.

Use a separate extraction environment with APSW 3.54.1.0 (SQLite 3.54.1) and
the ORFdata runtime/test dependencies. Neither normal imports nor quantitative
dataset queries call this tool or acquire the vendor database.

From this checkout, an explicit schema audit is:

```sh
python -m tools.inspect_dog_results pdResultView /path/to/audit/pdResultView \
  --transfer-mib 64 --cache-mib 32 --seconds 600
```

Supported pinned assets are `pdResultView`, `msfView` and `pdResult` in public
MassIVE MSV000099744. The small views contain subset/sorting IDs; the main
database contains the biological evidence. The source size and MDTM are checked
before and after range acquisition. Each new chunk is checked against an
independent shifted overlapping request before entering the checksum cache.
Cached bytes are checked on every reuse. SQLite performs all page, index and
overflow traversal through a read-only APSW VFS; no second page parser is used.

Limits cover received data payload, cached chunk payload, request count, elapsed
time and SQLite allocation. Receipts measure process peak RSS and current disk
usage. Received payload includes retries and overlap verification, but excludes
TLS/control overhead and bytes sent ahead by the server. Range hashes and the
header/metadata pin do not constitute an authoritative whole-database checksum.
The prototype's transfer limit is per invocation; production extraction still
needs a durable cumulative transfer ledger and overall scratch/output quotas.

`orfdata.sources.pfizer_vendor.DistributionDecoder` reads column types, map IDs
and ordered boxes from the native database metadata. It keeps numeric presence
flags, null vectors and nonfinite values explicit. It does not equate native
reporter height, PSM abundance, normalized abundance or the separately released
iBAQ summary. Protein assignment must retain the reported members and the full
set of sequence candidates in the checksum-pinned searched FASTA separately.

The offline tests include SQLite index/interior/overflow pages, failed/short and
wrong-offset reads, source changes, corrupted cache, quotas, real 11-channel
PSM and 44-sample protein arrays, missing quantification and malformed metadata.
The native fixture's provenance is in `tests/data/pfizer/vendor/NOTICE.md`.

Source references:

- https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=062b109e9adf461ba619369dd0d85aa9
- https://doi.org/10.1038/s41597-026-07512-9
- https://rogerbinns.github.io/apsw/vfs.html
- https://www.sqlite.org/fileformat.html
- https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1297839115v1
- https://github.com/BenBruyneel/proteinDiscover/blob/main/R/transformRaw.R
