# Canine vendor evidence extraction

This work implements ORFdata #23. The current code supports bounded native-table
export with durable checkpoints, native numeric-vector decoding, searched-reference
peptide matching and checked dog pool/run lineage in the common sample/assay API.
The full normalized dataset export, source-count reconciliation and release are pending.

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

The selected biological tables are pinned in `dog_result_plan.json` against the
inspected schema hash. They retain target and decoy identities, peptides/PSMs,
protein groups and every selected membership/modification link. Only acquisitions
referenced by those links are selected, using indexed composite-key lookups.
Spectral peak arrays and workflow/UI payloads are excluded. Existing protein-summary
abundances remain in their separate released dataset.

```sh
python -m tools.export_dog_results /path/to/audit/pdResult \
  --transfer-mib 4096 --cache-mib 1024 --scratch-mib 8192 --output-mib 2048 \
  --seconds 3600 --requests 30000 --read-ahead-kib 256
```

`native-evidence.partial.sqlite` is a staging file, not the consumer artifact.
Every batch commits exact native values, a source cursor and a typed row digest
in one transaction. Repeating the command resumes committed batches. Native
row IDs can be negative or noncontiguous; composite acquisition keys retain
int64 precision. Missing referenced acquisitions are recorded explicitly. An
output marked `selection_complete` has completed only the selected native facts;
it has not passed normalized dataset acceptance. Output and rollback journal
capacity are included in the scratch reservation. Do not open another range
reader on the same cache while an export is live.

Sequential scans use bounded read-ahead; index lookups use individual 4 KiB
pages. The source control connection can be reused only after a verified FTP
completion (including MassIVE's specifically observed bounded-transfer reset),
a successful NOOP and fresh SIZE/MDTM checks. Each acquired chunk still gets
independent offset/overlap verification. Permission errors, other transfer
failures and unsynchronized control responses fail the read.

Dog lineage is separately reproducible from the small deposited study XML,
the publication's supplementary workbook and the existing protein summary:

```sh
python -m tools.export_dog_lineage study.pdStudy supplementary-tables.xlsx \
  t04_protein_postbridge_ibaq_pd.txt native-evidence.partial.sqlite lineage.json
```

The command verifies pinned hashes and requires all 97 `WorkflowInputFiles`
rows already exported. Table S2's 44 dog channels comprise 38 tissue pools,
four 131C bridges, and empty 131N channels in Plex1 and Plex4. All 38 summary
samples match the complete published channel/tissue profiles. The native named
P1–P4 batches are associated with publication Plex1–Plex4 by their labels; that
basis is explicit rather than presented as a shared source identifier. Native
run IDs join via StudyFileID, with basename, size and fraction checks for all
96 raw files. The source's generic `Sample` role remains distinct from the
publication's tissue/bridge/empty role.

`orfdata.sources.pfizer_lineage.lineage_records` exposes 44 `Sample` records and
44 aggregate plus 1,056 run/channel `Assay` records. Fractions share their parent
pool, and no individual-dog positives or additional biological replicates are
inferred. All 96 native instrument fields report Orbitrap Fusion Lumos. The
publication's Eclipse report and the deposit's Lumos report remain separately
cited, with the discrepancy explicit.

The supplementary workbook was retrieved through Europe PMC's supplementaryFiles
endpoint for PMC13482810. Its SHA-256 is
`a95b7a2514679ba6f7283573059f27e6182b5620f6d0beead601e47daf85e740`;
its MD5 matches the article XML's `e354f58534b5c41c7819285605a2ad97`.
Article/supplement terms are CC BY-NC-ND 4.0, distinct from the MassIVE deposit's
CC0 terms. The original workbook is not included in this source tree or ordinary
dataset artifact; the fixture retains the dog channel facts and their cell
locations, provenance and hashes.

`orfdata.sources.pfizer_vendor.DistributionDecoder` reads column types, map IDs
and ordered boxes from the native database metadata. It keeps numeric presence
flags, null vectors and nonfinite values explicit. It does not equate native
reporter height, PSM abundance, normalized abundance or the separately released
iBAQ summary. Protein assignment must retain the reported members and the full
set of sequence candidates in the checksum-pinned searched FASTA separately.

`NativeQuantification.peptide_abundances` attaches each source PSM or peptide-group
array cell to the checked common sample/assay records. It retains native quality
codes, modification text, original record IDs, source column definitions and
numeric presence flags. PSM acquisition IDs must agree with the checked run;
unknown runs/channels fail instead of acquiring tissue labels from a filename.
Raw, normalized and scaled peptide-group abundances have separate units. The
common missing status carries explicit native null/nonfinite provenance, and a
measured zero remains zero. Offline native examples round-trip through
`QuantitativeDataset.create` and its unit/status/target query selectors.

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
