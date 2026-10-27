# Canine vendor evidence extraction

This work implements ORFdata #23. The current code supports bounded native-table
export with durable checkpoints, lossless numeric-vector storage, searched-reference
peptide matching and checked dog pool/run lineage in the common sample/assay API.
A completed selected export can be normalized and assembled into a common
measurement/source-evidence artifact with the packaged importer. Full native
extraction, full artifact execution/acceptance and release remain pending;
passing a fixture is not full acceptance.

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
Every disk-loaded chunk is checked against its full hash. One immutable verified
chunk (at most 16 MiB) can serve repeated SQLite page reads while the cached file's
identity, size, mtime and ctime remain unchanged. A changed file or an evicted
buffer requires full disk verification again; resumption starts with no memory
buffer. Receipts report this buffer's residency, peak size and hits.
An offline replay of 4,096 page reads from one acquired 4 MiB chunk returned
identical 16 MiB output hashes in 7.85 seconds with repeated disk verification
and 0.048 seconds with this buffer. That measures cached page reads only;
remote transfer and complete extraction performance remain separate.
SQLite performs all page, index and
overflow traversal through a read-only APSW VFS; no second page parser is used.

Transfer and request limits are cumulative for a cache. Each range is durably
reserved before opening the stream; a crash leaves its full reservation charged
until actual receipt is known. Completed and failed attempts settle to their
received payload. An asynchronous interrupt during a receive preserves its full
reservation too, because the interruption can precede a received-byte counter
update. Reopening a cache does not reset these counters. Migrating an
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

`--batch-rows` controls rows per atomic checkpoint, with a separate 16 MiB
serialized-row ceiling. It can change on resumption without changing the selected
facts or logical digest. A local replay of all 962,755 native group/PSM links
at 500 and 5,000 rows produced the same typed digest, with 1,927 and 194 commits
respectively. Larger batches reduce repeated scratch-tree accounting; they do
not change transfer/cache/output limits. Receipts retain the requested table and
batch options, and progress reporting measures rows added since resumption.

The 42 selected biological tables (482 columns) are pinned in `dog_result_plan.json` against the
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

The selection includes direct `QuanSpectrumInfoTargetPeptideGroups` links and
their native `QuanUsage`, as well as PSM-to-quantification links. Both sets drive
referenced quantification-spectrum and reporter selection. Protein identification
group metadata preserves the source's processing-node/group names.

An earlier 40-table checkpoint can adopt these additional facts with the explicit
`--extend-selection` option. Existing table columns, source identity, cursors,
row counts and hashes remain unchanged. The acquisition selections may expand
only before their queues or rows exist. The prior configuration and all progress
records are retained in `selection_history`, atomically with the extension.
Removing facts, changing source identity, or silently relabeling a begun
acquisition selection is rejected. Finish or inspect the live process before
running an extension; never open a second cache owner.

After all selected tables finish, normalize the completed staging file explicitly:

```sh
python -m tools.normalize_dog_evidence /path/to/audit/pdResult/native-evidence.partial.sqlite \
  /path/to/normalized --output-mib 2048 --scratch-mib 8192
```

The normalizer refuses an incomplete selection or a changed extraction plan. It
retains every selected scalar, opaque byte field, metadata definition and numeric
vector, with independent source-row and restored-row digests. Referenced acquisition
tables also reproduce the extractor's composite-key digest order. Missing links
remain explicit. The output's full checksum covers the selected evidence component,
not the unacquired vendor database. Receipts fingerprint the normalizer and its
codec/budget modules, and report runtime, peak process RSS and bounded disk use.

`tools.build_dog_dataset` assembles the normalized evidence, checked `lineage.json`,
complete searched FASTA and original protein summary. It requires the entire pinned
42-table selection before creating a release manifest, and pins the already
released summary's dataset ID/version, source hash and counts. Run `python -m
tools.build_dog_dataset --help` for the explicit local input paths, intended stable
component URLs, dataset version and output/scratch byte limits. An optional
`--compressed-evidence` input pins both gzip bytes and expanded SQLite bytes.
It creates a new output directory containing `manifest.json`, `dataset.sqlite`
and a successful `build-receipt.json`; it does not publish them or add a registry entry.

The packaged `pfizer-canine-evidence-v1` importer verifies all four inputs before
and after building, checks native record counts and run identities against lineage,
retains every candidate and reported membership, and installs all native abundance
vectors. Its scratch reservation covers inputs, optional expansion and three times
the SQLite output cap. SQLite enforces the output cap from the first schema write.
Interrupted/invalid builds remove their owned new output and expanded temporary;
existing outputs are never overwritten. The receipt retains code/input hashes,
counts and limits. A version-6 consumer dataset without its final successful build
receipt is rejected. No vendor or original supplementary workbook is needed to
rebuild from these pinned components.

`dataset.source_evidence().summary_group(summary_id)` returns the pinned baseline,
its reported representative, all matching native protein records and their original
group membership links. Missing native records and missing/unresolved links remain
explicit. These associations neither make the summary an independent experiment
nor turn a shared group into unique protein abundance. Every summary row is checked
for sample identity, duplicates, original iBAQ/rIBAQ missing/zero/positive counts
and absent matrix cells; abundance values remain in the baseline dataset.

The completed source-group audit reconciled all 15,726 membership rows against
12,654 native groups: every declared `ProteinCount` equals the membership count.
All 9,929 summary representatives join one native group each; 1,052 of those
groups have multiple members. This closes the group/summary reconciliation only;
the remaining peptide/PSM extraction and full artifact acceptance are still required.

After assembly, audit the consumer artifact against the original selected native
SQLite export, full searched FASTA and protein summary:

```sh
python -m tools.verify_dog_dataset /path/to/built/dataset.sqlite \
  /path/to/audit/pdResult/native-evidence.partial.sqlite \
  /path/to/lineage.json /path/to/search.fasta /path/to/protein.tsv /path/to/acceptance.json
```

This read-only command requires the full pinned selection. It compares **every
selected native field** with the artifact, including original BLOB/float bytes,
signed record IDs, confidence/modification definitions and target/decoy partitions.
It checks identity/sequence indexes, all reported membership endpoints (including
missing endpoints), native protein/search identity, all quantitative projection
targets/units/sample axes, original missing-acquisition reports and logical counts.
All common sample/assay records must match the pinned compiled lineage input.
It also verifies every retained searched protein, candidate-index coverage and
digest, summary representatives/native associations, original summary counts and
native declared group sizes versus membership counts. Source-declared group-count
inconsistencies are reported for review, without erasing the source facts.

Public measurement expansion and independent exhaustive substring searches of the
complete FASTA use a deterministic spread of at most 100 records per quantitative
table by default, including the first and last records. `--sample-records` changes
that explicitly (1–1,000). The receipt states full-scan versus sampled scope,
source/artifact hashes, compared counts, peak RSS and tool hash. Inputs are checked
for changes before completion. The command acquires no vendor data and writes only
its new receipt. The audit does not claim a checksum of unacquired vendor bytes.
Sampled numeric values and missingness also undergo independent `struct` decoding
of the original number/presence bytes, without the builder's NumericVector decoder.

Prepare the compressed archive and registry metadata only after that audit passes:

```sh
python -m tools.package_dog_dataset /path/to/built/dataset.sqlite \
  /path/to/acceptance.json /path/to/new-publication-directory \
  --url https://example.org/REVIEWED-RELEASE/dataset.sqlite.gz \
  --data-release REVIEWED-RELEASE --max-compressed-mib 1024
```

The URL above is a placeholder for the actual reviewed distribution destination.
The command does not publish or update the packaged registry. It requires the
full source selection and a matching successful receipt from the current verifier.
It writes a gzip archive with fixed header metadata, checks every decompressed byte
against the original hash without creating a second expanded file, and checks
for changed inputs. Compressed output has an explicit byte cap. Failure removes
only its owned new archive; existing archives are never overwritten.
`registry-entry.json` records the installed checksum/size used by DataCache and
the separate compressed checksum/size, data-release identity and source access
terms. The full searched FASTA's UniProt attribution remains explicit. The
archive, source components, acceptance receipt and registry entry still require
publication and real-release installation/consumer acceptance before #23 is done.

Offline tests install a compressed version-6 artifact through DatasetRegistry,
remove its private build inputs, and replay measurements and native group/search
queries without acquiring those inputs or the separately cached FASTA.

Native evidence encoding version 2 stores scalar values in declared column order
with an explicit presence mask. Column names remain in the per-table schema;
null scalar values remain different from omitted binary/vector fields. The
reader also accepts the original version-1 object encoding. All public native
records remain dictionaries, and bulk verification uses the original typed-row
digests. No selected fields are discarded to reduce storage.

With the production writer, the same 10,000-PSM slice and its native metadata
occupied 22,323,200 bytes as version 1 and 7,340,032 bytes as version 2. A
10,000-group spread occupied 28,856,320 and 23,498,752 bytes respectively. Every
selected native field, including numeric and opaque bytes, reconstructed to the
original typed digests. The rebuilt group measurement artifact retained all
1,320,000 logical cells and the identical full-field digest of 30,831 measured
records. These are bounded source slices, not full-artifact acceptance or size
guarantees.

`NumericVector` stores declared int32/int64/float64 values in little-endian bytes
with separate presence bits. A null vector, a missing element, a measured zero and
a nonfinite source value remain distinct. IEEE-754 payloads and signed integer IDs
round-trip exactly. `NativeDefinitions` uses each source column's own enum domain;
it never interprets protein confidence with a peptide confidence enum, or converts
a native modification position into a residue offset without defined semantics.

`EvidenceStore(path).records(table, key=..., sequence=...)` reads the normalized
component locally. It returns the original selected field names/values; `vector`
provides their normalized numeric encoding. No APSW, remote range access or vendor
database is required for these reads.

The version-6 quantitative projection stores one index entry per abundance vector
and shared sample/assay axes. Public `measurements`, mapped/paired queries and
Parquet export reconstruct ordinary `Measurement` records on demand, including
missing cells. Counts describe all logical cells. Native PSM abundance and the
three peptide-group abundance columns keep separate units and record IDs. Source
protein membership and searched-FASTA candidates remain separate attribution facts.
The internal flat observation tables are absent in version 6, so direct SQL cannot
misreport an empty observation set. Ordinary version-5 dataset writes are unchanged.

After `build_attribution` has completed against a checksum-verified complete
searched FASTA, `dataset.source_evidence()` exposes the indexed source facts:

```python
native = dataset.source_evidence()
group = list(native.records("TargetPeptideGroups", key=(8,)))
memberships = list(native.linked("TargetPeptideGroups", (8,),
                               via="TargetPeptideGroupsTargetProteins"))
candidates = list(native.candidates("AAAAAAAAAK"))
```

Each link retains its exact original fields and both endpoints, including an
explicit unresolved endpoint. Incoming queries retrieve group/peptide records
mentioning a native protein. A resolved link means that both source records are
present; confidence and exclusion policies are not inferred. The searched-reference
index retains every matching accession and repeated position, labels exact versus
I/L-compatible matches, and distinguishes no match from a sequence never searched.
Only target peptide sequences are searched by default. Decoy source records and
links remain accessible with their native partition; they are not assigned as
biological detections. All original searched protein entries are available via
`native.search_protein(accession)`, separate from any detection claim.

The complete target peptide membership audit contains 1,145,689 peptide groups,
1,262,557 peptide/protein links and 170,084 protein-group/peptide links. Every link
has both endpoints, and every peptide/protein pair has compatible sequence:
1,258,153 exact-residue matches and 4,404 I/L-equivalent matches. These are source
relationships, including excluded records, not a filtered set of detections.

Displayed counts do not necessarily count link rows. Across all peptide groups,
`NumberOfProteinGroups` and `MasterProteinAccessions` agree with the distinct
linked master proteins (`IsMasterProtein=0` in native enum 14); 458 peptide groups
have more group-link rows than distinct masters. Across all 12,654 protein groups,
`PeptideGroupCount` agrees with distinct linked peptide sequences after treating
I/L as equivalent, also with distinct native `UniqueSequenceID` values. Raw link
counts differ for 5,190 groups; counting distinct literal sequences leaves 18
differences. These are measured equalities for this source, not inferred rules
that replace the deposited counts or collapse the original peptide sequences.

Four peptides have a native parent link absent from their displayed parent lists:
160079 (`DLYTGLK`), 516803 (`IGDFGLAR`), 697816 (`MPLIGLGTWK`) and 993593
(`TGTFIAIDR`). In each case the additional protein has `ExcludedBy=3`, while the
listed proteins have `ExcludedBy=-1`. The `parent-status.json` regression fixture
retains all four cases, their complete native parent links and full-reference
sequence candidates. Neither the link index nor the candidate index rewrites
the displayed count/list. This audit covers the completed target memberships;
full PSM, control, acquisition and consumer-artifact acceptance remain separate.

The completed target PSM table contains 4,493,585 records and 953,980 distinct
sequences. All 952,678 distinct peptide-group sequences occur among them, but
7,575 PSMs supply 1,302 additional residue strings. Each has an I/L-equivalent
peptide-group sequence; none is replaced by that group's literal sequence.
Independent exhaustive searching of all 1,302 additions against the complete
FASTA found 9,380 placements. The combined target sequence audit therefore has
1,171,840 placements: 1,132,603 exact and 39,237 I/L-equivalent, with no unmatched
sequence. Candidate indexing must include PSM-only sequences and retain the
match basis for each original string. Complete native PSM links and quantities
are verified separately during full artifact acceptance.

All 4,493,585 native PSM/peptide-group links have both endpoints; each PSM links
to exactly one group, with no duplicate edges. Their native `UniqueSequenceID`
values agree on every link. Literal sequences agree on 4,485,319 links; the other
8,266 are I/L-equivalent and retain both original strings.

All 5,289,013 target PSM/protein links have both endpoints and exact residue
matches, with no duplicate edges; every target PSM has at least one parent.
Their parent sets differ from the linked peptide groups: 114 PSM/protein edges
across 103 PSMs are absent from their group's parent set, while 45,450 implied
group-parent edges across 17,221 PSMs are absent from the native PSM links.
The latter includes identical PSM/group sequence strings and high-confidence
PSMs. These measured differences do not establish a vendor assignment rule.
Consumers must preserve each native relationship rather than copy a group's
parents to its PSMs or replace a group's parents with an individual PSM's set.
The `psm-parent-sets.json` fixture retains three complete examples, including
PSM 15032 (`LTEEQAR`) and group 606096 (`ITEEQAR`), whose single native parents
are different proteins. Independently enumerated full-FASTA candidates retain
both proteins and their match bases separately from those native assignments.

The source's displayed `PsmCount` is distinct from an unfiltered link count.
Across all 1,145,689 peptide groups, every displayed count agrees with this
comparison: for native group `Confidence=3`, count linked PSMs with
`MatchConfidence=3`; for other groups, count all links. Thus 131,596 displayed
counts differ from raw link counts. This is a measured equality in this deposit,
not a default evidence filter or a general rule for other vendor analyses.
All 962,755 protein-group/PSM links, by comparison, have `MatchConfidence=3`
and `ExcludedBy=-1`; their counts agree with every one of the 12,654 groups'
displayed `PsmCount` values. Lower-confidence source PSMs and their other links
remain available unchanged.

The complete modification tables contain 10,063,345 target links and 7,644,492
decoy links. Every link resolves to a native PSM and one of the 1,787 retained modification definitions;
there are no duplicate `(workflow, peptide, modification, Position)` tuples,
and every target and decoy PSM has at least one modification link. The same five
definition IDs occur in both partitions:

| Native ID | Name | Native specificity label | Target links | Decoy links |
| --- | --- | --- | ---: | ---: |
| 2 | Acetyl | Protein N-Terminus | 57,330 | 44,868 |
| 8 | Carbamidomethyl | Any | 1,437,923 | 1,274,637 |
| 45 | Oxidation | Any | 1,372,286 | 1,036,298 |
| 626 | TMT6plex | Any | 2,702,221 | 1,968,300 |
| 628 | TMT6plex | Any N-Terminus | 4,493,585 | 3,320,389 |

Specificity labels come from native `PositionType` enum 3. IDs 626 and 628
share a name and Unimod accession but remain distinct definitions. All native
`Position` values in both partitions are present; their integers are retained
without converting them to residue offsets or inferring localization certainty.
These are unfiltered source-link counts, not independent detections. Full
consumer-artifact verification remains separate.

The completed control tables contain 20,596 decoy proteins, 1,521 decoy protein
groups and 1,076,028 decoy peptide groups. The peptide groups' native confidence
labels are Low for 1,065,181, Medium for 9,126 and High for 1,721 records. All
1,721 High records and all 1,521 decoy protein groups have `ExcludedBy=-1`.
They retain their decoy partition alongside those source quality fields.
Decoy protein `ProteinFDRConfidence` uses distribution map 36 with one `Combined`
cell: 18,857 null vectors and 1,037/574/128 vectors containing native codes
1/2/3. These codes retain their own distribution metadata; peptide confidence
enum 12 does not label them. The complete typed row digests and declared
protein/peptide-group sequence lengths reconcile.

All 3,320,389 decoy PSMs also reproduce their complete typed source digest.
Their native confidence labels are Low for 3,268,321, Medium for 43,681 and
High for 8,387 records. Of the High records, 3,992 have `ExcludedBy=-1` and
4,395 have `ExcludedBy=3`; all Low and Medium records have `ExcludedBy=3`.
Every decoy PSM has native `PsmAmbiguity=5` (`Unambiguous`, enum 11). These
quality labels do not change the control partition. The selected PSM schema
has no declared sequence-length field, so no such length comparison is claimed.

The six complete decoy membership tables have no missing endpoints or duplicate
edges. All 3,637,761 PSM/protein links match exact residues; the 1,179,103
peptide-group/protein links comprise 1,175,569 exact and 3,534 I/L-equivalent
matches. Each PSM links to exactly one peptide group with the same native
sequence ID; 3,880 links have different I/L-equivalent literal strings. All
original strings and relationships remain separate. There are 171,739 PSMs
with multiple native protein parents, despite the `Unambiguous` PSM label;
the largest observed parent count is 256. This is not evidence of a vendor cap.
All 1,521 protein groups' declared counts reconcile with distinct linked
proteins, native peptide sequence IDs and PSMs. Their 4,003 PSM links all have
native `MatchConfidence=3` and `ExcludedBy=-1`. Acquisition records and full
consumer-artifact acceptance remain separate audits.

On an explicit 100-group benchmark, the public source-evidence interface preserved
all 13,200 quantitative cells and all 20,649 searched proteins. All 120 candidate
placements agreed with an independent exhaustive search of the complete FASTA.
This benchmark contained no full native membership tables and is not a release
or a substitute for their full-count and representative reconciliation.

Sequential scans use `--read-ahead-kib`; index lookups default to individual
4 KiB pages. `--lookup-read-ahead-kib` can independently increase the indexed
window when measured page locality warrants it. Both options accept multiples
of 4 KiB up to 16 MiB, stop before existing cached segments, and retain the same
cumulative transfer/cache limits and exact record selection. Receipts record
both settings. The source control connection can be reused only after a verified FTP
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
