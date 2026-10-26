# Canine peptide and protein-group extraction (#23)

Audit the public processed-file inventory before choosing an input. The smaller
pdResultView and msfView databases may expose the required facts without the
86 GB pdResult. Query their real schemas and records through SQLite, not a
second parser inferred from the prefix.

Use an opt-in, read-only FTPS range reader and APSW VFS. Bound transfer, cache,
requests, elapsed time and SQLite allocation; persist source identity and chunk
hashes, count failed/retried payloads, and refuse a whole-file fallback. Verify
REST through independently requested overlapping bytes. This detects ignored or
inconsistent offsets, but is not an authoritative complete-file checksum.
Record this integrity limit explicitly. Reject identity changes and corrupted
cache entries, with no mixing of remote versions.

First extract the schema and a bounded set of native records. Choose columns and
joins only from that evidence. Preserve peptide/PSM/modification/confidence and
control definitions, all reported group memberships, run/plex/channel/bridge
and pooled-tissue lineage. Resolve against the exact deposited searched FASTA,
without claiming it covers all isoforms. Retain missing values and the unresolved
publication/deposit instrument reports. Reconcile selected original records and
source summary counts before producing the compact dataset artifact.

Offline tests must exercise real SQLite index/interior/overflow traversal,
interrupted and wrong-offset reads, source changes, budget exhaustion, cache
corruption, shared mappings and missing quantities. Normal library imports and
dataset queries must not initiate vendor acquisition.
