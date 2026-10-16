# RefSeq protein snapshots

Add `Proteome.from_refseq(protein_fasta, feature_table, reference=..., expected_sha256=...)`.
Read verified bytes offline, require matching assembly/provider/release, and join
versioned protein and RNA accessions to NCBI Gene IDs from the NCBI download
feature table. Retain all occurrences and the complete FASTA inventory. Missing
or contradictory mappings fail. Report CDS rows without proteins separately.
Do not infer spliced coordinates from genomic spans. Persist both source hashes
in the annotated layer and pin them in ReferenceContext. Verify small malformed
fixtures, multi-locus proteins, deterministic replay and the real dog annotation.
