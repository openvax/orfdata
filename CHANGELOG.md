# Unreleased

- `orfdata data fetch` repairs an installation whose files fail verification,
  reporting `acquired: repaired`, as `setup` already did; `--force` still
  replaces valid files. An unreadable installation is reported, not replaced,
  unless `--force` is given, in both `fetch` and `setup`, which now share one
  policy; a multi-name fetch checks this for every name before downloading (#152).
- `orfdata data status` lists damaged and unreadable installations, each with its
  exact `data fetch` command, instead of hiding them. Scripts that read
  `status --json` as "usable only" should filter on `installed` or use
  `data list --installed`. `list_data(local=True)` selects the same rows. An
  unreadable viral collection is `inaccessible`, not `corrupt`. Fetch commands
  in listings carry `--cache-root` when one is given (#151).
- `orfdata data list` and `status` share their 25-row preview among resource
  kinds, so datasets, genomes, search references and viral collections appear
  next to catalogues and chains, and report what each kind left out (#153).
- `orfdata data fetch` summaries name the resource, its version and what the
  command did (`downloaded`, `built`, `already installed`) for every kind,
  including viral collections (#154).
- Exit with a status scripts can branch on: 1 when a command fails (network,
  checksums, damaged data or files, unknown names), 2 for a command line the
  parser or the command's own checks reject, and 3 when a resource the command
  needs isn't installed. Those checks run before anything is read; a few option
  values are checked only by library code and exit 1. Only command-line
  mistakes point to `--help`. Up to 0.44.0 nearly every failure exited 2; a partial
  `setup` exited 1 (#155).
- Raise `orfdata.NotInstalledError`, also a `FileNotFoundError` and a
  `ValueError`, for catalogues, datasets, search references, genomes, chains,
  viral collections and PyEnsembl annotations that aren't installed; opening one
  now names the command that installs it, with its kind, version and cache root,
  instead of only its cache path.
- Report unexpected internal errors with their traceback instead of as a
  one-line user error. A malformed catalogue file is reported by name. An
  unknown dataset `--reference` or `--asset` lists the choices, and a restricted
  chain without `--accept-terms` is a usage error in every command.

# 0.44.0

- Show a command guide for bare `orfdata` and command groups; add `help COMMAND`
  and an explicit `query` subcommand with discoverable examples.
- Default to concise human summaries and bounded, wrapped tables across every
  command. Preserve native measurement units, reference scopes and pairing basis.
- Require explicit machine output: `--json` for complete JSON documents and
  dataset `--format jsonl` for the former streaming representation. Explicit
  FASTA/TSV/GTF formats and Python APIs are unchanged. Human `--limit` settings
  never truncate exports. Proteome/lift summaries remain on stderr.
- Run every CLI leaf in human and JSON modes against tiny offline fixtures in CI,
  with exhaustive help, error and output-contract checks. See the
  [CLI migration guide](https://openvax.github.io/orfdata/cli/#migrating-scripts-from-042-or-043).

# 0.43.1

- Use matching Ensembl chromosome DNA in the CHM13 installation example and
  explain the consortium FASTA's different contig names (#149).

# 0.43.0

- Resolve pinned annotations' taxonomy and genetic codes from their declared species,
  including T2T-CHM13. Preserve explicitly saved reference metadata on replay.
- Add explicit chromosome FASTA installation/loading through PyEnsembl to the
  pinned-genome Python API and `orfdata genomes install --genome-fasta`. Loading
  remains offline, and annotation-only installations remain supported (#146).
- Require PyEnsembl 2.32.3 for native DNA and saved genetic-code metadata.

# 0.42.4

- Reconstruct mouse RibORF source ORFs from pinned GENCODE M25 RNA and exon
  paths before placing them on Ensembl 107. Require exact full RNA, CDS and
  protein preservation; retain changed or ambiguous native records explicitly.
- Publish the verified 24.9 MB generation-3 catalogue through the common data
  commands; rebuilding from 198 MB of source references remains available.
- Keep native reference/evidence identity through typed `ReferenceMapping`
  provenance on placed ORF records, including JSON and TSV replay. Generation 3
  replaces unvalidated mouse generations 1/2 for new builds.

# 0.42.3

- Import canine mammary-stroma RNA, protein-group and peptide measurements from
  GSE135183/PXD023023 through the common data commands and dataset API. Preserve
  14 shared specimens, 28 paired compartments and two missing protein mates,
  native units, controls, group ambiguity and source-bound incomplete references.
- Add explicit `same_specimen` pairing for separate preparations without asserting
  identical aliquots or adjacent sections. Preserve MS zero/NaN missingness and
  measured RNA zeros; do not infer gene/protein equivalence or ncORF translation.

# 0.42.2

- Publish the existing 16,484-record mouse ncORF catalogue of Chang et al.
  as a checksum-pinned download. `orfdata data fetch ncorf-mammals-mouse`
  downloads the 2.22 MB artifact; explicit source builds remain available.
  Source/reference provenance, aggregate Ribo-seq evidence and the unresolved
  source row are unchanged.

# 0.42.1

- Include the actual carrier RNA's contig and strand in every `TranscriptORF`
  export, including annotated CDSs and `at` placements without a catalogue
  model. Preserve unknown source-RNA locus values as `null` (#135).
- Keep older v1 snapshots readable while requiring any supplied carrier locus
  fields to match during replay. Preserve the original source model separately
  when its ORF is placed on another transcript.

# 0.42.0

- Add `orfdata data info/fetch/path` for resources named in `data list`;
  `downloads` and `download` address the same command group. Fetch automatically
  chooses a published artifact or a build from pinned sources, supports multiple
  explicit names, and checks all names and restricted terms before acquiring.
- Add matching `data_info()`, `fetch_data()` and `data_path()` Python functions.
  Preserve version selection, cache roots, checksums, existing generations and
  explicit downloads. Bundled catalogues now expose their installed file path.
- Document verified processed canine paired-study sources and the outstanding
  dog Ribo-seq/ncORF evidence gap. No new experimental dataset or inferred
  sample pairing is silently added to the registry.

# 0.41.0

- Include PyEnsembl, spreadsheet readers, Parquet exports and sequence matching
  in normal installations. Remove the redundant `pyensembl`, `build-data`,
  `columnar`, `data` and `viral` extras and their obsolete runtime install
  branches. Reference and measurement downloads remain explicit and pinned.
- Use Biotite for the experimental dog search's BLOSUM62 local alignments,
  preserving scores, affine gaps, free local overhangs and optimal-alignment
  tie order. Refuse incomplete tie enumeration. Biopython remains a development
  dependency for independent comparisons; evidence states are unchanged.
- Retain `canine-extraction` only for APSW's vendor-database VFS. Published
  dataset queries and normalized-source rebuilds work with a normal install.
  Python 3.11 uses compatible Biotite 1.6; Python 3.12+ supports Biotite 2.

# 0.40.0

- Add bounded, resumable FTPS extraction and lossless normalization of Pfizer
  canine peptide, PSM, protein-group and acquisition evidence. Preserve source
  identity checks, cumulative transfer accounting and explicit resource limits
  across interrupted runs (#23).
- Add native evidence queries, reported memberships, separate searched-FASTA
  candidates and links to the existing same-study protein summary. Preserve
  confidence/modification definitions, decoys, missingness, ambiguity and
  tissue-pool/run/channel provenance. A compact quantitative driver retains
  native numeric vectors and exposes their cells through the measurement API.
- Expose reported assignments and sequence-compatible candidates through the
  common `mappings()` API and protein-group memberships through `relations()`.
  Preserve all candidates when filtering by a destination; incomplete source
  assignments remain explicitly `unresolved`. Source-specific details stay in
  the lower-level adapter.
- Support pinned split tar.gz source components and dataset artifacts through
  DataCache >=1.29.3, with verified logical-file installation, offline inspection
  and reuse, and preservation of existing cached generations on failed refresh.
- Publish the complete canine consumer dataset: 70,280,798 selected source rows
  across 42 tables, 953,980 target peptide sequences and 1,171,840 compatible
  placements. Full-source acceptance verifies the native facts, common queries
  and frozen independent candidate digest. Download 6.90 GB; installed size
  24.71 GB. The bounded build reserves 34.4 GB rather than 103.1 GB of free space.

# 0.39.1

- Allow `ScopedTarget.from_orf` to identify annotated and source-defined
  transcript ORFs without requiring a catalogue record (#61).
- Link measured source regions to candidate or partial ORFs with an explicit
  `region_orf_context` relation. Retain region measurements, mapping evidence
  and candidate ambiguity; the relation does not establish ORF translation,
  transfer abundance or authorize ORF-specific ratios (#64).

# 0.39.0

- Keep bounded source ORFs on their full parent RNA with explicit frame, separate
  unknown/incomplete start and stop boundaries, and leading/trailing partial
  codons retained in coding coordinates. Protein sequences contain translated
  complete codons only; no missing sequence is padded (#64).
- Expose coding origin/phase, coding genomic positions and reference identity
  directly on TranscriptORF. Coordinate queries retain unavailable residues for
  partial codons and unknown frames. Catalogue provenance remains optional for
  annotated CDSs and source-defined placements (#61).
- Add checked TranscriptORF reconstruction against a supplied or archived RNA,
  preserving source evidence, reference context and exact sequences. Replay
  source-defined translation to reject inconsistent phases or proteins.
- Pin native RPL18 and KAZN partial-CDS source rows and verify both strands,
  non-AUG starts, exact boundaries and annotation/source snapshot round trips.

# 0.38.2

- Accept explicit source-transcript genetic-code declarations in OpenProt
  import manifests. Mitochondrial rows without a declaration still fail;
  unsupported codes and declarations naming absent transcripts are rejected.
  Verify the pinned human IP_306399 row under table 2 while preserving its
  coordinates and aggregate evidence counts (#57). Whole human/yeast catalogue
  publication remains separate from this qualified import path (#18).

# 0.38.1

- Accept native ORF IDs reused by distinct sources in directly constructed
  catalogues. Qualify colliding catalogue IDs automatically while preserving
  each original record and its evidence in `members`. Reject duplicate IDs
  within one source and raise on ambiguous bare or combined-string lookups.
- Add explicit `source=` selection to `ORFCatalog.get()` and
  `ExtraORFs.for_orf()`. Preserve colliding occurrences through binding,
  proteome construction and offline round trips (#71).

# 0.38.0

- Add an experimental, resumable Sage pipeline for re-searching the 132
  available Canine PeptideAtlas mzML runs against Ensembl 110 dog proteins and
  candidate ncORFs. It verifies remote Deflate64 members, preserves per-run
  receipts, applies class-specific target-decoy FDR and conservative sequence
  filters, and exports candidate peptide and ORF tables without treating them
  as released ORFdata evidence. One real PXD052932 run now covers the complete
  fetch, Sage, parquet, aggregation and restart path.
- Score ties receive one q-value for the complete tie; Ensembl selenocysteine
  and other non-BLOSUM62 residues are normalized explicitly for homology; and
  an incomplete batch exits nonzero. Regression tests cover all three cases.

# 0.37.0

- RibORF catalogue release 2 preserves every short source ORF that passes the
  same transcript, sequence, start and stop checks as longer records. It adds
  11,902 mouse, 11,679 human, 5,732 zebrafish, 1,186 C. elegans and 241 yeast
  ORFs shorter than five residues; failed rows remain individually traceable.
  Release 1 remains reproducible and addressable by version.
- ORF length is now an explicit consumer policy. `ORFCatalog.query()`,
  `extra_orfs()`, `with_extra_orfs()`, `with_nc_orfs()`, `Proteome.build()` and
  `ReferenceTranscriptome.proteome()` accept `min_protein_length`; both query
  and proteome CLIs expose `--min-protein-length`. Proteome layers and cache
  identities record the threshold and the number filtered, while source
  catalogues stay lossless.

# 0.36.0

- Public names now say what their values are. `ProteinSequence` is one protein
  sequence with provenance, `ProteinSet` is a named set of protein sequences
  from any source, and `list_data()` returns `DataStatus` values for bundled,
  downloadable and buildable data. Fetcher and viral-cache implementation
  classes are no longer exported from the package root. A `Proteome` exposes
  `.proteins`, `from_proteins()` and `with_proteins()`; producer links target a
  `protein` or `record`, so the public API no longer calls proteins “entries.”
  Downloaded and selected viral files share the `ViralReferenceFiles` return
  type and the literal verbs `open_viral_references()` and
  `select_viral_references()` replace “bundle” and “index” in the public API.
- `ProteinSet` is not virus-specific. `with_protein_set()` adds an explicitly
  selected donor, transgene or other component while retaining its source
  taxon, exact reference, checksum, role and sourced target assertion.
  `with_viral()` accepts the same type plus a sourced host association; curated
  `ViralReference` objects continue to provide that association themselves.
- Producer identity no longer includes the wrapper that exported it. Direct
  and wrapped imports of the same source protein or evidence record merge once
  and retain every path in `exporters`; distinct samples and conflicting
  biological content remain distinct or fail. Supplied proteins may keep an
  unknown source taxon, which remains `None` in snapshots and `unknown` in
  FASTA and requires `compatibility="allow_unknown"` before admission to a
  known host. FASTA headers now include producer version, provider, export
  paths, exact reference, source checksum and explicit true/false/unknown
  boundary states with injection-safe escaping. Exact-sequence deduplication
  writes a neutral group header and keeps every contributing protein in the
  FASTA index instead of borrowing one producer's metadata.
- Viral collection selection is one API everywhere. `viral_collection()` now
  accepts a collection object, name, alias, panel name, or a species with one
  registered collection; the download/open functions derive the target species
  from it instead of defaulting independently to human. `viral_collections()`
  can filter by species or source without reading or downloading data. A
  collection's `source()` and `source_summary()` both default to its declared primary
  source, and an undeclared source is an error rather than an empty answer.
- The checked-TSV descriptor, expected-sequence check and fetched-genome value
  are internal implementation types rather than package-root API. Users work
  with collections, ORF catalogues, quantitative datasets and proteomes; file
  layout and verification do not add another public vocabulary.
- Reference tables are pandas data frames. `viral_panel_coverage()` and
  `viral_contexts()` return a `DataFrame` instead of a tuple of dataclasses,
  and `collection.source_summary(source)` returns one keyed by entry abbreviation, so
  an inventory or a source summary can be filtered, joined and aggregated with pandas
  rather than through bespoke selectors. The packaged table's own header lines
  travel on `frame.attrs`. `ViralPanelEntry`, `LedgerEntry` and `ViralContext`
  are gone: a row of a frame is a row of a frame. A source summary's required
  sequence is checked internally, and the per-row helpers are plain functions —
  `published(row)`, `same_molecule(row, sequence)` and
  `verify_record(row, record)`. Source-summary construction also returns a frame. pandas 3 is a new
  runtime dependency. On-disk formats and checksums are unchanged: the
  packaged tables round-trip byte for byte through the new reader and writer.
- One way to carry checked metadata, whatever its source or species.
  Every TSV ORFdata ships — a collection's inventory and one source summary
  per source — is read through the same internal `CheckedTSV`: a tab-separated file with a named
  schema and pinned SHA-256, read and written in one place
  (`orfdata.reference_table`) instead of once per file. A collection declares
  the tables it reads rather than naming resources and checksums in module
  constants, so a second panel, species or source is a registry entry and not
  another loader.
- One source-summary format for every source. The existing on-disk
  `orfdata.reference_ledger.v2` schema records,
  per entry and per source, the identifier that source resolved it to, the
  molecule's length, MD5 and SHA-256, how many coding features it annotates,
  its NCBI taxon and its own description. The GenBank side of the GDC
  collection, which used to be a separate summary table in its own format, is
  now `gdc-grch38-d1-vd1-genbank.tsv`, and the alternate source summary is
  `gdc-grch38-d1-vd1-pave.tsv`; `collection.source_summary(source)` reads either.
  `tools/build_reference_ledger.py --source NAME` builds any of them from one
  pass over the inventory, replacing the alternate-only builder.
- Check every downloaded record against the pinned summary for its source, not only
  the alternate ones. A primary record whose length, checksums or coding
  feature count no longer match what the packaged summary pins now fails the
  download loudly instead of being written. The record's source and identifier
  must match the summary too, and an expected sequence must occur in the summary for the
  source it names; identical bytes from the wrong source cannot satisfy either
  check. The GenBank summary's checksum is also the collection's cache identity,
  so a re-checked source invalidates caches by itself.
- Tell a source that does not publish a record from one that cannot be read.
  `SourceUnavailable` is raised for a refusal the server states, so a source summary
  records `unavailable` and `inconsistent` from what happened rather than from
  matching an error message.
- Read a GenBank record's `DEFINITION` and its `source` feature's taxon, so a
  fetched genome carries its own title and NCBI taxon and no second lookup is
  needed to describe it. The parsed document remains a `GenBankRecord`; the
  fetchers share the same internal genome value.
- Read an IUPAC ambiguity codon as the residue it determines. A codon carrying
  an ambiguity code translates to a residue whenever every base it stands for
  gives the same one — `ATY` is isoleucine either way — and only to `X` when
  they genuinely differ, which is how sequence databases publish them. This
  corrects the two source-summary rows previously recorded as `inconsistent`:
  HPV171 and HPV172 each carry exactly one `Y` in their published genome, in
  an `ATY` and a `TGY` codon that cannot change the protein, so PaVE's
  coordinates do reproduce its published proteins after all. The alternate
  source is now `published` for 187 of the 189 papillomavirus rows, the other
  two being types it does not publish at all.
- Record how each coding sequence begins. `transcripts.tsv` carries
  `start_codon`, the three bases read first, and `initiation`, one of `atg`,
  `alternative` (another initiation codon of the same genetic code),
  `declared` (the source names the initiator with `transl_except`),
  `non_canonical`, `partial` or `unknown`; `viral_reference_summary()` counts
  them and names the non-canonical ones. Across the GDC collection that is
  `atg` 1,686, `partial` 66, `alternative` 3, `non_canonical` 2 and
  `declared` 1 — the last CMV's `UL30A`, an `ACG` start with
  `/transl_except=(pos:complement(38131..38133),aa:Met)` and ribosome-profiling
  evidence, which ORFdata previously dropped. The two non-canonical ones are
  HPV63 and HPV65, whose submitters note that the protein conflicts with the
  conceptual translation. The residue a start yields stays in the proteome,
  where the source publishes it.
- Record translational recoding the source claims. A valueless GenBank
  qualifier is no longer discarded, so `ribosomal_slippage` reaches a new
  `recoding` column and the summary's `recoded_transcripts`. Two transcripts
  carry it: HTLV-1's `gag-pro` and `gag-pro-pol`, written as joins that read
  the slippage base twice. Their sequences and proteins are correct and
  varcode reads them per product, but a variant at the reused base is reported
  as an exon loss, because a GTF exon is a genomic interval and effect
  predictors assume exons do not overlap.
- Survey how completely a collection is annotated. `viral_annotation_survey()`
  and `orfdata viral completeness` report, per genome, its coding transcripts,
  how many are complete, how many carry no gene name, how many are spliced,
  and — where a collection names an alternate source publishing the identical
  molecule — how many more coding features that source publishes, read from
  the packaged source summary without a network request. No gene list is built in, so
  the survey is true of any collection, species and source.
- `--alternate-annotations` (`prefer_alternate_annotation=True`) takes the
  alternate source's coding features wherever it publishes the identical
  molecule and annotates at least as much of it, in place of each submission's
  own. The genome FASTA is byte-identical either way: this changes only how
  the same bases are read. Over the GDC collection it is the difference
  between 1,758 and 1,993 coding transcripts, between 40 and 194 genomes with
  a spliced model, between 49 and 3 with a transcript whose reading frame does
  not close, and between 216 and 14 papillomavirus coding features with no
  gene name at all — at the cost of replacing each submission's naming with
  one source's. The remaining gaps are the two papillomavirus types PaVE does
  not publish, and the herpes, hepatitis, retrovirus and polyomavirus rows,
  which have no alternate source. The default still changes nothing it does
  not have to.

- Add the cached GDC `GRCh38.d1.vd1` source-reference collection, named
  `gdc-viruses`, as three layers that are named the same way everywhere:
  `genome` (200 viral genomic DNA records, the molecules a DNA reference such
  as GRCh38.d1.vd1 carries as extra contigs), `transcriptome` (one spliced
  coding transcript per source CDS feature, with `annotation.gtf` placing it on
  the genome) and `proteome` (the source translation of each transcript).
  `open_viral_references(download_if_missing=True, species="human")` acquires
  it; `bundle.layers`, `bundle.layer("transcriptome")` and
  `viral_reference_summary()` describe it, and `orfdata viral summary`
  prints that summary. Each cached GenBank record is stored exactly as NCBI
  sent it, so `SHA256SUMS` attests to the source record itself.
- The GTF, transcript FASTA and protein FASTA of a collection share
  identifiers, so `viral_genome(bundle_or_index)` opens it as a PyEnsembl
  genome and varcode, isovar, topiary and vaxrank can give a variant on a viral
  contig a coding consequence, including under a non-standard genetic code the
  source record declares (`transl_table`). A frameshifted CDS, or one the
  source marks as running past the record, carries no start or stop codon rows;
  `transcripts.tsv` records that as `annotation_complete` and `partial`
  (`5_prime`, `3_prime` or `both`), which covers 69 of the 1,758 transcripts:
  EBV's repeat-region genes, HIV-1's `pol`, which is translated as a Gag-Pol
  frameshift fusion, and the E4 and E5 ORFs of papillomavirus rows still on
  their inventory accession, which carry no start codon of their own because
  they are expressed from spliced transcripts.
- Carry the molecule the reference a collection follows was actually built
  from. A reference inventory can name an entry at two sources whose sequences
  differ, and GDC's decoy FASTA carries the PaVE reference rather than the
  GenBank accession beside it: PaVE's `HPV16REF` reproduces the GDC HPV16 DNA
  MD5 this repository already pinned (`f2463721818c29998ef5ab917e6df74a`)
  exactly, while `K02718.1` is a different sequence two bases shorter with the
  sparser 1985 annotation. 24 of the 189 papillomavirus rows differ this way,
  among them HPV16, HPV18, HPV56 and HPV59 — four IARC Group 1 agents — and
  the inventory's own `PaVE Notes` column marks 18 of them `Revised`.
  HPV16 now carries 7,906 bases and ten complete coding transcripts, including
  the spliced E1^E4, E8^E2 and E6\* and an 83-residue E5, where the GenBank
  accession gave eight with E4 and E5 5'-partial and no spliced models.
- Sources are a named, published thing rather than a special case. The internal
  source registry knows `genbank` and `pave` today, and each fetch returns the
  molecule, its coding features and the source's own document. A collection stores every document
  under `records/`, named for the entry, the source and the identifier, so
  `SHA256SUMS` attests to whichever source served it, and `genome_source` and
  `annotation_source` say where each genome and transcript came from.
- Collections are a registry rather than constants. `viral_collections()` and
  `viral_collection(name)` describe what a collection is for — its target
  species and taxon, inventory, contexts, primary source and alternate — and
  every entry point takes a collection name, so a second species or source set
  is a registry entry rather than a code path. `install_species` reads the
  collections a species has from the same registry.
- Source summaries are general: the internal reader handles a packaged,
  checksummed record of what each entry looks like at any collection source
  source, with states `published`, `unavailable` and `inconsistent` — the last
  for an alternate whose own annotation does not reproduce its published
  proteins, which no longer applies to any of the 189 rows. A fetched alternate is verified
  against the pinned length and both checksums before anything is written.
  `--inventory-molecules` keeps each row on the accession its inventory names
  first, and `--no-annotation-fill` keeps strictly to each entry's own record;
  either way the choice, and any source that could not be had, is recorded and
  reported by `orfdata viral summary`.
- Fill the annotation gap in three GDC rows whose GenBank submissions carry no
  CDS feature at all — HPV2, HPV8 and HPV51, the last an IARC Group 1 agent in
  the `cancer-causing` subset — from the Papillomavirus Episteme reference the
  GDC inventory already names for every papillomavirus row. PaVE publishes the
  identical sequence for all three, so the 26 transcripts and proteins it adds
  (including spliced E1^E4 and E8^E2 models) need no coordinate transfer. The
  fill applies only where the source record annotates nothing, only on an exact
  sequence match, and only when a feature's coordinates reproduce the protein
  PaVE publishes for it; `annotation_source` records it on every transcript and
  genome, and `--no-annotation-fill` (`fill_annotation_gaps=False`) keeps
  strictly to the source records. The fill adds to those records rather than
  standing in for them, so a PaVE that is unreachable or no longer publishing
  the same sequence leaves one genome without transcripts, with its reason
  reported, instead of discarding a download of 200 source records.
- Read GenBank flat files with `orfdata.genbank` instead of Biopython, adding
  no dependency in its place: the subset NCBI EFetch emits is small, and the
  parsers available either lose what the file says — biotite keeps a feature's
  locations in a `frozenset`, so the order of a spliced CDS is gone — or are
  too little used to put under a published package. The reader keeps a
  location's structure and its `<`/`>` partial markers, and is checked against
  an independent parser (`gb-io`, a development dependency) on all 200 records
  and 1,678 CDS features with no difference. Where the two disagree, ORFdata
  keeps the more faithful reading: `gb-io` reads `order(...)` as a join and
  cannot parse the `complement(order(>104408,...))` on an EBV gene at all. Translation uses ORFdata's own genetic codes, so nothing
  silently falls back to the standard code. Viral support now needs nothing
  beyond the base install; only the data is downloaded.
- Join a CDS that a circular genome's origin splits into two partial halves.
  HBV's polymerase is written as `2307..>3215` and `<1..1623` in `X04615.1`,
  and reads as one 2,532-base gene of 843 amino acids. The halves are joined
  only when the result proves it — whole codons, a start codon, a terminal
  stop, no stop in between, and the two published translations meeting end to
  end — and the joined transcript records `genbank:origin_join` as its
  `annotation_source` and keeps both halves' own descriptions, which name the
  ranges `(AA 1-303)` and `(AA 304-843)` of the one protein. A CDS whose record
  describes its product only in a note, as HBV's does, now keeps that note as
  its product rather than nothing.
- `Proteome.write_annotation(directory)` and `orfdata proteome --annotation
  DIRECTORY` write the same three files for noncanonical ORFs, and
  `orfdata --format gtf` writes the ORF models alone. The ncORF GTF writes no
  `gene` rows, so it can be concatenated with the assembly's own annotation.
  See the new [Downstream tools](https://openvax.github.io/orfdata/downstream/)
  page.
- `orfdata download` provides one discovery, status and acquisition hierarchy
  for viral collections, catalogues, datasets, native search references, pinned
  genomes and chain files; the existing resource-specific commands remain
  available. An unknown kind, species or context is an error rather than an
  empty listing.
- `select_viral_references()` creates reusable union or intersection selections
  for nonexclusive contexts and named subsets without giving an
  accession/version a context-specific identity. It rebuilds a damaged index
  from the local bundle instead of failing, and excluding a context without
  naming one now starts from every context.
- The package YAML registry exposes three named selectors:
  `cancer-integrated-driver` (15 references), `cancer-causing` (20 IARC Group 1
  references) and `cancer-associated` (37 IARC Group 1, 2A and 2B references).
  Every abbreviation is checked against the GDC inventory when the registry
  loads, so a misspelled one is an error instead of a silently smaller subset,
  and the listed count is what the subset resolves to. Being outside one tier
  does not classify a reference as non-carcinogenic.
- `orfdata setup --species human` installs the default viral collection;
  `--collections none` (or `collections="none"`) sets up catalogues and chains
  alone, and a collection whose optional dependency is missing is reported as
  skipped rather than failing the setup. DataCache, PyYAML and gb-io are base
  dependencies, so a regular install can download viral references.
- `human` is an analysis target; source-host taxon IDs remain empty when the
  supplied source inventory does not establish them. The published context is
  the 200-reference GDC cancer-transcriptome screen (also selectable as
  `cancer`), and contexts are nonexclusive metadata that never determine cache
  identity.

- Resolve PyEnsembl 2.31 reference strains using `species_taxon_id`, while
  preserving registered organism identities (including domestic dog), explicit
  strain checks and compatibility with earlier PyEnsembl releases.

- Replace the proteome constructor's separate ncORF and catalogue arguments with
  `nc_orf_sources`: a string or list of dataset names/evidence presets, a composed
  evidence rule, or loaded `ORFCatalog` objects. Empty selections add no ncORFs;
  lists take the union. `"all"` and evidence presets search every compatible
  registered dataset, including optional ones. Missing resources fail with
  installation instructions instead of silently omitting datasets.
- `ReferenceTranscriptome.proteome`, `Proteome.build`, and `with_nc_orfs` share
  this selection. `nc_orfs`, `catalog`, `catalogs` and legacy constructor
  selectors are removed. Use `nc_orf_sources=catalog` for loaded data. The
  lower-level `with_extra_orfs(evidence, catalog=...)` API remains available.
- The proteome CLI uses `--nc-orf-sources` for presets, dataset names and local
  files; repeat it or use comma-separated strings. Separate catalogue and
  evidence flags are removed from this command.

# 0.35.0

Proteomes can now select reference CDS proteins, transcript-bound ncORFs and
host-associated viral proteins in one build. `ReferenceTranscriptome` separates
the transcript models from the protein sources. `proteome()` defaults to
reference CDS proteins, with `isoforms="all"` or `"longest"` to select every
isoform or one protein per gene. ncORFs are opt-in: `nc_orfs="ms"`, `"ribosome"`,
`"all"` or an evidence rule selects them; `nc_orfs=True` means `"ms"`.
`proteome(reference_cds=False, nc_orfs="all")` builds an ncORF-only proteome
through the same API. FASTA headers and summaries expose `protein_source`,
while the existing `kind` labels remain compatible with saved snapshots.
`Proteome.build` and the CLI use the same reference-only default; the CLI adds
`--isoforms longest` and `--nc-orfs NAME`. The older selectors remain supported
on `Proteome.build`, and `with_extra_orfs` is an alias for `with_nc_orfs`.
Supplied sample-derived
host sequences can be added with `DerivedProtein`, preserving parent accessions,
sample and change IDs without treating them as reference translations. Snapshot
schema v2 also reads earlier v1 snapshots.

DLA ligand release 2026-10-16.2 marks iRT spike-in spectra that the first
release missed.

- **iRT reads.** `kaabinejadian-dog-dla-ligands` now defaults to release
  2026-10-16.2 (importer `kaabinejadian-dla-v2`), which marks 10,629 rows as
  `irt_standard`, up from 2,242, and names the rule in `irt_rule`:
  - `irt_sequence`: a kit peptide or fragment, with I and L equal;
  - `kit_only`: PEAKS attributes the row only to the kit;
  - `irt_precursor`: a de novo, homology or kit-less database-search read at
    an iRT precursor;
  - `de_novo_irt`: a de novo sequence matched to the kit, or sharing 6
    residues with a kit peptide.
- **Who's affected.** Half of all de novo rows were iRT reads. The
  database-search ligand lists barely change: 12,948 dog ligands across all
  databases (was 12,949); the UniProt and cell-line counts are unchanged.
  Marked rows' mappings are `unmapped`.
- **First release.** Release 2026-10-16.1 stays available and rebuilds
  byte-identically.
- **Missing runs.** The dataset docs now explain which runs behind 1,695
  rows' best spectra are missing from the PRIDE deposit, and why.

# 0.34.1

- Require Python 3.11 or later and test through Python 3.14.

# 0.34.0

- Viral reference queries now return immutable `ViralReference` objects with
  `.dna` and `.orfs` properties. Use `viral_reference(id)` to select an exact
  sequence by stable ID; browsing by host, group and type still returns every
  matching reference. The properties load on access, so discovery does not
  construct ORFs or copy DNA into every result.
- This replaces the dictionary rows returned by `viral_references()` in
  0.33.0. The existing `viral_reference_sequence(id)` and
  `viral_orfs(reference_dna, ...)` functions remain available; alignment
  annotation still checks actual reference DNA. `ref.to_dict()` preserves
  the old JSON metadata shape for persisted coverage reports.

# 0.33.0

- Discover the bundled viral coding references by documented host, convenience
  group and virus type, while keeping viral taxonomy and host taxa separate.
  HPV16's human association comes from PaVE; the type's ICTV species and
  source URLs remain in each reference record. An unlisted host is not an
  exhaustive host-range conclusion.
- Include both checksum-verified HPV16 reference sequences offline. Call
  `viral_reference_sequence(id)` or `viral_orfs(reference_id=id)` without
  obtaining genome DNA first. The existing `viral_orfs(reference_dna, ...)`
  continues to verify actual alignment-reference sequence identity.

# 0.32.0

Dog peptide data: MHC class I ligands and the Canine PeptideAtlas.

- **`kaabinejadian-dog-dla-ligands`.** Every peptide row of PRIDE PXD074485's
  PEAKS tables (Kaabinejadian et al., iScience 2026; CC0): DLA class I ligands
  from four dog tumours (osteosarcoma, two anal sac adenocarcinomas, oral
  melanoma), each searched against four dog protein databases, and from three
  human HCT 116 lines expressing one DLA-88 allele each. Samples keep breed and
  DLA typing; rows keep PEAKS area, intensity, spectrum count, modified form and
  identification type (database search, de novo or homology). iRT standards
  are flagged in every search, and de novo sequences aren't mapped to the
  proteins PEAKS lists. Published as a data release;
  `orfdata datasets download kaabinejadian-dog-dla-ligands`.
- **`peptideatlas-dog-peptides`.** The Canine PeptideAtlas 2024-01 build
  (Kusebauch et al., Sci Data 2026): PSM counts for 149,546 peptides in 18
  samples, each peptide's protein locations with ROS_Cfam_1.0 coordinates, and
  17,795 protein calls, 17,110 of them mapped to Ensembl 110 proteins by
  identical sequence.
  PeptideAtlas states no licence, so this release is built locally from the
  original download (`orfdata datasets build peptideatlas-dog-peptides`) and
  never redistributed.
- **Build-only datasets.** A dataset release with no published artifact
  (`"assets": {}`) is built from its pinned sources. `DatasetRegistry.download`
  explains how to build it, `open` and `inspect` use the local build, and
  `DatasetRegistry.published(name)` tells the two apart.
- **Storage version 5.** Target mappings no longer repeat their reference for
  every candidate; the DLA dataset is a third smaller. Older files remain
  readable.
- Neither dataset searched for noncanonical ORFs, so dog still has no ncORF
  catalogue; the catalogue guide lists what was checked.

# 0.31.0

- Experimental pinned RefSeq proteome inputs (#102).

# 0.30.0

Ribo-seq ORF catalogues for zebrafish, C. elegans and yeast, and two optional
human ones.

- **Zebrafish, C. elegans and yeast.** Yang et al.'s RibORF 2.0 ORFs
  (Nat Commun 2024) for each species: `riborf-zebrafish` (9,759 ORFs on
  GRCz11, Ensembl 115), `riborf-worm` (2,849 on WBcel235, Ensembl 105) and
  `riborf-yeast` (1,128 on R64-1-1). `orfdata setup --species zebrafish` (or
  `worm`, `yeast`) builds them from their sources' files. `zebrafish`, `worm`
  and `c_elegans` are known species names.
  - Yeast ORFs are placed from their genomic spans, because the source's yeast
    transcripts carry UTRs that Ensembl's lack. A record names an Ensembl
    transcript only when it carries the ORF; only those 171 place on a
    PyEnsembl genome, and the rest stay queryable in the catalogue.
- **Human, optional.** `ncorf-mammals-human` (Chang et al., 11,621 ORFs) and
  `riborf-human` (Yang et al., 51,174) are used only when asked for:
  `load_catalog(species="human", catalogs="all")`, `extra_orfs(genome, evidence,
  catalogs="all")`, `Proteome.build(..., catalogs="all")` or `--catalogs all`.
  Human's default stays the bundled GENCODE catalogue, so loading it stays
  offline and fast, and its selections are unchanged. Install them with
  `orfdata setup --species human --catalogs all`.
- **Choosing catalogues.** `catalogs=` takes `None` (the species' default set),
  `"all"` or catalogue names (`"bundled"` for the bundled one), in
  `load_catalog`, `extra_orfs`, `Proteome.build`, `Proteome.with_extra_orfs`
  and `install_species`. A registry entry marked `"optional": true` is used
  only when named or with `"all"`; `catalog_status()` reports `optional`.
- **Dog.** No dog catalogue: no public dog Ribo-seq exists, and every dog MS
  dataset searched canonical proteins only. The catalogue guide lists what
  was checked.
- The RibORF importer keeps C. elegans transcript names with dots whole, and
  lists ORFs over ambiguous bases as unresolved. The mouse catalogue it builds
  is byte-identical to 0.28.0's.
- "No ORF catalogue covers" errors list the catalogues by species.

# 0.29.0

`viral_references()` lists the bundled, sequence-qualified viral coding
references. `viral_orfs(reference_dna, contig=...)` returns existing
`SourceRNA`/`TranscriptORF` models offline for the nine RefSeq NC_001526.4
HPV16 CDSs, including the two joined coding paths. Native and exact GDC
circular-rotation placements retain separate sequence identities and source
pins. These are CDS fragments, not complete transcript models or direct
translation evidence. Other viruses and strains remain explicit coverage gaps.

Biopython is needed only in the build-data/dev extras to rebuild and verify
the pinned source annotation, not to load the released catalogue.

# 0.28.0

Mouse gets two Ribo-seq ORF catalogues, merged with OpenProt's, and human gets
T2T-CHM13.

**Upgrading:** mouse ORFs now come from three catalogues, merged, and all three
must be installed. Run `orfdata setup --species mouse --accept-terms` (with
`orfdata[data,build-data]`; about 130 MB of sources). With only
`openprot-mouse` installed, as 0.27.1 left it, mouse proteomes raise and name
what is missing. The default `"ms"` evidence still selects OpenProt's 2,976
MS-supported ORFs.

## Evidence names

- The MS evidence kinds have species-neutral names: `"proteomics"` (peptides
  from digested cellular proteins; was `"shotgun"`) and `"immunopeptidomics"`
  (peptides eluted from MHC molecules; was `"mhc"`, and `"hla"` for human).
  - `MS(mhc=True)` counts only immunopeptidomics peptides and `MS(mhc=False)`
    only proteomics peptides; `MS()` counts every MS peptide, including those
    whose source doesn't say which. `methods=` still names kinds directly.
  - The old names still work as evidence names and in `MS(methods=...)`, with
    a `DeprecationWarning`. Catalogues, evidence summaries and saved rules that
    use them load with the new names, so the bundled human catalogue is
    unchanged. Catalogues saved by 0.28.0 use the new names, which 0.27.1
    can't read.
  - Proteome FASTA headers say `evidence=ribosome,immunopeptidomics`, not
    `ribosome,mhc`. Saved proteome snapshots keep the labels they were written
    with.
- `evidence="all"` selects every catalogue ORF, like `None`.
- `"riboseq"` is another name for `"ribosome"`.
- A list, tuple or set of names selects ORFs with any of them:
  `["proteomics", "immunopeptidomics"]` is
  `MS(mhc=False) | MS(mhc=True)`.

## Mouse

- Ribo-seq ORFs of Chang et al. (eLife 2025): `ncorf-mammals-mouse`, 16,484
  ncORFs called by at least two of PRICE, RiboCode and Ribo-TISH in at least
  two libraries, on Ensembl 107 GRCm39, built from the authors' pinned files,
  Table S3 and Ensembl 107's GTF.
  - Coding paths come from each ORF's transcript coordinates, because the
    source BED garbles multi-exon blocks.
  - All 16,484 bind exactly on Ensembl 107, and 16,470 place on Ensembl 115.
  - Near-cognate starts read as methionine and are declared on each record.
  - Library and caller counts are aggregate Ribo-seq summaries.
- RibORF ORFs of Yang et al. (Nat Commun 2024): `riborf-mouse`, 40,089 ORFs of
  at least 5 residues that RibORF 2.0 found translated in 224 mouse Ribo-seq
  samples.
  - The source gives mm10 transcript coordinates and no proteins, so each ORF
    is placed on the same Ensembl 107 GRCm39 transcript and translated there,
    without a restricted mm10 chain.
  - 12,222 shorter ORFs (a bare start-stop, for one) are left out. 1,817 rows
    on transcripts that changed, have no sequence or are gone are listed as
    unresolved.
  - 39,983 place on Ensembl 115.
  - IDs start `riborf_`. The source's ouORF, odORF and noncoding types take
    Chang et al.'s names: uoORF, doORF and lncORF.
- `load_catalog` merges every registered catalogue for a species and assembly,
  rather than refusing to choose between two.
  - The merge follows the registry's integer `priority`, so a new catalogue
    never renames existing ORFs.
  - Mouse merges into 48,805 ORFs. `RiboSeq()` adds 46,094 to the Ensembl 115
    proteome, and `RiboSeq() & MS()` adds 371.
  - A catalogue from a source already used on the same assembly raises, naming
    both catalogues. An ORF that two merged catalogues failed to lift is
    reported once.
- In merged catalogues, a biotype query matches any source's label, as gene
  and transcript queries do. TSV export has a `start_codon` column.
- `orfdata setup --species S` (`install_species`) installs everything a
  species needs in one step: every catalogue release any reachable assembly
  uses, and the chain files.
  - It repairs damaged installs.
  - It downloads a source several catalogues share once.
  - A download that stalls for `--timeout` seconds is retried, then fails.
  - It reports each item. If any failed, it exits 1, and `install_species`
    raises `InstallError`, after installing the rest.

## T2T-CHM13

- `orfdata genomes install t2t-chm13` installs Ensembl's rapid-release CHM13
  annotation (2022_07), which needs `orfdata[pyensembl]`. It is pinned by URL,
  size and SHA-256, and `load_genome("t2t-chm13")` opens it offline as a
  PyEnsembl genome.
- The T2T Consortium's public-domain GRCh38↔CHM13 chains are pinned, so the
  bundled catalogue lifts to CHM13 automatically: 7,208 ORFs lift and 6,767
  place.
- `Proteome.build(load_genome("t2t-chm13"))` gives 104,957 annotated proteins
  plus 2,039 MS-supported ORFs, and `orfdata proteome --genome t2t-chm13` does
  the same.
- `pinned_genomes`, `install_genome`, `load_genome` and `orfdata genomes
  list|install` handle annotations that `EnsemblRelease` doesn't cover.
  - The index belongs to the pinned files, and `--force` rebuilds a damaged
    one.
  - The genome keeps its species and assembly accession when pickled. Its
    paths are absolute, so it works from any directory.

## Other changes

- Varcode no longer fails on extra ORFs for a lifted genome.
  `ExtraORFs.excluded` holds only records of its catalogue, and ORFs a lift
  left out move to `ExtraORFs.unlifted`. Proteomes still list them, counted as
  `unlifted`.
- An ORF whose coding path is on a transcript that encodes another protein is
  excluded with that reason and transcript, not "Transcript absent". When its
  start codon reads differently there, the reason says so.
- `merge_catalogs` reads each catalogue's metadata once rather than per
  record, so merging the mouse catalogues takes seconds, not minutes.
- `load_catalog` reuses a lifted or merged result within a process, handing
  each caller its own copy.
- Install hints name the `--cache-root` they were given.

# 0.27.1

- `orfdata.reference.species_name_for_taxon(taxon_id)` resolves a known NCBI
  species taxonomy ID to its canonical name offline, without an annotation
  (#83, #84).

# 0.27.0

- Catalogues on other assemblies are lifted automatically. `load_catalog`,
  `extra_orfs(genome)` and `Proteome.build(genome)` now work on every human
  and mouse assembly PyEnsembl knows, once the chains are downloaded:
  - Human GRCh37 places 6,857 ORFs and NCBI36 4,195.
  - Mouse GRCm38 places 2,856 and NCBIM37 1,620, the latter in two hops
    through GRCm38.
  - Each catalogue contributes one release; one per assembly is used, and
    catalogues from several assemblies are merged.
  - Missing catalogues or chains raise an error with the command that
    installs them. Nothing is downloaded while loading.
- `orfdata chains download` with no assemblies downloads every pinned chain
  (`download_chains()` returns what it downloaded and skipped). UCSC's
  restricted ones are skipped without `--accept-terms`; `chains list` shows
  what is downloaded.
- Explicit catalogues on another assembly are lifted too, in `extra_orfs` and
  `Proteome.build`. Missing catalogues and chains are reported together,
  under the registry's or `--cache-root`'s root. Catalogues on assemblies no
  chain reaches are ignored, a source already used from a nearer assembly
  counts once, and `hg38`/`mm39`/`GRCh38.p14` and other species' assemblies
  get errors that say what to use. `catalog_status` follows exactly what
  `load_catalog` uses.
- Ensembl's GRCh38↔NCBI36 and GRCm38↔NCBIM37 chains are pinned.
- `chain_routes`, `lift_to` and multi-hop `orfdata lift`.
- `catalog_status()` reports the assemblies each catalogue reaches, and
  whether their chains are downloaded.

# 0.26.0

- Mouse noncanonical ORFs from OpenProt 2.2 (#18, #55):
  `orfdata catalogs build openprot-mouse` builds 2,976 GRCm39 records from
  OpenProt's pinned TSV, BED and FASTA downloads. That's every protein with at
  least two unique MS peptides on an Ensembl 106 transcript. orfdata does not
  redistribute OpenProt's files.
  - `load_catalog(species="mouse")`, `extra_orfs` and `Proteome.build` then use
    it. Mouse Ensembl 115 gains 2,854 ORFs; 2,864 bind exactly on Ensembl 106.
  - RefSeq rows are listed as unresolved.
  - MS and translation-event scores are aggregate summaries.
  - Every row is checked against the verified coordinate offset.
- A registry release with no published artifact is built locally from its
  pinned sources. `load_catalog`, `CatalogRegistry.open`/`inspect` and
  `orfdata catalogs query` open such a build. A missing or damaged
  installation gives the build or download command, `--force` included.
  Without `orfdata[data]` the error says to install it, instead of raising
  ImportError.
- A source accession split across several coding paths (`IP_1@ENSMUST…`)
  answers to its bare ID, which names all of its records; an ID that names
  several records always raises rather than picking one.
- A placement on the source transcript is never repeated, and a lifted ORF is
  on its source transcript under the version its new assembly gives it.
- `catalog_status()` and `orfdata catalogs status` list every catalogue,
  bundled included, with species, assembly, record count, installation and
  which is selected (#70).
- `extra_orfs` places a record pinned to another release of the assembly by
  path and protein on its own transcript too, so exclusions give the real
  reason (for example "sequence unavailable") instead of the release
  mismatch. On another version of that transcript, `source_binding` says
  "ORF transcript version mismatch".
- Tests run with a private `ORFDATA_CACHE`.

# 0.25.0

- `orfdata lift --to BUILD` lifts the bundled catalogue, a species'
  catalogue or a local `--catalog FILE`, writing the result and a summary of
  what did not lift (#72). `orfdata chains list` and `orfdata chains download
  FROM TO` manage four pinned chain files: Ensembl's GRCh37/GRCh38 pair and
  UCSC's GRCm38/GRCm39 pair, which need `--accept-terms` for UCSC's
  non-commercial licence. In Python, `pinned_chains()`, `download_chain()` and
  `load_chain()`. A chain is checked against its pinned size and SHA-256 and
  never downloaded implicitly; `--source-path` installs a local copy the same
  way.
- `ChainFile` checks a pinned file's direction and species, so the wrong
  file can't silently relabel a catalogue; `orfdata lift` also refuses a lift
  to the catalogue's own assembly.
- `ORFCatalog.save(path)` writes JSON or gzipped JSON deterministically. TSV
  exports of merged or lifted catalogues include `members_json`,
  `liftover_json` and `lifted_from_json` instead of failing.
- The `orfdata` query command takes `--min-ribo-seq-studies`,
  `--min-ms-peptides`, `--ms-methods` and `--either`, like `orfdata proteome`
  (#73). `orfdata catalogs query` takes them too. `--genome-build` alone picks
  the catalogue for that assembly.

# 0.24.0

- `Proteome.build(genome, evidence, catalog=, viral=, cache=True)` builds the
  usual proteome in one call and keeps it in the OpenVax DataCache root, keyed
  by reference, evidence rule, catalogue and viral identities, snapshot schema
  and orfdata version; an identical later build loads the verified snapshot
  (#74). The key also covers the pyensembl version and, for a custom genome,
  its local annotation files' sizes and SHA-256. `orfdata proteome --cache
  [--cache-root DIR]` does the same. `cache_root` is also where registered
  catalogues are looked for, and builds killed partway leave no files behind
  for more than a day. A species without a catalogue says to pass
  `extra_orfs=False` for annotated proteins only.
- Snapshots save deterministically, with no gzip timestamp, and faster: about
  8 s instead of 15 s for human Ensembl 115 with MS-supported ORFs. A cache
  hit loads in about 7 s; a first cached build takes about 30 s.

# 0.23.0

- Lift catalogues between assemblies of one species with `ChainFile` and
  `lift_catalog` (#72). Records lift only when every coding block lies in one
  ungapped alignment block on one contig and strand, preferring the
  same-named chromosome over fix patches; the rest are listed with the
  reason. Lifted records keep their original in `lifted_from`, the chain's
  identity in `liftover`, and any strain pin. 7,247 of the 7,264 bundled ORFs
  lift to GRCh37 and 6,857 place on Ensembl 75; all that lift back return to
  their GRCh38 coordinates. `contig_names=` renames RefSeq-named targets, and
  UCSC scaffold names read as Ensembl's. ORFs that don't lift are reported
  among `extra_orfs` and proteome exclusions under the rule that selects them.
  A lifted record can lead a merge. When only another assembly has a
  catalogue, `load_catalog`'s error says which to lift.
- `extra_orfs` places records pinned to another annotation release of the
  same species, strain and assembly, by their exact coding path and protein.
  It also tries the transcript the source names when strict binding fails,
  for example on another version of it; strict `record.bind` still requires
  the pinned reference.

# 0.22.0

- `merge_catalogs` and `catalog=[...]` combine ORF catalogues for one species
  and assembly (#71, #72). Records from different catalogues with the same
  spliced coding path, protein, genetic code and declared start are one ORF:
  it keeps the first catalogue's ID, every source's record in `members`, and
  the union of their evidence. The others' IDs appear as `source:id` in
  `equivalent_ids` on placements and proteome entries and in FASTA headers.
  `ORFCatalog.get` and `ExtraORFs.for_orf` accept any source's ID; an ID two
  sources use for different ORFs is qualified, and the bare ID raises. Merges
  require one strain, record ORFs that share a path but differ in protein or
  code as conflicts, and can be repeated. Merged catalogues round-trip through
  JSON. `orfdata proteome` and `orfdata` queries take repeated `--catalog`
  files, and queries find records by any source's ID. A study two
  catalogues label differently counts once when their PMIDs agree. Queries
  and placements match any source's transcript and gene. A redundant
  declared ATG start doesn't split an ORF, and records join an ORF on their
  own transcript first.
- `ORFCatalog.source` gives the release label that proteome layers record.

# 0.21.0

- Set Ribo-seq and MS stringency separately with evidence rules (#73):
  `RiboSeq(min_studies=2)` and `MS(min_peptides=2, unique=True, methods=...)`,
  combined with `&` (both) or `|` (either). `with_extra_orfs`, `extra_orfs`,
  `ORFCatalog.query` and `has_evidence` accept them. Method names such as
  `"ms"` and `"hla"` are the simplest rules and select exactly as before.
  Peptides count once, including I/L variants, and nested peptides count as
  one: the count is the largest set of supporting peptides none of which lies
  inside another, so stricter rules always select subsets. A source's
  aggregate count stands for one observation. Rules have one canonical form
  (`RiboSeq() & MS() == MS() & RiboSeq()`), and `EvidenceRule.from_dict` reads
  stored ones back. Aggregate-only MS support is labelled `ms` in proteome
  entries. Proteome
  layers store the rule as data, replacing the separate `filters` field, so
  snapshot identities change from 0.20.0; `ExtraORFs.evidence` holds the rule.
- `orfdata proteome` takes `--min-ribo-seq-studies`, `--min-ms-peptides`,
  `--ms-methods` and `--either`.

# 0.20.0

- Describe ORFs by RNA coordinates with `TranscriptORF.at`: declared non-AUG
  starts, fragments without a known start, recoded stop codons such as
  readthrough and selenocysteine, and ORFs across circular-RNA junctions on a
  `SourceRNA` (#64, #65, #66). Catalogue models can declare `start_codon`.
- Annotated CDSs report `start_codon`, `stop_codon`, decoded internal stop
  codons and sub-codon introns, such as programmed frameshifts. Ensembl omits
  start-codon rows for GUG, ACG and AUU starts, so a CDS that begins after its
  RNA's first base now counts as started whatever its codon. One at the first
  base without a row is still reported incomplete, because PyEnsembl does not
  keep the `cds_start_NF` tag that would tell the two apart (20 human Ensembl
  115 transcripts, such as GNAS `ENST00000306120`).
- `cds_relation` relates any ORF to its transcript's annotated CDS (upstream,
  internal, N-terminal extension and so on), alongside the source's category.
  For a CDS that begins with a partial codon, the frame is read back from its
  stop codon; protein residues count from Ensembl's leading `X`.
- `ProteinProduct` and `products_from_uniprot` map precursor products such as
  POMC's hormones to codons and genomic positions (#67).
- `RNAFeature` and `rfam_ires` locate reported IRES elements from the bundled
  Rfam 15.1 human cellular seeds (CC0) on transcripts (#63).
- An examples page covers bicistronic RNAs, ORF categories, non-AUG starts,
  recoding, frameshifts, ERV proteins, circular RNAs, precursors and IRESs.
- `load_catalog(species=..., genome_build=...)` picks a species' ORF catalogue:
  the bundled one, or else one downloaded through `CatalogRegistry`. It defaults
  to human, as PyEnsembl does, never downloads, and opens each registered
  catalogue's default release unless the assembly needs another. `extra_orfs`
  and `with_extra_orfs` use it, so a mouse or dog catalogue added to the packaged
  registry and downloaded needs no `catalog=`. A species without a catalogue
  gets an error naming those that exist. The ORF query CLI takes `--species`.

# 0.19.0

- Build layered reference proteomes from a PyEnsembl genome: annotated proteins,
  evidence-selected noncanonical ORFs and host-restricted viral proteins. Every
  entry lists each transcript encoding it; layers record role, source,
  selection and catalogue identity; source-binding history and exclusions are
  kept. Checksummed snapshots reload offline (#68). Adds `orfdata proteome`.
- `TranscriptORF` (formerly `ReferenceORF`, still available) represents annotated
  CDSs and catalogue ORFs on a transcript, with Python-slice `start`/`end`
  coordinates and source transcript/gene identity. `extra_orfs(genome)` places
  each ORF on every transcript carrying its exact coding path, with
  reference-checked transcript and gene lookups and a genomic `overlapping()`
  query; `annotated_orfs()` finds annotated CDSs carried by a transcript, such as
  a bicistronic RNA's upstream protein (#61, #62). A CDS without a marked
  start or stop is placed on another transcript only where that RNA continues
  identically beyond it. `orf.initiation_context` reads start flanks from the
  transcript each placement is on. Transcript IDs that contain dots, such as
  C. elegans `Y74C9A.3.1`, are compared whole.
- Bind ORFs under declared NCBI tables 2 and 3 (vertebrate and yeast
  mitochondrial codes), and normalize contig names as PyEnsembl does.
- Include aggregate evidence in generic TSV exports (#58).
- Publish documentation with MkDocs on GitHub Pages.

# 0.18.0

- Manage optional ORF catalogues through explicit, checksum-verified discovery,
  download, source-build and offline query APIs backed by OpenVax DataCache.
- Retain the existing ORFCatalog/PyEnsembl model with declared source/reference
  contexts and independent accepted/unresolved coverage.

# 0.17.0

- Preserve source-scoped aggregate ORF evidence counts alongside resolved peptide
  and ribosome observations, with explicit resolution filters and shared identity.
- Keep counts distinct from abundance, peptide uniqueness, spectrum review and
  evidence for modified variant/RNA products.

# 0.16.0

- Import native yeast rich/starvation RNA and ribosome-footprint region counts,
  with source-reported read densities, normalization inputs and sample metadata.
- Preserve historical source-window identity and unknown genomic coordinates;
  do not infer complete ORFs or unqualified cross-assay aliquot matches.

# 0.15.0

- Select paired observations by scoped left/right targets and native per-side
  quantity, status, unit and scale filters in Python and the CLI.
- Filter relationships and assay observations in SQLite before reconstructing
  unrelated records, retaining pairing strength, missing mates and attribution.

# 0.14.0

- Support source-defined region targets for RNA/Ribo measurements over UTRs,
  introns and quantification windows, through the existing scoped query API.
- Preserve the distinction between a measured region and a complete translated
  ORF; retain explicit mapping, reference and comparison requirements.

# 0.13.0

- Import the six native Wang healthy-human-tissue abundance series with historical
  Ensembl 83 references, adjacent-specimen pairing and source calibration.
- Preserve censored copy estimates, native numeric text and protein-group
  ambiguity; keep gene-level comparisons distinct from individual ORF abundance.
- Read complete pinned Ensembl peptide FASTAs, including gzip sources, exact
  versioned identifiers, alternate-reference entries and native stop characters.
  Header taxonomy remains unknown when absent; explicit source species is retained.

# 0.12.0

- Reopen and inspect local source builds offline through `open_build` and
  `inspect_build`, preserving source provenance and validating build identity.
- Add explicit `--built` selection to CLI queries, pairing views and inspection;
  published artifacts remain the default, with no download or fallback on reads.

# 0.11.0

- Import native Pfizer/Oziolor dog and mouse RNA gene counts, with measured zeros,
  sample/run metadata and exact Ensembl 110 gene/reference validation.
- Keep mouse biological strains separate from the reference genome, retain
  differing source tissue labels, and leave unresolved RNA/protein pairing explicit.
- Deduplicate three equivalent GEO count deposits; do not relabel nonhuman
  ortholog matrices as human measurements. Missing human counts remain tracked in #37.

# 0.10.0

- Expose complete pinned Pfizer human/mouse/dog search FASTAs through explicit
  reference downloads and offline Python/CLI queries, independently of abundance.
- Preserve all native sequences and headers, source terms and scoped identities,
  including search entries without reported protein-group representatives.
- Distinguish search membership from detection and unresolved group assignment;
  keep quantitative artifact hashes and historical cache identities unchanged.

# 0.9.0

- Select exact native `unit`, `scale` and `log_base` in measurement queries,
  mapped views, CLI queries and Parquet exports, including storage v1-v4.
- Preserve unfiltered behavior, values, zeros, missingness and source-row lineage.
  An unavailable series returns no observations; no conversion or fallback occurs.

# 0.8.0

- Compact new quantitative datasets with integer reference/target/assay keys
  and lossless reconstruction of repeated measurement fields (storage v4).
  Python/CLI records and scoped identities remain unchanged.
- Keep published storage v1/v2/v3 readable; validate against databases created
  by the corresponding released versions, including censoring and input lineage.
- Populate index statistics so narrow target queries use the target index.
- On a real 100,000-row RNA slice, reduce installed size from 121.6 to 38.6 MB
  and preserve every observation field. Existing data artifacts are unchanged.

# 0.7.0

- Import CoPY yeast SRM protein/peptide quantities, distinct label-free estimates,
  aggregate RNA and four RNA libraries through the same species-scoped dataset API.
- Preserve upper bounds, missingness, quality status, standards, uncertainty and
  calibration. RNA copies/cell retain their 60,000-mRNA assumption and input lineage.
- Add explicit cohort RNA/protein comparisons; individual replicate pairing,
  grouped-protein attribution and unspecified historical references stay unresolved.
- Register a checksum-verified 69,125-observation data release and source rebuild,
  using stable versioned PMC objects instead of a regenerated ZIP (#33).
- Retain four standard-sequence mismatches and duplicate label-free rows with
  visible provenance; further source reconciliation is tracked in #34.

# 0.6.0

- Import the complete public Pfizer/Oziolor human, mouse and dog protein-group
  summary tables through one species-scoped API, retaining native iBAQ/rIBAQ,
  measured zeros, source missingness and unresolved group representatives.
- Preserve source sample identities, TMT channels, pools, normalization and
  pinned search FASTAs. No inferred RNA pairing or protein copies-per-cell.
- Add explicit source rebuilds with `DatasetRegistry.build` and
  `orfdata datasets build`, using DataCache checksums and local source overrides.
- Register versioned compressed data artifacts, including the complete canine
  summary in an 84.4 MB download, with source counts and redistribution terms.

# 0.5.0

- Represent sequence-only and incomplete source references with pinned inputs
  and explicit null metadata, without inventing a genomic assembly (#28).
  Existing complete-reference keys remain unchanged.
- Add source-to-annotation mapping records and views that retain original
  measurements, all candidate targets, validation limits and unmapped reasons.
  Shared groups stay groups; gene RNA cannot map directly to ORF abundance.
- Support mapped Python/CLI queries with explicit ambiguity and sequence
  validation filters. SQLite storage version 3 adds mapping tables while older
  datasets remain readable without inferred mappings.
- Pair RNA and protein assays with different native references within a species
  (#30). Derived values retain explicitly scoped input observations, including
  colliding local IDs. Gene-to-protein-group links remain ambiguous context.

# 0.4.0

- Add explicit assay pairings and within-reference target relationships,
  including missing mates, temporal metadata and source-backed attribution.
- Query paired RNA, ribosome and protein measurements with caller-selected
  pairing bases, preserving native values and technical/quantity multiplicity.
  Parent RNA and shared groups cannot become ORF-specific ratios.
- Add named native-unit ratios with input lineage and explicit unresolved zero
  or censored denominators, plus opt-in pinned orthology comparisons that retain
  both species' original observations.
- Extend SQLite storage to version 2 for relationships while reading original
  version 1 datasets without inferring pairings. Add `orfdata datasets pairs`.

# 0.3.0

- Add immutable quantitative manifests, samples, assays, measurements and target
  mappings, with native units, missingness, provenance and explicit derivations.
- Store streaming datasets in reference-partitioned SQLite tables with validated
  sample/assay/target relationships; species and dataset releases never share IDs.
- Add optional DataCache downloads and source-pinned import materialization,
  offline discovery/inspection/query APIs and `orfdata datasets` commands.
- Add exact target construction from bound reference ORFs and PyEnsembl features,
  and an optional Parquet export. The initial quantitative registry is empty;
  source-specific data releases follow separately in the roadmap.

# 0.2.0

- Add explicit taxonomy, strain, assembly and annotation contexts and scoped
  targets, using the same API for arbitrary species. Preserve source names and
  resolve known dog aliases without conflating assemblies or breeds.
- Pin optional ORF source contexts and retain resolved contexts and sequence
  fingerprints in bound-reference exports. Custom genomes missing species or
  taxonomy now require an explicit `ReferenceContext` when binding.
- Validate real human, mouse, dog and budding yeast reference slices offline.
  These fixtures do not add nonhuman experimental protein evidence.
- Include reference fixtures, docs and check scripts in source distributions,
  and validate the extracted source package in CI (#26).

# 0.1.0

- Move the versioned Isovar reference catalogue into a disease-independent package.
- Separate ORF models, source evidence, PyEnsembl bindings and variant modifications.
- Add a short introduction to the PyEnsembl evidence layer, including noncoding
  RNA and alternative ORF coverage, evidence interpretation and current limits.
