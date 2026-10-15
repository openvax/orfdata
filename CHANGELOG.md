# 0.29.0

- Mouse Ribo-seq ORFs from Chang et al. (eLife 2025): `ncorf-mammals-mouse`.
  It has 16,484 ncORFs called by at least two of PRICE, RiboCode and
  Ribo-TISH in at least two libraries, on Ensembl 107 GRCm39, built from the
  authors' pinned files, Table S3 and Ensembl 107's GTF.
  - Coding paths come from each ORF's transcript coordinates, because the
    source BED garbles multi-exon blocks.
  - All 16,484 bind exactly on Ensembl 107, and 16,470 place on Ensembl 115.
  - Near-cognate starts read as methionine and are declared on each record.
  - Library and caller counts are aggregate Ribo-seq summaries.
- `load_catalog` merges every registered catalogue for a species and
  assembly, rather than refusing to choose between two, in the registry's
  `priority` order so a new catalogue never renames existing ORFs. Mouse now
  combines OpenProt's and Chang et al.'s ORFs: 19,406 in all, 54 in both,
  named by OpenProt where shared. **Upgrading:** mouse needs both installed;
  run `orfdata setup --species mouse` (with `orfdata[build-data]`).
- A catalogue from a source already used on the same assembly raises rather
  than being skipped; one ORF that two merged catalogues failed to lift is
  reported once. Reference keys are cached, which speeds up opening large
  catalogues.
- `orfdata setup --species S` (`install_species`) builds or downloads every
  registered catalogue for the species and its chain files in one step, and
  repairs a damaged one. Errors about missing pieces suggest it.

# 0.28.0

- T2T-CHM13v2.0 (#72): `orfdata genomes install t2t-chm13` installs Ensembl's
  rapid-release CHM13 annotation (2022_07). It is pinned by URL, size and
  SHA-256, and `load_genome("t2t-chm13")` opens it offline as a PyEnsembl
  genome. The T2T Consortium's public-domain GRCh38↔CHM13 chains are pinned,
  so the bundled catalogue lifts to it automatically: 7,208 ORFs lift and
  6,767 place. `Proteome.build(load_genome("t2t-chm13"))` gives 104,957
  annotated proteins plus 2,039 MS-supported ORFs, and `orfdata proteome
  --genome t2t-chm13` does the same.
- `pinned_genomes`, `install_genome`, `load_genome` and `orfdata genomes
  list|install` handle annotations that `EnsemblRelease` doesn't cover. The
  index belongs to the pinned files, `--force` rebuilds a damaged one, and the
  genome keeps its species and assembly accession when pickled.
- Varcode no longer fails on extra ORFs for a lifted genome:
  `ExtraORFs.excluded` holds only records of its catalogue, and ORFs a lift
  left out move to `ExtraORFs.unlifted` (proteomes still list them, counted
  as `unlifted`).
- An ORF whose coding path is on a transcript that encodes another protein is
  excluded with that reason and transcript, not "Transcript absent".
- `merge_catalogs` reads each catalogue's metadata once rather than per
  record: merging the two mouse catalogues takes seconds, not minutes.
  `load_catalog` reuses a lifted or merged result within a process.

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
