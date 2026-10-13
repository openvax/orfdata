# 0.19.0

- Build layered reference proteomes from a PyEnsembl genome: annotated proteins,
  evidence-selected noncanonical ORFs and host-restricted viral proteins, with
  per-layer provenance and reported exclusions. Adds `orfdata proteome`.
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
