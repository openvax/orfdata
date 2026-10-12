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
