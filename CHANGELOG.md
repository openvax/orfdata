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
