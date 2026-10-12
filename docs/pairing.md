# Comparing RNA, ribosome and protein abundance

Use declared experimental relationships to compare assays. A tissue label alone
does not make an RNA sample and a protein sample a matched pair.

```python
for match in dataset.matched_measurements(allowed_bases=["same_material"]):
    row = match.to_table_row()
    print(row["rna_abundance"], row["protein_abundance"])
```

Rows include named RNA, ribosome occupancy, translation-efficiency, protein,
peptide and spectrum-count columns, each containing original measurement
records with units, missingness and provenance. Missing columns are empty lists.
Multiple values remain multiple values; the API does not average them or create
a Cartesian product of technical replicates. Each row retains its assay pairing
and target relationship. A three-assay experiment can declare several pairwise
comparisons; sharing an RNA measurement does not create independent replicates.

Pairing bases are `same_material`, `adjacent_donor_tissue`,
`biological_replicate`, `independent_cultures`, `tissue_aggregate` and
`condition_aggregate`. The caller must choose the allowed bases. By default,
only explicitly concurrent assays are returned; opt into `longitudinal` or
`unknown` with `allowed_temporal`. Missing mates stay explicit. Technical runs,
TMT channels and bridges remain in the assays and do not become biological
sample counts.

The CLI streams the same table rows:

```sh
orfdata datasets pairs DATASET --version RELEASE --basis same_material
```

## Which targets are being compared?

`TargetRelation` links exact, reference-scoped targets with mapping evidence
and cardinality. `same_target` compares one target across assays/samples.
`canonical_translation` links a gene/transcript to its annotated protein.
`parent_rna_context` links parent RNA to an ORF but keeps the RNA value attached
to its original gene/transcript. `shared_group_member` records an ambiguous
protein-group assignment; the group value stays a group value.

Importers pass `AssayPair` and `TargetRelation` iterators as `pairings=` and
`relations=` to `QuantitativeDataset.create`. IDs belong to the dataset release
and reference context. Invalid assay references, contradictory known donor
identities and conflicting one-to-one assignments fail the import. Original
storage-version-1 datasets remain readable and have no declared relationships.

## Ratios are derived measurements

```python
ratio = match.ratio(
    selected_protein_measurement, selected_rna_measurement,
    measurement_id="ratio-1", method="native protein/RNA ratio",
    normalization={"protein": "source bridge normalization", "rna": "source TPM"},
)
```

Select actual measurements from opposite sides of the match. Ratios require
one-to-one target attribution and linear inputs. Parent-RNA/ORF and shared-group
links cannot produce target-specific ratios. There is no implicit log conversion,
unit conversion, replicate averaging or uncertainty propagation.

An iBAQ/TPM ratio keeps `(iBAQ)/(TPM)` units. Copies-per-cell values also retain
their native composite units and the calibration recorded on their source
assays. Missing, censored or nonpositive denominators produce an unresolved
`not_measured` result with its reason, never infinity or a measured zero. RPF/RNA
ratios are source-defined translation-efficiency measures, not measured flux.
All derived values retain their input IDs, pairing, target attribution and
normalization assumptions.

## Across species

`compare_orthologs(left_dataset, right_dataset, mappings)` is separate and
opt-in. Each `OrthologyMapping` pins its source file, release, both reference
contexts, cardinality, confidence definition and evidence. The default accepts
one-to-one canonical gene mappings; one-to-many relationships require explicit
permission through `allowed_cardinalities`.

The result contains each species' original observations. It neither merges
identities nor transfers abundance or protein-existence evidence. ORF comparisons
additionally require their own sequence, frame and positional evidence and
`allow_orfs=True`; a gene ortholog cannot establish an orthologous ncORF.
Human, mouse, dog, yeast and arbitrary other taxa use this same API.
