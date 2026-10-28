# Human, mouse and dog protein abundance

The [Oziolor et al. atlas](https://doi.org/10.1038/s41597-026-07512-9) provides
parallel TMT protein measurements. ORFdata keeps each species in a separate
dataset, using the same samples, assays, measurements and query API.

The [dog peptide and protein-group evidence](pfizer-peptides.md) extension adds
native memberships and peptide/PSM quantities while retaining this summary as
same-study context.

| Dataset | Reported protein groups | Source channels | Source rows |
|---|---:|---:|---:|
| `pfizer-atlas-human-protein` | 10,352 | 37 | 325,702 |
| `pfizer-atlas-mouse-protein` | 12,124 | 63 | 615,991 |
| `pfizer-atlas-dog-protein` | 9,929 | 38 | 313,263 |

Release `2026-10-12.1` contains the complete March 23, 2026 summary tables from
MassIVE MSV000099729, MSV000099754 and MSV000099744, respectively. Each source
row contributes two observations: native `iBAQ` and `rIBAQ`. Neither is calibrated
copies per cell. These are pooled tissue samples; rodent condition labels also
include sex. Channels and technical runs do not count as independent donors.

```sh
pip install 'orfdata'
orfdata datasets list --species dog
orfdata datasets download pfizer-atlas-dog-protein
orfdata datasets query pfizer-atlas-dog-protein --quantity protein_abundance --unit iBAQ --scale linear
```

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()
dataset = registry.open("pfizer-atlas-dog-protein")  # After explicit download.
reference, = dataset.manifest.references
group = ScopedTarget(reference, "protein_group", "A0A0A0MPB0", "MSV000099744:ProteinGroup")
for value in dataset.measurements(target=group, unit="iBAQ", scale="linear"):
    print(value.sample_id, value.unit, value.value, value.status)
```

Omit the unit selector to return both native series. They retain the same source
row and sample identity; two units do not imply two independent measurements.

Every target retains its species and the checksum of the deposited UniProt
search FASTA. The original genomic assembly and UniProt release are unknown.
They stay unknown; a future annotation mapping must be explicit and validated.
Group representatives can be inspected through `reference_mappings()`, but their
lookup does not establish all group members. Opting into unresolved mapped
views keeps the original group measurement visible, without assigning its value
to an individual protein. No RNA pairings are inferred in this release.

The [complete search references](search-references.md) are available through
`registry.download_reference(name)` and `registry.open_reference(name)` without
rebuilding abundance. This includes proteins without reported representatives;
search membership is separate from detection or protein-group assignment.

Source `NA` values remain missing and numeric zeros remain measured zeros.
Absent group/channel combinations have no observation. For dog, 64,039 combinations
are absent and eight reported rows contain `NA`, reproducing the existing
64,047 missing iBAQ matrix cells. Each value retains its source line, column,
original numeric text, representative, normalization sum and native `Abundance`
field. Nonfinite transformed strings such as `-Inf` are kept as provenance only.

These searches used UniProt sequences without spliced isoforms. Coverage does
not establish that an unreported ncORF was searched or absent, nor does bulk
tissue MS establish MHC presentation. The article and deposit disagree on the
instrument name; that discrepancy remains explicit in assay metadata.

Downloads are 87.3 MB (human), 154.8 MB (mouse) and 84.4 MB (dog). Installed
SQLite files occupy 1.46 GB, 2.72 GB and 1.40 GB, respectively. Queries stream
observations; they do not load the entire file or a dense matrix into memory.

Source tables are CC0 according to the MassIVE deposits; UniProt search records
retain UniProt Consortium attribution and CC-BY-4.0 terms. Both sources have
pinned hashes and sizes. The data artifact is separate from the Python package.
Discovery and installed queries work offline. `download` verifies and installs
the compressed artifact; it does not fetch the large vendor search database.

To reproduce a source import, use `registry.build(name, source_paths=...)`,
`orfdata datasets build NAME`, or `orfdata.sources.pfizer.build_protein("dog")`.
The source overrides are `protein.tsv` and `search.fasta`; DataCache verifies
their original checksums before running the importer. Source slices with actual
positive, zero and missing values test the same importer for all three species.
