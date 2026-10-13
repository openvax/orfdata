# Matched human tissue RNA and protein abundance

The [Wang et al. atlas](https://doi.org/10.15252/msb.20188503) measures adjacent
cryosections from 29 healthy human tissue specimens. `wang-human-tissue-abundance`
keeps six original table series in the same species-scoped API as dog, mouse and
yeast data.

| Source table | Targets | Native quantity |
| --- | ---: | --- |
| EV1 protein groups | 15,210 | iBAQ |
| EV1 genes | 13,640 | Author-normalized iBAQ |
| EV2 transcripts | 93,674 | FPKM |
| EV2 genes | 18,072 | Author-normalized FPKM |
| EV5 protein groups | 15,210 | Estimated copies/cell |
| EV5 transcripts | 60,519 | Estimated copies/cell |

There are 6,273,425 source cells. The same specimen has several processing views;
these are not additional biological replicates. Native specimen IDs, RNA/MS run
IDs, units, source cells and processing methods remain attached to the data.
Bone marrow has no matched RNA column and is outside the 29-tissue dataset.
The download is approximately 397 MB and installs to a 3.86 GB SQLite database.
Every source cell and RNA-copy input link was independently reconciled before
publication. Missing or censored cells remain explicit observations.

```sh
orfdata datasets download wang-human-tissue-abundance --version 2026-10-12.1
```

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()
dataset = registry.open("wang-human-tissue-abundance", "2026-10-12.1")
reference, = dataset.manifest.references
gene = ScopedTarget(reference, "gene", "ENSG00000004059", "ensembl")  # ARF5
for value in dataset.measurements(target=gene, sample_id="V122"):
    print(value.quantity, value.value, value.unit, value.status)
```

Select the same gene on both sides for a specimen's native gene-level comparison:

```python
match, = dataset.matched_measurements(
    allowed_bases=["adjacent_donor_tissue"],
    pair_ids=["rna_genes:protein_genes:V122"],
    left_target=gene, right_target=gene,
    left_filters={"unit": "FPKM"}, right_filters={"unit": "iBAQ"},
)
print(match.left[0].value, match.right[0].value)  # 51.09, 112559620.9
```

With Ensembl 83 already installed in PyEnsembl, the same target can be built with
`ScopedTarget.from_pyensembl(gene_object, kind="gene", namespace="ensembl",
context=reference)`. The context pins the historical GRCh38.p5 assembly accession,
GTF and complete Ensembl 83 protein FASTA. Alternate-reference protein entries
absent from the GTF remain in the native search set. Search lookups preserve
versioned FASTA IDs; they do not silently strip versions or imply detection.

Protein-group values remain group observations. Representative proteins and their
parent genes/transcripts are context, with unresolved individual attribution.
The source's gene-level summaries permit gene-level RNA/protein comparisons;
they do not quantify a particular isoform, uORF or ncORF. Pairings explicitly use
`adjacent_donor_tissue`, independently of the target relationship. See the
[comparison API](measurements.md) for selecting pairings and deriving ratios.

Protein copies are the authors' **proteomic-ruler estimates**, not isotope-standard
measurements. Raw intensity inputs and per-sample calibration constants are absent
from these tables; corresponding iBAQ observations are linked as associated data,
not falsely recorded as the ruler's inputs. RNA copies link their source FPKM
cells and retain the assumption that mRNA constitutes **2% of total RNA mass**.
Unreported calibration uncertainty stays unknown.

Protein zeros are unquantified values; gene-RNA and transcript-copy zeros retain
the source censoring semantics. They are not measurements of biological absence.
The methods state a 1-FPKM threshold, but gene values undergo median normalization;
2,677 positive gene values are below 1. Their values remain intact, and no numeric
detection limit is invented in normalized units.
Unfiltered transcript FPKM retains its native zeros and tiny numeric text values;
a nonzero estimate does not establish expression above the study's 1-FPKM cutoff.
Original numeric text and zero sentinels remain in provenance. Stable identifiers
join the source tables even when their display gene labels differ.

The [Eraslan reanalysis](https://doi.org/10.15252/msb.20188513) selects major
isoforms and normalizes exon/intron counts differently. It is a separate planned
view of the same specimens, not silently mixed into this release. These tables
cover canonical coding-gene expression; they add RNA/protein context without
expanding the experimental noncanonical-ORF catalogue.

Source supplements retain CC-BY-4.0 attribution, and Ensembl annotations retain
their source terms. Discovery and installed queries work offline; acquisition and
source rebuilding are explicit through DataCache.
