# Dog and mouse RNA abundance

The [Oziolor et al. atlas](https://doi.org/10.1038/s41597-026-07512-9) deposits
native gene counts alongside its [protein measurements](pfizer-protein.md).
ORFdata preserves each species separately, using the same measurement API.

| Dataset | RNA libraries | Measured genes | Observations | Reference |
| --- | ---: | ---: | ---: | --- |
| `pfizer-atlas-dog-rna` | 434 | 30,941 | 13,428,394 | ROS_Cfam_1.0, Ensembl 110 |
| `pfizer-atlas-mouse-rna` | 606 | 56,051 | 33,966,906 | GRCm39, Ensembl 110 |

The dog download is about 534 MB and installs to 3.41 GB; the mouse download is
1.27 GB and installs to 8.66 GB. Queries select observations from the installed
SQLite database without loading the whole matrix into memory.

The original `rcnt` values are integer gene counts after the authors' modified
nf-core/Salmon pipeline. They are neither TPM nor calibrated RNA copies per cell.
Zeros stay measured zeros: 4,218,202 dog observations and 18,019,555 mouse
observations. Every value retains the native numeric text and source line.

```sh
pip install 'orfdata[data]'
orfdata datasets download pfizer-atlas-dog-rna --version 2026-10-12.1
orfdata datasets query pfizer-atlas-dog-rna --sample 1 --unit counts
```

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()
dataset = registry.open("pfizer-atlas-dog-rna", "2026-10-12.1")
reference, = dataset.manifest.references
gene = ScopedTarget(reference, "gene", "ENSCAFG00845000002", "ensembl")
for observation in dataset.measurements(target=gene, unit="counts"):
    print(observation.sample_id, observation.value, observation.status)
```

All source gene IDs were checked against the checksum-pinned Ensembl 110 GTF.
The import retains the source's coding and noncoding genes without a
protein-coding filter.
With that annotation already installed in PyEnsembl, construct the same target
from a gene object:

```python
from pyensembl import EnsemblRelease

genome = EnsemblRelease(110, species="dog")
feature = genome.gene_by_id("ENSCAFG00845000002")
gene = ScopedTarget.from_pyensembl(feature, kind="gene", namespace="ensembl",
                                 context=reference)
```

An incompatible species, assembly or annotation release fails explicitly.
There are 180 dog and 890 mouse annotation genes outside the measured count
universe; querying one yields no observation, not a synthesized zero. Gene RNA
abundance is context for a parent RNA and does not establish the abundance or
translation of any particular uORF, ncORF or transcript isoform.

Samples retain original numeric SIDs, GEO/BioSample accessions, SRA runs, batches,
native tissue labels and the public source metadata. Mouse libraries include
308 source-labelled `c57` and 298 `cd1` samples; these are biological strain
labels, independent of the annotation reference. Donor IDs and sex are unknown
at the deposited sample level. More specific GEO tissue labels remain alongside
the count table's broader labels. There is no verified RNA-library-to-protein-pool
crosswalk, so this release creates no assay pairs or protein-per-RNA ratios.

GSE325522, GSE325523 and GSE325525 each carry the same decompressed count matrix.
ORFdata pins one input and records the equivalent mirrors instead of counting
the data three times. The full matrix contains 2,183 nonhuman libraries. Although
the metadata lists 44 human libraries, neither the native CSV nor the deposited
human-ortholog RDS contains their counts. That [source gap](https://github.com/openvax/orfdata/issues/37)
remains explicit. Rat, pig and macaque rows are present in the original matrix
and remain candidates for separate releases through the same importer.

`registry.build(name, source_paths=...)` or `orfdata datasets build NAME` rebuilds
from pinned counts, GEO metadata, factual Table S5 sample records and the
historical GTF; workbook rebuilding needs `orfdata[build-data]`. Python returns
a `QuantitativeDataset` ready to query. The CLI returns its path; explicit CLI
queries of local builds are tracked in [#43](https://github.com/openvax/orfdata/issues/43).
Normal registry queries open the independently verified published artifact.

GEO data and Ensembl annotations retain their source reuse terms and attribution.
The paper supplement has separate CC-BY-NC-ND-4.0 terms; the original workbook
is not bundled in the quantitative artifact. Source hashes, terms and processing
details are in each manifest. All acquisition is explicit; installed queries
and metadata discovery work offline.
