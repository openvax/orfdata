# Released data and gaps

ORFdata keeps datasets in separate species/reference contexts and exposes them
through the same discovery, download and query API. This table describes actual
data coverage, independently of what the API can represent.

| Species | RNA abundance | Protein abundance | RNA/protein pairing |
| --- | --- | --- | --- |
| Human | Wang: 29 specimens, gene/transcript FPKM and estimated copies/cell | Wang: 15,210 groups, iBAQ and estimated copies/cell; Pfizer: 10,352 groups, 37 channels, iBAQ/rIBAQ | Wang: adjacent cryosections; Pfizer: native RNA missing |
| Mouse | Pfizer: 606 libraries, 56,051 gene targets, counts | Pfizer: 12,124 groups, 63 channels, iBAQ/rIBAQ | Unresolved library-to-pool correspondence |
| Dog | Pfizer: 434 libraries, 30,941 gene targets, counts | Pfizer: 9,929 groups, 38 channels, iBAQ/rIBAQ | Unresolved library-to-pool correspondence |
| Budding yeast | CoPY: four RNA libraries plus a source aggregate, FPKM and estimated copies/cell | CoPY: calibrated SRM and separately calibrated label-free estimates | Explicit source-cohort pairing; individual cross-assay replicate matches unknown |

Gene targets include measured zeros; their presence in a matrix does not imply
expression in every sample. Protein-group representatives do not establish
unique protein attribution. TMT channels are not independent donors. CoPY RNA
copies/cell assume 60,000 total mRNAs; its 1,903 SRM summaries contain 1,167 point
estimates, 668 upper bounds and 68 unquantifiable results. Calibration and
uncertainty remain attached to observations. Human Wang protein copies use the
proteomic-ruler estimate; RNA copies assume mRNA is 2% of total RNA mass. These
calibrations remain distinct from yeast isotope-standard measurements. Wang
processing views do not add specimens, and gene-level comparisons do not resolve
individual protein isoforms or noncanonical ORFs.

The separate existence-evidence catalogue contains 7,264 human noncanonical
ORFs, including 2,186 with reported MS support. It has no released mouse, dog or
yeast experimental ncORF catalogue yet. Its ribosome evidence and pooled PSM
counts are not sample-level ribosome occupancy or protein abundance. No species
currently has a released quantitative Ribo-seq dataset in this registry.

See the short guides for [matched human tissues](human-tissue.md),
[protein measurements](pfizer-protein.md),
[native dog/mouse RNA](pfizer-rna.md), [yeast RNA/protein calibration](copy-yeast.md)
and [complete search sequences](search-references.md). The
[roadmap](https://github.com/openvax/orfdata/issues/3) tracks further human
comparators, quantitative Ribo-seq, ncORF evidence and additional species.
Cat, cattle, rat and other taxa use the same model when source datasets are
added. Orthology comparisons are explicit and never transfer observations or
protein-existence evidence between species.
