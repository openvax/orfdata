# Yeast RNA and ribosome occupancy

`ingolia-yeast-rna-ribo` contains the processed RNA and ribosome-footprint
measurements from [Ingolia et al. (2009)](https://doi.org/10.1126/science.1168978),
[GSE13750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13750).
It uses the same dataset API as human, mouse and dog.

```sh
pip install 'orfdata[data]'
orfdata datasets info ingolia-yeast-rna-ribo
orfdata datasets download ingolia-yeast-rna-ribo
```

The download is about 24 MB compressed, 210 MB installed. It contains 329,952
observations: counts, reads/kb and RPKM for 13,748 source regions in eight
libraries. Four libraries measure RNA and four measure ribosome footprints,
in rich medium or after 20 minutes of amino-acid starvation. Source sample
titles identify the modality, condition and nominal replicate. In particular,
GSM346114 is `Footprints-rich-2` despite its `fp_rich1` filename.

```python
from orfdata import DatasetRegistry, ScopedTarget

dataset = DatasetRegistry().open("ingolia-yeast-rna-ribo")
reference, = dataset.manifest.references
target = ScopedTarget(reference, "region", "HRA1", "Ingolia2009:feature")
for measurement in dataset.measurements(target=target, unit="RPKM"):
    print(measurement.sample_id, measurement.quantity, measurement.value)
```

Queries work offline after download. `quantity="rna_abundance"` or
`quantity="ribosome_occupancy"` selects one modality. `unit="counts"` returns
the original aligned read counts. The other two native series have
`origin="derived"` and link to their count observation. Source formulas,
feature lengths and total CDS-aligned reads remain attached; the reported
rounded values are preserved. Zeros remain measured zeros. GSM346118 has no
zero-count rows; ORFdata does not add or remove a pseudocount.

## What these regions establish

The source includes 3,523 5′ UTRs, 4,592 3′ UTRs, 225 intron windows and 5,408
other feature labels, including coding, ncRNA and mitochondrial labels.
Ribosome occupancy over a UTR does not identify a particular uORF. These targets
therefore retain `kind="region"`; they do not expand the experimental ORF
catalogue or establish stable protein abundance, copies/cell or translation flux.

The yeast reference has taxon 4932 and the source annotation date, SGD
2008-06-22. Exact genomic windows and assembly remain unresolved, including
the UTR definitions attributed to Nagalakshmi et al. (2008). The biological
strain BY4741 is sample metadata. Current PyEnsembl coordinates cannot be
substituted for these historical windows without a validated mapping.

The paper reports RNA and footprints from the same samples, but the exact
GSM-to-GSM aliquot crosswalk is unresolved. The release retains eight separate
libraries and creates no assay pairs or translation-efficiency ratios. It also
does not pair these cultures with the separate [CoPY protein study](copy-yeast.md).
[Issue #11](https://github.com/openvax/orfdata/issues/11) tracks these source gaps.

To rebuild, install `orfdata[data,build-data]` and run
`orfdata datasets build ingolia-yeast-rna-ribo`, or
`DatasetRegistry().build("ingolia-yeast-rna-ribo")`. All 17 source inputs have
pinned URLs, sizes and checksums. Verified local inputs can be supplied through
`--source-path ASSET=PATH` or `source_paths={...}`. Query a local build using
`DatasetRegistry().open_build(...)` or `orfdata datasets query ... --built`.
