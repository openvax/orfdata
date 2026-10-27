# Yeast protein and RNA abundance

`copy-yeast-abundance` imports [Lawless et al., CoPY (2016)](https://doi.org/10.1074/mcp.M115.054288)
through the same dataset API as human, mouse and dog. It includes calibrated
protein and peptide quantities, label-free estimates and RNA from budding yeast
grown in glucose-limited chemostats. The strain is EUROSCARF Y11335 BY4742,
with the reported genotype retained on each sample.

```sh
pip install 'orfdata'
orfdata datasets info copy-yeast-abundance
orfdata datasets download copy-yeast-abundance
orfdata datasets pairs copy-yeast-abundance --basis condition_aggregate
```

The release contains 69,125 observations in a roughly 14 MB compressed download
(168 MB installed). Downloads are explicit and checksum-verified through
DataCache. Queries work offline after installation.

| Source measurements | Coverage | Interpretation |
|---|---:|---|
| QconCAT SRM protein summaries | 1,903 rows | 1,167 quantified; 668 upper bounds; 68 unquantifiable |
| SRM peptide observations | 3,836 rows × four cultures and author median | Protein-specific isotope standards; copies/cell, classification and uncertainty |
| Aggregate RNA | 6,581 labels × FPKM and estimated copies/cell | Counts aggregated before FPKM; copies assume 60,000 total mRNAs/cell |
| Individual RNA libraries | Four × 7,137 tracking IDs | Native FPKM, 95% confidence interval and Cufflinks status |
| Label-free estimates | 3,166 rows × two methods | Q-Exactive internal-standard Top3 and SAX total-protein scaling remain separate |

The 55 joined SRM protein groups retain their combined attribution. Label-free
identities also remain groups because full membership is unresolved. Two rows
share the label YNL239W; both survive unchanged. Upper bounds have
`status="below_detection"`, `value=None` and a `detection_limit`; missing and
unquantifiable values never become zeros. Author-calibrated SRM estimates have
`origin="observed"`; the RNA copies/cell conversion has `origin="derived"`
with its original FPKM input and calibration assumption.

The QconCAT FASTA contains 93 synthetic calibration standards, not native yeast
proteins. Four reported peptides do not match their assigned standard; their
measurements retain an explicit validation flag. See [source reconciliation #34](https://github.com/openvax/orfdata/issues/34).
Raw transitions and light/heavy ratios are not supplied in the processed
workbook. SRM raw data are identified by PASS00717; PXD002694 is the label-free
deposit. RNA is from GSE73898. Cufflinks LOWDATA rows retain their original
positive/zero estimates and confidence bounds in provenance, while the query
value remains `not_measured` because reliable quantification is unresolved.

## Comparing RNA and protein

The paper supports same-cell measurements at the source-cohort level. The
individual MS-to-RNA replicate crosswalk is unresolved, so comparisons require
explicit `condition_aggregate` selection. Replicate numbers alone never join
samples. Label-free methods have different calibrations and no inferred RNA pair.

```python
from orfdata import DatasetRegistry

dataset = DatasetRegistry().open("copy-yeast-abundance")
match = next(
    m for m in dataset.matched_measurements(allowed_bases=["condition_aggregate"])
    if m.relation.left.identifier == "YML035C"  # AMD1, yeast source context
)
protein = match.right[0]  # 3203.943997 copies/cell; author SRM estimate
rna = next(m for m in match.left if m.unit == "copies/cell")  # ~2.11665
ratio = match.ratio(
    protein, rna, measurement_id="AMD1:protein-per-RNA",
    method="native calibrated protein / estimated RNA",
    normalization={"assumed_total_mrna_per_cell": 60000, "level": "source cohort"},
)
print(ratio.value, ratio.unit, ratio.derivation)
```

The ratio is conditional on the RNA calibration. It is neither a translation
rate nor ribosome occupancy. Group relationships cannot yield an individual
protein ratio. All quantities keep their source units and assumptions.

## References and other species

Both yeast contexts carry taxon 4932 and pinned source hashes. The historical
SGD/Bowtie genome and annotation releases are unspecified; assembly, annotation
version and reference strain remain null. These source labels support canonical
gene/protein relationships but cannot bind a current PyEnsembl transcript until
a reference mapping is validated. Sample BY4742 strain is separate from the
unknown reference-genome strain. This dataset does not establish which yeast
noncanonical ORFs were searched, expressed or absent.

For cross-species comparisons, supply independently validated, pinned orthology
records rather than joining gene names. For example, with an installed mammal
RNA dataset and a curated JSON list of `OrthologyMapping` records:

```python
import json
from orfdata import OrthologyMapping, compare_orthologs

with open("curated-yeast-mammal-orthology.json") as handle:
    mappings = [OrthologyMapping.from_dict(row) for row in json.load(handle)]
for comparison in compare_orthologs(dataset, mammal_rna_dataset, mappings):
    print(comparison.left, comparison.right)  # Separate original measurements
```

The mapping must name both exact references, its release, source checksum,
cardinality and confidence. No expression or protein-existence evidence transfers
between species, and different cultures/calibrations remain different experiments.
See [pairing and orthology](pairing.md).

To rebuild the release from the pinned public workbook, standards, metadata and
four RNA files, install `orfdata` and run
`orfdata datasets build copy-yeast-abundance`. Python uses
`DatasetRegistry().build("copy-yeast-abundance")`. Local verified inputs can be
supplied with repeated `--source-path ASSET=PATH` or `source_paths={...}`.
The workbook and standards use stable versioned PMC objects; the dynamically
generated EuropePMC ZIP is not used as a content identity.
