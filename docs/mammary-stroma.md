# Paired canine mammary-stroma RNA and protein

`canine-mammary-stroma` imports the processed measurements of
[Amini et al. (2021)](https://doi.org/10.1016/j.neo.2021.03.001),
[GSE135183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135183) and
[PXD023023](https://www.ebi.ac.uk/pride/archive/projects/PXD023023).
It uses the same data commands and `QuantitativeDataset` API as other species
and sources:

```sh
orfdata data info canine-mammary-stroma
orfdata data fetch canine-mammary-stroma
orfdata data path canine-mammary-stroma
```

```python
from orfdata import DatasetRegistry

registry = DatasetRegistry()
dataset = registry.open("canine-mammary-stroma")  # Offline; no implicit download.
rna_reference, protein_reference = dataset.manifest.references
pairs = list(dataset.pairings())
complete_pairs = [pair for pair in pairs if pair.right_assay_id is not None]
pair = complete_pairs[0]
rna = dataset.measurements(
    reference=pair.reference, assay_id=pair.left_assay_id, unit="estimated_count"
)
protein = dataset.measurements(
    reference=pair.right_reference, assay_id=pair.right_assay_id, unit="LFQ intensity"
)
```

The release contains 993,560 native measurements in a 52 MB download and a
514 MB installed SQLite file. Rebuilding uses approximately 8 MB of pinned processed inputs; no raw MS
files, sequencing reads, alignment index or vendor database is required.

| Layer | Native scope | Quantities |
| --- | --- | --- |
| RNA | 18,582 Ensembl gene IDs × 30 samples, 15 cases | Integer estimated counts from Kallisto/tximport |
| Protein | 1,583 native MaxQuant groups × 28 samples, 14 cases | Intensity, LFQ intensity, iBAQ |
| Peptide | 5,413 native peptide rows × 28 samples | Intensity and the separately reported LFQ intensity columns |

Different quantitative columns do not add biological replicates. The 58
assays retain their original sample and run identities. Native MS rows include
77 reversed protein groups and 52 contaminant groups; these counts are not
counts of confidently detected canine proteins.

## Specimen pairing

The authors reused archived FFPE mammary specimens for separate laser-capture
RNA and protein preparations. The importer verifies a unique match between
clinical Table 1 and GEO by breed, age and histological subtype, recording
spelling normalizations and both original sets of labels. This is a
**metadata-derived crosswalk within an author-declared shared cohort**, not an
author-provided case-ID mapping.

Fourteen shared cases give 28 complete assay pairs: normal stroma and
cancer-associated stroma separately. Two RNA samples from GEO case 127 have
explicitly missing protein mates; the paper reports exclusion of proteomics
case 7 because its paraffin block was unavailable. No protein sample is invented.

These pairs use `same_specimen`. This basis means separate preparations from the
same original specimen, without claiming identical aliquots or adjacent
sections. Temporal labels refer to the original biological specimen, not the
dates of RNA sequencing and mass-spectrometry acquisition. Missing mates retain
an unknown temporal relationship. Pairing does not itself establish a
relationship between RNA genes and protein groups: this release has no
cross-modal target relations or current-reference mappings.

The normal compartment belongs to a cancer-bearing mammary specimen. It is not
a healthy-animal normal-tissue atlas. Source case IDs identify individual
specimens, and repeated compartments or quantities never imply additional dogs.

## Missingness, native records and reference limits

- RNA zeros remain measured zero counts. Counts are estimated gene-level values,
  with no conversion to TPM, copies per cell or translation efficiency.
- MaxQuant MS zeros denote unquantified cells. Their numeric values are `None`
  with status `missing`; the original `"0"` survives in provenance. The 2,156
  native `"NaN"` iBAQ cells also remain missing, with a separate reason.
- Every native MS field survives once in the corresponding mapping's
  `provenance["native_record"]`. Measurements retain source row/column, original
  numeric text, identification type and the original control/decoy flags.
- Peptide mappings retain every reported native protein-group ID. A unique
  assignment to one *group* is not unique protein or isoform detection. Group
  records retain reported protein labels, including contaminants, without
  assigning them to a canine biological protein reference.
- No source QC threshold, aggregation, imputation, differential-expression model
  or RNA/protein ratio is applied. Consumers can inspect the native flags and
  scores explicitly.

RNA was quantified against Ensembl transcripts on CanFam3.1; its exact release
and transcript sequence bytes are not reported. The protein search used a
canonical UniProt dog reference dated 2019-07-29 with decoys and contaminants;
the searched FASTA is absent from the deposited processed files. Both references
remain source-bound and incomplete. Source membership and uniqueness flags are
preserved, but no full-background sequence uniqueness is asserted. These data
provide neither ncORF translation evidence nor MHC-presentation evidence.

The [earlier source audit](dog-data-sources.md) records the acquisition and
crosswalk evidence. The PRIDE project declares CC0; GEO numerical data and
factual metadata use [GEO's public-data terms](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
Factual clinical curation cites Table 1 of the CC BY-NC-ND article. Its original
text and figures are excluded from the dataset artifact.
