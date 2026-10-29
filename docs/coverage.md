# Released data and gaps

ORFdata keeps datasets in separate species/reference contexts and exposes them
through the same discovery, download and query API. This table describes actual
data coverage, independently of what the API can represent.

| Species | RNA abundance | Protein abundance | RNA/protein pairing |
| --- | --- | --- | --- |
| Human | Wang: 29 specimens, gene/transcript FPKM and estimated copies/cell | Wang: 15,210 groups, iBAQ and estimated copies/cell; Pfizer: 10,352 groups, 37 channels, iBAQ/rIBAQ | Wang: adjacent cryosections; Pfizer: native RNA missing |
| Mouse | Pfizer: 606 libraries, 56,051 gene targets, counts | Pfizer: 12,124 groups, 63 channels, iBAQ/rIBAQ | Unresolved library-to-pool correspondence |
| Dog | Pfizer: 434 libraries, 30,941 gene targets, counts; mammary stroma: 18,582 genes × 30 samples | Pfizer: 9,929 groups, 38 channels, iBAQ/rIBAQ; mammary stroma: 1,583 native groups × 28 samples | Mammary stroma: 14 shared specimens, 28 paired compartments and two RNA samples without protein mates; Pfizer crosswalk remains unresolved |
| Budding yeast | CoPY: four RNA libraries plus a source aggregate, FPKM and estimated copies/cell; Ingolia: four RNA libraries, counts and densities | CoPY: calibrated SRM and separately calibrated label-free estimates | CoPY: explicit source-cohort pairing; individual cross-assay replicate matches unknown |

## Peptides and MHC ligands

Dog has four peptide-level datasets:

- **Mammary stroma** ([guide](mammary-stroma.md)): 5,413 native peptide rows, intensity and LFQ columns, reported group memberships and original control/decoy flags across 28 samples.
- **Pfizer native peptide evidence** ([guide](pfizer-peptides.md)): 953,980 target
  peptide sequences, native peptide/PSM reporter quantities and reported protein-group
  memberships. The 96 runs retain 38 tissue pools, four bridge channels and two
  empty channels; source memberships and searched-FASTA candidates stay separate.
- **Canine PeptideAtlas** ([guide](dog-peptides.md), build 573): 149,546 tryptic peptides with PSM counts
  in 18 Labrador tissue, plasma and urine samples, every protein location, and
  17,795 protein calls, 17,110 of them mapped to Ensembl 110 proteins with an
  identical sequence. PeptideAtlas states no licence,
  so this dataset is built locally from its download and never redistributed.
- **DLA class I ligands** ([guide](dog-peptides.md), Kaabinejadian et al.): 76,550 PEAKS rows from four
  dog tumours, each searched against four protein databases, and from three
  human HCT 116 lines that each express one DLA-88 allele.

All four searched known proteins only. Other species have no peptide-level
measurement dataset yet; the human and mouse ORF catalogues carry pooled
peptide evidence for their noncanonical ORFs.

## Reading the table

- Gene targets include measured zeros; their presence in a matrix does not imply
  expression in every sample.
- Protein-group representatives do not establish unique protein attribution.
- TMT channels are not independent donors.
- Wang processing views do not add specimens, and gene-level comparisons do not
  resolve individual protein isoforms or noncanonical ORFs.

Calibration and uncertainty remain attached to observations:

- CoPY RNA copies/cell assume 60,000 total mRNAs. Its 1,903 SRM summaries
  contain 1,167 point estimates, 668 upper bounds and 68 unquantifiable results.
- Human Wang protein copies use the proteomic-ruler estimate; RNA copies assume
  mRNA is 2% of total RNA mass. These calibrations remain distinct from yeast
  isotope-standard measurements.

## ORF evidence and ribosome profiling

The separate existence-evidence catalogue contains 7,264 human noncanonical
ORFs, including 2,186 with reported MS support. Mouse has 60,707 ORFs: Ribo-seq ORFs
from Chang et al. and Yang et al., and 2,976 OpenProt 2.2 ORFs with at least
two unique MS peptides, built from their sources' files. Yang et al.'s Ribo-seq ORFs
cover zebrafish (15,491), C. elegans (4,035) and yeast (1,369), and with Chang et al.'s
are optional human catalogues; dog has no experimental ncORF catalogue, because no
dog Ribo-seq or ncORF MS search is public ([details](catalogs.md#dog)). Its
ribosome evidence and pooled PSM counts are not sample-level ribosome
occupancy or protein abundance.

The separate Ingolia yeast release contains four Ribo-seq libraries across
13,748 source regions, alongside four RNA libraries. Exact region coordinates
and RNA/RPF aliquot pairing remain unresolved. Other species do not yet have
released quantitative Ribo-seq datasets in this registry.

## More detail

See the short guides for [matched human tissues](human-tissue.md),
[protein measurements](pfizer-protein.md),
[native dog/mouse RNA](pfizer-rna.md), [yeast RNA/protein calibration](copy-yeast.md),
[yeast RNA/ribosome occupancy](ingolia-yeast.md)
and [complete search sequences](search-references.md). The
[roadmap](https://github.com/openvax/orfdata/issues/3) tracks further human
comparators, quantitative Ribo-seq, ncORF evidence and additional species.
Cat, cattle, rat and other taxa use the same model when source datasets are
added. Orthology comparisons are explicit and never transfer observations or
protein-existence evidence between species.
