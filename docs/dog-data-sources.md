# Canine paired measurements and translation-evidence audit

Source check: **28 October 2026**. This audit identifies importable measurement
files and their limitations. It does not add a dataset to the registry.
The [machine-readable audit](data/dog-source-audit-20261028.json) contains source
URLs, hashes, counts, search queries and a candidate sample crosswalk.

## Paired RNA and protein: a verified study

[Amini et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8042244/)
studied cancer-associated and normal stroma from archived canine mammary
carcinomas. The authors explicitly reused the cases of their earlier RNA study.
Original case 7 was excluded from proteomics because its paraffin block was
unavailable. Measurements come from separate laser-capture preparations of
the same tissue cases, not identical RNA/protein aliquots.

| Source | Verified processed data | Download |
| --- | --- | --- |
| [GSE135183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135183) | 18,582 genes × 30 RNA samples, 15 cases with two stromal compartments | [Raw counts, 606 KB compressed](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135183/suppl/GSE135183_raw_counts.csv.gz) |
| [PXD023023](https://www.ebi.ac.uk/pride/archive/projects/PXD023023) | 1,583 native protein-group rows with 28 LFQ sample columns | [Protein groups, 2.8 MB](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD023023/proteinGroups.txt) |
| PXD023023 | 5,413 peptide rows with native quantitative fields | [Peptides, 4.5 MB](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD023023/peptides.txt) |

The downloaded protein and peptide files match PRIDE's deposited SHA-1 checksums;
the audit also records SHA-256 for all three processed files. RNA values are
nonnegative integer counts. Native protein-group counts include controls and
contaminants and are not counts of confidently detected biological proteins.

PRIDE filenames use case numbers 1–15; GEO uses case IDs such as 109 and 135.
Matching the article's Table 1 to GEO by breed, age and histological subtype
uniquely reconciles all 14 proteomics cases. The resulting 28 sample links are
**metadata-derived**, not an author-supplied identifier crosswalk. Breed spelling
normalizations are recorded. GEO case 127 is the unmatched RNA case.

An importer must retain that pairing basis, separate stromal compartments,
protein-group ambiguity, controls, missing quantities and native units. The
deposited search used a dated canonical UniProt reference, but its FASTA is not
among the 31 PRIDE files. Sequence/reference qualification remains necessary;
these data do not establish ncORF translation. Author-derived contrasts in the
supplementary workbooks must stay distinct from sample-level measurements.
Import work remains tracked in [#20](https://github.com/openvax/orfdata/issues/20).

## Dog Ribo-seq and ncORF translation

The date-stamped dog-specific SRA, GEO and ENA queries in the audit returned no
matching sequencing records. This is a bounded repository search, not proof
that no unpublished or differently indexed experiment exists.

[Pfeffer et al., 2015](https://www.nature.com/articles/ncomms9403)
reports ribosome profiling and total RNA sequencing of canine pancreas rough-ER
preparations aligned to CanFam3.1/Ensembl 77. Its accession section lists EMDB
and PDB deposits, without a sequencing accession. No downloadable reads or
processed ORF-level counts were verified in this audit. Even if recovered, the
preparation would describe an enriched pancreatic compartment, not a broad
normal-tissue canine translatome.

Existing canine PeptideAtlas, Pfizer and DLA datasets retain their deposited
canonical-search scope. The experimental dog ncORF search pipeline is a route
to new analysis, not a released experimentally supported catalogue. Dog ncORF
translation coverage therefore remains unavailable pending qualified data or
a separately validated reanalysis. See [#12](https://github.com/openvax/orfdata/issues/12).
