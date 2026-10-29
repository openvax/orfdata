# Canine mammary-stroma factual regression slices

Source: Amini et al. (2021), DOI https://doi.org/10.1016/j.neo.2021.03.001,
GEO GSE135183 and PRIDE PXD023023. Acquired and verified on 2026-10-29.

The fixture retains eight native protein-group rows, four native peptide rows,
and three native RNA-count rows, with all original quantitative sample columns.
Source IDs and all selected MS cell values are unchanged; local line numbers
refer to these subsets. Selected rows include real contaminants, reversed decoys,
NaN iBAQ cells, unquantified zero sentinels and shared peptide-to-group mappings.
GEO metadata retains only factual identity and processing fields for all 30
samples. The small `paper.xml` is a synthetic XML wrapper around the factual
clinical Table 1 values, not the article XML or its prose/figures. `project.json`
is a factual accession/license subset. `manifest.json` pins these fixture bytes
under a separate regression dataset ID, with reduced row counts.

Original proteinGroups.txt SHA-256:
67360353510fdb77cc970f671ca9ec17b07686f2dc5771bb6804d182d5c17423
Original peptides.txt SHA-256:
8964af6adbe8ac6fce6897f5eeaaf5da3d7c72fe2eb78ee2653901876a2f0172
Original GSE135183_raw_counts.csv.gz SHA-256:
de3fc2a57fa05682996b328929c94fc51eadd126ec2299870b8f46b24be46a67
Original GSE135183_family.soft.gz SHA-256:
a199bafcc51fefdd3a4231710d0ad6b157a3f9d6c25ae42c26b0bac0c44667c2
Original PMC8042244.xml SHA-256:
c933efb05e263a4d6106c93b2e5344f9ecec08e7875633ec38eef069dfad171a

PRIDE explicitly declares CC0 for PXD023023:
https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD023023
GEO numerical data and factual metadata are publicly distributed under:
https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html
The article is CC BY-NC-ND 4.0; its original text and figures are not included.
Clinical facts are attributed to Table 1. No unique isoform abundance, reference
sequence assignment, ncORF translation or MHC presentation is established.
