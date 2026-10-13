These are small, real-data slices of Wang et al. (2019),
https://doi.org/10.15252/msb.20188503 (CC-BY-4.0). Original versioned source URLs,
SHA-256 values and byte sizes are in `expected.json` and each fixture asset's
`original_source` metadata. Zip/workbook bytes are newly encoded fixture slices;
selected source cells, identifier strings and numeric text are unchanged.

The fixture retains adrenal gland, colon and duodenum columns; quantitative
sample metadata also includes bone marrow, which lacks matched RNA. Native rows
cover ARF5, KRT75, a multi-protein KRT13 group, INO80B/INO80B-WBP1's differing
copy-table gene label, an alternate-reference group, and very small numeric
FPKM text. Original source row numbers and column labels are in `expected.json`.
The protein-gene sheet retains its 24 trailing empty rows. The fixture has 237
abundance cells, not the full study's coverage.

Ensembl release 83 GTF and peptide FASTA slices retain complete native records
for selected genes and representative proteins, including alternate-reference
entries. Attribute Ensembl/EMBL-EBI and preserve its source terms:
https://www.ensembl.org/info/about/legal/disclaimer.html
Historical reference counts in the fixture manifest describe only the slice.

Native copy values are source estimates: proteomic-ruler protein calibration
and RNA calibration assuming mRNA is 2% of total RNA mass. Zero copy values
are censored/unquantified source values, not biological absence. Group values
are not individual protein or noncanonical ORF abundance.
