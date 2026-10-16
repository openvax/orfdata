# Kaabinejadian et al. DLA ligand fixtures

Real source slices for `tests/test_dog_peptides.py`, for testing the importer;
not a dataset release.

- `Results.zip` holds four of PRIDE PXD074485's PEAKS peptide tables, with
  their original headers and member names, and 22 of their data rows (numbered
  from the first after the header):
  `Lola-BB76-UniProt9615` (rows 1, 3, 6, 13, 15, 42, 64, 69),
  `Lola-BB76-MergedCanFams` (5, 10, 11, 56, 60, 145),
  `Lola-H58A-UniProt9615` (4, 9) and
  `HCT116-DLA8850101-BB76-SwissProt9606` (3, 14, 20, 55, 57, 98). They cover
  database-search, DeepNovo and Homolog rows, modified forms, Biognosys iRT
  standards (including a kit fragment that the merged database, which lacks
  the kit, attributes to dog proteins), a de novo row listing an iRT entry among
  dog proteins, a row without an accession, rows without a PEAKS area, a
  peptide in both of Lola's BB7.6 searches, and Lola's two immunoprecipitations.
- `sdrf.tsv` is the deposit's complete community-annotated SDRF.

Both are CC0 (PRIDE PXD074485; Kaabinejadian et al., iScience 2026,
doi:10.1016/j.isci.2026.115975). The full archive's SHA-256 is pinned in
`orfdata/data/datasets.json` and in `manifest.json`'s `upstream` metadata.
