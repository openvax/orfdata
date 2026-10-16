# Kaabinejadian et al. DLA ligand fixtures

Real source slices for `tests/test_dog_peptides.py`, for testing the importer;
not a dataset release.

- `Results.zip` holds six of PRIDE PXD074485's PEAKS peptide tables, with
  their original headers and member names, and 33 of their data rows (numbered
  from the first after the header):
  `Lola-BB76-UniProt9615` (rows 1, 3, 6, 13, 15, 25, 42, 64, 69),
  `Lola-BB76-MergedCanFams` (1, 5, 10, 11, 56, 60, 112, 145),
  `Lola-H58A-UniProt9615` (4, 9),
  `HCT116-DLA8850101-BB76-SwissProt9606` (3, 14, 20, 55, 57, 98),
  `Tumor163828A-BB76-UniProt9615` (20, 643, 1203, 1587, 1988, 2043) and
  `Tumor163828A-BB76-MergedCanFams` (670, 1757). They cover
  database-search, DeepNovo and Homolog rows, modified forms, Biognosys iRT
  standards (including a kit fragment that the merged database, which lacks
  the kit, attributes to dog proteins), a de novo row listing an iRT entry among
  dog proteins, a row without an accession, rows without a PEAKS area, a
  peptide in both of Lola's BB7.6 searches, and Lola's two immunoprecipitations.
  The rows added for importer `kaabinejadian-dla-v2` are misread iRT spectra
  with their kit peptides (a de novo read and a homology read at a kit
  peptide's precursor, and a kit-less database-search read), and
  database-search rows at kit precursors that stay unmarked.
- `manifest.json` describes the slice for importer `kaabinejadian-dla-v2`;
  `manifest-v1.json` describes the same slice for `kaabinejadian-dla-v1`.
- `sdrf.tsv` is the deposit's complete community-annotated SDRF.

Both are CC0 (PRIDE PXD074485; Kaabinejadian et al., iScience 2026,
doi:10.1016/j.isci.2026.115975). The full archive's SHA-256 is pinned in
`orfdata/data/datasets.json` and in `manifest.json`'s `upstream` metadata.
