# Canine PeptideAtlas fixtures

Real source slices for `tests/test_dog_peptides.py`, for testing the importer;
not a dataset release. They are factual extracts; PeptideAtlas states no
licence, and orfdata never redistributes the full export.

- `atlas_build.tsv.zip` holds the tables the importer reads from the build 573
  export (`atlas_build_573.tsv.zip`, 2024-04-30), with their original headers:
  six peptides (PAp00000001, PAp16545961, PAp00494376, PAp07488729,
  PAp16560428 and PAp07498639), all their protein locations, their PSM counts
  in three search batches (brain, liver and timsTOF depleted plasma), seven
  protein calls, and the proteins those name. They cover a peptide shared by
  UniProt, Ensembl and contaminant proteins with exon-spanning coordinates,
  one unique to a UniProt protein, one only in contaminants, one placed with
  isoleucine for leucine, one only in decoys, and calls that map to one, two
  or no Ensembl 110 proteins.
- `ensembl.pep.fa.gz` holds the Ensembl 110 dog proteins
  (`Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz`) whose sequences equal
  those proteins', plus the file's first entry, under Ensembl's
  [data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).

Source: Kusebauch et al., Sci Data 2026, doi:10.1038/s41597-026-06647-z;
<https://peptideatlas.org/builds/canine/>. The full files' SHA-256 are pinned in
`orfdata/data/datasets.json` and in `manifest.json`'s `upstream` metadata.
