# Yang et al. RibORF mouse fixtures

Real source slices for `tests/test_riborf.py`, for testing the importer; not a
catalogue release.

- `SD3.xlsx`: nineteen rows of the `mouse` sheet of Supplementary Data 3 of
  Yang et al., Nat Commun 15:1932 (2024), doi:10.1038/s41467-024-46240-9 (CC
  BY 4.0): a uORF, an N-terminal truncation and an extension, ACG and GTG
  starts, ouORF, odORF and noncoding ORFs (one on a transcript in the ncRNA
  FASTA), a 5-residue ORF, two uORFs on different paths with the same start
  and stop, and ORFs that are not imported: an annotated CDS, a readthrough
  ORF, a 1- and a 4-residue ORF, and three unresolved ones (a transcript that
  changed length since GENCODE vM25, a TEC transcript with no sequence, and a
  transcript whose sequence changed at the same length, on which another ORF
  still places).
- `ensembl107.gtf.gz`, `ensembl107.cdna.fa.gz` and `ensembl107.ncrna.fa.gz`:
  Ensembl 107 GRCm39 GTF lines and sequences for their sixteen transcripts and their genes,
  from <https://ftp.ensembl.org/pub/release-107/>, under Ensembl's
  [data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).

The full sources' SHA-256 are pinned in `orfdata/data/catalogs.json`.
