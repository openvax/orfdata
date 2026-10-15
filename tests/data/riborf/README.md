# Yang et al. RibORF mouse fixtures

Real source slices for `tests/test_riborf.py`, for testing the importer; not a
catalogue release.

- `SD3.xlsx`: eight rows of the `mouse` sheet of Supplementary Data 3 of Yang
  et al., Nat Commun 15:1932 (2024), doi:10.1038/s41467-024-46240-9 (CC BY
  4.0): a uORF, an N-terminal truncation and an extension, ACG and GTG starts,
  a uORF on a transcript that changed length since GENCODE vM25 (unresolved),
  and an annotated CDS and a readthrough ORF (not imported).
- `ensembl107.gtf.gz`, `ensembl107.cdna.fa.gz` and `ensembl107.ncrna.fa.gz`:
  Ensembl 107 GRCm39 GTF lines and sequences for their eight transcripts,
  from <https://ftp.ensembl.org/pub/release-107/>, under Ensembl's
  [data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).

The full sources' SHA-256 are pinned in `orfdata/data/catalogs.json`.
