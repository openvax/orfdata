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

## Zebrafish, C. elegans, human and yeast

`SD3.xlsx` also has `zebrafish`, `worm`, `human` and `yeast` sheets with real
rows of the same Supplementary Data 3 sheets, used by `tests/test_more_species.py`:

- **zebrafish, worm, human:** for each, an annotated CDS (not imported), a
  short ORF, and a uORF, a near-cognate start, an lncRNA ORF, an iORF and an
  N-terminal truncation, each on the shortest such transcript; for zebrafish an
  ORF over an ambiguous base, and for human one on a transcript without
  sequence, both unresolved. C. elegans has no unresolved row.
- **yeast:** six rows on chromosome XIV: a truncation inside its Ensembl
  transcript, a uORF in a UTR Ensembl's transcript lacks, an ORF on an unstable
  transcript (CUT808), a spliced ORF (unresolved), a short ORF and an annotated
  CDS.

`zebrafish-ensembl115.*`, `worm-ensembl105.*` and `human-ensembl105.*` hold the
GTF lines (with their genes' lines), cDNA and ncRNA sequences of those
transcripts from Ensembl 115 GRCz11, 105 WBcel235 and 105 GRCh38.
`yeast-ensembl115.gtf.gz` and `yeast-ensembl115.cdna.fa.gz` hold Ensembl 115
R64-1-1 records of the yeast rows' transcripts, and `yeast-ensembl115.dna.fa.gz`
chromosome XIV. All from <https://ftp.ensembl.org/pub/>, under Ensembl's
[data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).
