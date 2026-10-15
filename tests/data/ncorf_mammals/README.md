# Chang et al. mouse ncORF fixtures

Real source slices for `tests/test_ncorf_mammals.py`, for testing the importer;
not a catalogue release.

- `mouse.bed` and `mouse.pep.fa`: the rows and proteins of six ORFs from
  `gxelab/ncorf_mammals` at commit `9a7baa318fd53f4faf7fa8fe0030e1e2b0d44ee0`
  (`ncorfs/mouse.bed`, `ncorfs/mouse.pep.fa`; MIT License).
- `S3.xlsx`: the header and the same six mouse rows of Table S3
  (`S3_ncORF_annotation`) from Chang et al., eLife 2025, RP109128 v3
  (doi:10.7554/eLife.109128; CC BY 4.0).
- `ensembl107.gtf.gz` and `ensembl107.cdna.fa.gz`: Ensembl 107 GRCm39 GTF lines
  and cDNA/ncRNA sequences for their five transcripts, from
  <https://ftp.ensembl.org/pub/release-107/>, under Ensembl's
  [data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).

The six cover ATG, CTG, GTG and TTG starts, single- and multi-exon coding paths
(the source BED garbles multi-exon block sizes), and the one ORF in the source
without a stop codon (`13_74704301_74702220_-`), which the importer lists as
unresolved. The full sources' SHA-256 are pinned in `orfdata/data/catalogs.json`.
