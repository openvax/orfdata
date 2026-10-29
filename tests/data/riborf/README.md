# Yang et al. RibORF mouse fixtures

Nineteen exact Supplementary Data 3 rows from Yang et al., Nat Commun
15:1932 (2024), DOI 10.1038/s41467-024-46240-9 (CC BY 4.0). The canonical
and readthrough rows are excluded. All remaining 17 reconstruct against native
GENCODE M25; 13 have exact full-RNA placements on Ensembl 107 and four retain
native models with unresolved target RNA. The unchanged short CDS in the
changed chr14 RNA is not promoted to a uniquely validated placement.

`gencodeM25.gtf.gz` and `gencodeM25.rna.fa.gz` retain exact native annotation
lines and FASTA records for these source transcript versions. The existing
`ensembl107.*` slices include their stable-ID target counterparts and every
full-RNA-identical target candidate found in the complete target FASTAs. The
source workbook also contains the species fixtures documented below.

`v3-manifest.json` pins these factual subsets and their coverage under their
actual reference/version labels. Full source pins are in the packaged registry.
GENCODE source files are from its official release_M25 directory:
https://www.gencodegenes.org/mouse/release_M25.html
Native GTF SHA-256: e6c9704cb5832b8a1c71431ecac9ea30057377c72cdc9a12dc5d0405ab066d5a
Native RNA SHA-256: 5139cf409a680b60434716ca10ad4e474cbbe5d2e4d11af7553edfad8ed73670
GENCODE/EMBL-EBI data-use terms and Ensembl data-use terms apply to their
reference data. Citation and original release identities are retained.

## Zebrafish, C. elegans, human and yeast

`SD3.xlsx` also has `zebrafish`, `worm`, `human` and `yeast` sheets with real
rows of the same Supplementary Data 3 sheets, used by `tests/test_more_species.py`:

- **zebrafish, worm, human:** for each, an annotated CDS (not imported), a
  short ORF, and a uORF, a near-cognate start, an lncRNA ORF, an iORF and an
  N-terminal truncation, each on the shortest such transcript; for zebrafish an
  ORF over an ambiguous base, and for human one on a transcript without
  sequence, both unresolved. C. elegans has no unresolved row.
- **yeast:** seven rows on chromosome XIV: a truncation inside its Ensembl
  transcript, a uORF in a UTR Ensembl's transcript lacks, an ORF on an unstable
  transcript (CUT808), an extension across CPT1's intron, Ty1's frameshifted
  ORF (YNL284C-B, unresolved), a short ORF and an annotated CDS.

`zebrafish-ensembl115.*`, `worm-ensembl105.*` and `human-ensembl105.*` hold the
GTF lines (with their genes' lines), cDNA and ncRNA sequences of those
transcripts from Ensembl 115 GRCz11, 105 WBcel235 and 105 GRCh38.
`yeast-ensembl115.gtf.gz` and `yeast-ensembl115.cdna.fa.gz` hold Ensembl 115
R64-1-1 records of the yeast rows' transcripts, and `yeast-ensembl115.dna.fa.gz`
chromosome XIV. All from <https://ftp.ensembl.org/pub/>, under Ensembl's
[data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).
