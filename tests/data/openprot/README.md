# OpenProt mouse fixtures

Real source slices for `tests/test_openprot.py`. They are factual extracts for
testing the importer, not a catalogue release.

`mouse-openprot-2_2-altprots+isoforms-min_2_pep-GRCm39.106.{tsv,bed,fasta}.zip`
hold 16 rows, for eight accessions, from OpenProt 2.2's mouse
`altprots+isoforms-min_2_pep-GRCm39.106` downloads, with the original member
names, header and row order. They are TSV lines 2, 4, 11, 12, 69, 1992,
2981–2983, 3882, 3883 and 5277–5281, and the FASTA entries naming those
accessions. The full archives' SHA-256 are pinned in
`orfdata/data/catalogs.json` (the TSV's, `4b4f4660…`, matches the one recorded
in issue #18). The accessions cover:

- a plus-strand novel isoform (II_085990, Rab1a);
- a minus-strand isoform with 15 coding exons and translation evidence (II_912148);
- one AltProt at two chrY loci (IP_824771);
- one path on transcripts of three overlapping genes, plus RefSeq rows (IP_1064920);
- an accession whose sequence is under another FASTA entry's `PA=` alias (II_1000149);
- an accession with only RefSeq rows (IP_957586);
- an ORF on a TEC transcript that Ensembl's cDNA/ncRNA files leave out (IP_825804).

`ensembl106.gtf.gz` holds the Ensembl 106 GRCm39 GTF lines (gene, transcript,
exon, CDS, start/stop codon) for the ten Ensembl transcripts and their genes.
`ensembl106.cdna.fa.gz` holds their entries from Ensembl 106's cDNA and ncRNA
FASTAs; the TEC transcript has none in either. Sources:
<https://ftp.ensembl.org/pub/release-106/>, under Ensembl's
[data use terms](https://www.ensembl.org/info/about/legal/disclaimer.html).

OpenProt: Leblanc S, et al. Nucleic Acids Res. 2024;52(D1):D522–D528,
doi:10.1093/nar/gkad1050. OpenProt states no licence for its downloads.
