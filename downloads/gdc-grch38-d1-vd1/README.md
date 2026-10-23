# GRCh38.d1.vd1 viral reference bundle

Downloaded from NCBI GenBank on 2026-10-22 using the exact GenBank accession in each row of the GDC GRCh38.d1.vd1 viral inventory. The accompanying manifest retains any RefSeq ID listed by GDC.

Files:

- `genomes.fasta`: one nucleotide reference per panel row.
- `coding_sequences.fasta`: nucleotide sequences extracted from source CDS features.
- `proteins.fasta`: source translations of those CDS features.
- `cds.tsv`: CDS feature metadata, source locations and sequence lengths.
- `manifest.tsv`: source identifiers, resolved accession versions, viral taxon IDs, sequence hashes and CDS counts. Source host taxon IDs remain empty when unknown.
- `genbank/`: gzipped GenBank records, one per source reference.

This bundle contains 200 genomes and 1678 paired coding and protein sequences. Genome FASTA IDs are `viral:<accession.version>`; paired CDS/protein IDs add the source CDS ordinal and labels. Collection and context provenance do not alter IDs.

The collection is admitted for human analysis. That target does not assert a source host; `source_host_taxon_ids` is empty when the source inventory does not establish one.

This is a source-record download, not an ORFdata curated viral-reference release. It does not alter `viral_references()` or `Proteome.build(..., viral=...)`.

Source inventory: https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt
NCBI retrieval: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
