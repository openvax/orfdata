# Reference binding fixtures

These are small **canonical** transcript/sequence slices for testing the same
binding API across species. They contain no newly imported noncanonical ORF
or mass-spectrometry evidence.

Human TP53, mouse Trp53 and dog TP53 come from Ensembl REST records captured on
2026-10-12. Their annotation provider is `ensembl-rest` and their version is the
snapshot date: a live REST response is not assigned an unverified numbered
Ensembl release. Yeast ACT1 comes from the pinned Ensembl 115 R64-1-1 GTF/cDNA
and protein FASTAs, which import SGD annotations for S288C. Each fixture lists
source URLs, byte counts and SHA-256 hashes. The protein sequence is obtained
from the source protein FASTA/API independently of translation in the test.

Tests construct transcript/exon-only GTF slices and use actual PyEnsembl Genome
objects to check full spliced coding paths and translation. Source protein
annotations are omitted from these small GTFs; that omission is not evidence
that these genes lack annotated proteins. Three fixtures are minus-strand;
mouse is plus-strand. Yeast retains its SGD identifiers and strain taxon 559292
separately from species taxon 4932.

Sources: [Ensembl data use](https://www.ensembl.org/info/about/legal/disclaimer.html),
[yeast assembly/annotation](https://jun2026.archive.ensembl.org/Saccharomyces_cerevisiae/Info/Annotation),
[NCBI S288C taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=559292).
