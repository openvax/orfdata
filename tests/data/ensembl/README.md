# Ensembl 115 reference slices

Exact source records for RBMY1F, RBMY1J and PRAME from the full human Ensembl 115
GTF, cDNA, ncRNA and peptide files. These small gzip files were originally audited
in OncoRef #649. The manifest distinguishes fixture hashes from full source hashes.
The ncRNA slice is an empty gzip because these 19 selected transcript records are
in the cDNA file. Thirteen have annotated proteins; identical RBMY1F/J protein
sequences retain separate source identities.

Primary sources are under
`https://ftp.ensembl.org/pub/release-115/{gtf,fasta}/homo_sapiens/`.
GTF: `Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz`.
FASTA: `cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`,
`ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz`,
`pep/Homo_sapiens.GRCh38.pep.all.fa.gz`.

Tests label all modified records as synthetic mutations of these fixtures.
Reference membership does not establish expression or empirical translation.

`internal-stop/` contains one further exact source transcript/protein record,
ENSP00000488168.1, whose six-residue annotated protein contains two `*` symbols.
The reader must preserve the native sequence rather than discard or repair it.
