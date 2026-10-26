# Pfizer dog Proteome Discoverer source slices

`native.json` contains exact selected rows/columns retrieved on 2026-10-26 from
the public MSV000099744 / PXD070290 deposit's
`search/search/IQP_1269_Dog_P1toP4.pdResult` (86,023,450,624 bytes;
MDTM 20251028233329), using bounded explicit-FTPS range reads through SQLite.
The database's full-file checksum is unavailable. The fixture records the
source identity and schema fingerprint, not a full-file integrity claim.

Tables and columns are named as deposited. The samples include the first three
target PSMs/proteins/peptide groups; three quantified and two shared peptide
groups from the first 100 row IDs; the corresponding selected distribution
metadata; and all 97 WorkflowInputFiles rows (96 raw files and one processing
MSF). Numeric BLOBs retain their original bytes in base64 with a byte length
and SHA-256. Other metadata rows are retained exactly.

These factual slices are not a full extraction or a consumer release. The
fixture's candidate examples do not establish unique protein attribution,
individual-dog expression or DLA presentation. The full searched FASTA contains
more exact sequence candidates for peptide group 8 than the vendor's reported
protein list; a complete importer must retain those distinct facts.

Attribution: Oziolor et al., Joel D. Federspiel and Pfizer. MassIVE dataset terms:
CC0-1.0. Source article: https://doi.org/10.1038/s41597-026-07512-9.
Deposit: https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=062b109e9adf461ba619369dd0d85aa9.

The numeric value/presence storage convention is independently documented in
https://github.com/BenBruyneel/proteinDiscover/blob/main/R/transformRaw.R.
No implementation from that project is included here; native type/label
metadata and exact source bytes govern the decoding tests.
