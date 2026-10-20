# Viral coding references

Browse the packaged viral references by host, convenience group and virus type,
then load their source-annotated coding paths without fetching the genome:

```python
from orfdata import viral_orfs, viral_reference, viral_references

references = viral_references(host="human", group="HPV", virus_type="HPV16")
print([(r.id, r.coding_models) for r in references])
# [('NC_001526.4', 9), ('GRCh38.d1.vd1:HPV16', 9)]

ref = viral_reference("NC_001526.4")  # choose the native RefSeq placement
reference_dna = ref.dna
for orf in ref.orfs:
    print(orf.orf_id, orf.gene_name, orf.coding_blocks, orf.protein_sequence)
```

The small bundled catalogue loads offline without PyEnsembl, Biopython or a
human genome. It supplies the complete DNA for the native RefSeq NC_001526.4
and GDC HPV16 placements; both map to the same nine source CDSs. The two DNA
strings differ by a verified circular rotation. `viral_reference_sequence(id)`
returns the exact forward FASTA bases, with T rather than U. A missing ID raises
`ValueError`. `viral_references()` returns immutable `ViralReference` data
classes; `.dna` and `.orfs` load only on access, so browsing does not copy
genome bases or construct ORF models for every result. Choose by reference ID,
not by result position: the native and GDC sequences have different coordinates.
Use `ref.to_dict()` when serializing reference metadata to JSON; it does not
load DNA or ORF models.

For a combined protein FASTA, callers can pass the virus type or exact
reference ID directly to `Proteome.build` without constructing a separate
viral input object:

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

genome = EnsemblRelease(115, species="human")  # Already installed with pyensembl install.
combined = Proteome.build(genome, nc_orfs="all", viral="HPV16")
combined.write_fasta("human-hpv16.fasta")
```

`viral="all"` includes each distinct bundled viral protein set associated
with the host. When equivalent references differ only in DNA coordinate
placement, it chooses the native source ID. The viral FASTA header records
the chosen exact reference ID, source gene and product label; the snapshot
keeps its source checksum and host-range attribution. See
[Build a reference proteome](proteome.md#viral-proteomes).

`host` accepts a known species name such as `"human"` or an NCBI taxon ID such
as `9606`. The released HPV16 reference objects have viral `taxon_id=333760`
and `host_taxon_ids=(9606,)`; these are different identities. Their documented
human host association is sourced to PaVE. `ref.group == "HPV"`
is a convenient selection label, not a formal taxonomic rank. The
`ref.virus_type == "HPV16"` label is the virus
type; ICTV places it in the species *Alphapapillomavirus 9*. Each reference also
retains the association and taxonomy source URLs. Host associations can cover
several host taxa as new references are curated; an unlisted host is unknown,
not proof that infection is impossible. NCBI Virus's isolation-host metadata is
not an exhaustive host-range list. `ref.host_taxon_ids` and `ref.host_source` have the
same meaning as in [ViralProteome](proteome.md#viral-proteomes).

When annotating an existing BAM, continue to bind against its **actual** full
alignment-reference DNA:

```python
orfs = viral_orfs(alignment_reference_dna, contig="HPV16")
```

This checks the sequence hash and its qualified coordinates. A contig name or
related strain cannot transfer CDS annotations to mismatching DNA. It raises
`ValueError` for an unavailable reference rather than returning an empty ORF
set. `viral_references()` reports coverage and source checksums, and
`viral_orfs()` returns the existing `SourceRNA`/`TranscriptORF` models.

Each model uses a `SourceRNA` containing only the available CDS, including
ordered splice blocks. Full transcript ends are not asserted. The existing
`TranscriptORF.at()` validates its source protein and coding boundaries.
These models have `orf_type="reference_CDS"`. Their internal kind `extra_orf`
describes a source-defined ORF outside a loaded PyEnsembl annotation, not a
claim that a RefSeq CDS is biologically noncanonical. No Ribo-seq or MS
support is inferred from the reference annotation.

The nine annotations are not an exhaustive HPV16 ORF catalogue. Other HPV
types/lineages, EBV/KSHV, CMV, hepatitis viruses including HDV, HTLV-1/HIV,
and polyomaviruses remain explicit coverage work in
[issues #91–98](https://github.com/openvax/orfdata/issues/91).

Isovar owns reconstructed sample RNA and read evidence. It can use these
coding paths to anchor a frame on observed sequence, preserving variants
and incomplete starts/stops. ORFdata does not assemble sample RNA, replace
missing bases, infer protein abundance or establish viral tumor causality.
RNA and coding-reference matches are not direct protein evidence.

Source attribution and rebuilding instructions are in
[`orfdata/data/viral/NOTICE.md`](https://github.com/openvax/orfdata/blob/main/orfdata/data/viral/NOTICE.md).
