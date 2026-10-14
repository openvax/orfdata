# Examples

Real cases from human Ensembl 115, from ordinary upstream ORFs to circular RNAs
and processed precursors. The examples continue one Python session:

```python
from pyensembl import EnsemblRelease
from orfdata import ProteinProduct, SourceRNA, TranscriptORF, annotated_orfs, extra_orfs, rfam_ires

genome = EnsemblRelease(115, species="human")
orfs = extra_orfs(genome)  # Every bundled catalogue ORF on every transcript carrying it.
```

Coordinates are Python slices: `transcript.sequence[orf.start:orf.end]` is the
coding sequence, including its stop codon, and protein products slice the
protein sequence.

## ORF categories

Each catalogue ORF keeps its source category in `orf_type`. `cds_relation` is
computed against the transcript's annotated CDS in the annotation you loaded:

```python
examples = ["c10norep165", "c10norep53", "c10norep101", "c11riboseqorf10", "c11norep16", "c10norep110"]
for orf_id in examples:
    orf = orfs.for_orf(orf_id)[0]  # The transcript the source names comes first.
    print(orf.orf_type, orf.gene_name, orf.transcript_id, orf.start, orf.end, orf.cds_relation)
```

```text
uORF FGFR2 ENST00000457416 96 198 upstream
uoORF EPC1 ENST00000375110 151 367 upstream_overlapping
intORF ANAPC16 ENST00000299381 296 365 internal
dORF TIMM10B ENST00000254616 929 1259 downstream
doORF PHLDA2 ENST00000314222 511 640 downstream_overlapping
lncRNA NUTM2B-AS1 ENST00000483080 88 199 no_cds
```

On their source transcripts in Ensembl 115, the two agree for over 98% of
bundled ORFs:

| `orf_type` (source) | `cds_relation` (Ensembl 115) | ORFs | Other relations |
| --- | --- | ---: | ---: |
| `uORF` | `upstream` | 2,972 | 11 |
| `uoORF` | `upstream_overlapping` | 668 | 4 |
| `intORF` | `internal` | 712 | 3 |
| `dORF` | `downstream` | 487 | 11 |
| `doORF` | `downstream_overlapping` | 60 | 1 |
| `lncRNA` | `no_cds` | 1,892 | 4 |
| `processed_transcript` | `no_cds` | 278 | 5 |

The exceptions mostly follow annotation changes since Ensembl 101, which the
catalogue was built on. Keep both fields: the category is the source's, the
relation describes the annotation you are using.

## Two proteins on one RNA

`annotated_orfs(transcript)` lists every annotated CDS a transcript carries,
including other genes' proteins. MIURF's CDS lies in the "5′ UTR" of MIEF1's
mRNA, and SNURF's in SNRPN's:

```python
for transcript_id in ("ENST00000325301", "ENST00000390687"):  # MIEF1, SNRPN
    for orf in annotated_orfs(genome.transcript_by_id(transcript_id)):
        print(orf.gene_name, orf.source_gene_name, orf.orf_id, orf.start, orf.end, orf.cds_relation)
```

```text
MIEF1 MIURF ENSP00000490747.1 176 389 upstream
MIEF1 MIEF1 ENSP00000390029.1 487 615 C_terminal_truncation
MIEF1 MIEF1 ENSP00000327124.2 487 1879 annotated_cds
SNRPN SNURF ENSP00000342152.4 62 278 upstream
SNRPN SNRPN ENSP00000462211.1 466 644 C_terminal_truncation
SNRPN SNRPN ENSP00000463510.1 466 974 C_terminal_truncation
SNRPN SNRPN ENSP00000375105.4 466 1189 annotated_cds
```

Each protein keeps its own ID and gene; `gene_name` is the RNA's gene. The
`C_terminal_truncation` rows are shorter CDSs of other MIEF1 and SNRPN
transcripts whose coding paths these RNAs also carry. Identical placements are
reported once, with the other protein IDs in `equivalent_ids`.

## Starts other than AUG

Ensembl annotates some CDSs from CUG, GUG or AUU codons. Ensembl gives GUG and
AUU starts no `start_codon` row; orfdata still reports them as complete, because
the CDS begins inside its RNA:

```python
for gene, transcript_id in (("EIF4G2", "ENST00000339995"), ("TEAD1", "ENST00000527636"), ("MYC", "ENST00000621592")):
    cds = TranscriptORF.from_annotation(genome.transcript_by_id(transcript_id))
    print(gene, cds.start_codon, cds.start_complete, cds.protein_sequence[:5], len(cds.protein_sequence))
```

```text
EIF4G2 GTG True MESAI 907
TEAD1 ATT True MEPSS 426
MYC CTG True MDFFR 454
```

Isoforms initiated upstream of an annotated start can be read from the RNA with
`TranscriptORF.at`. FGF2's annotated CDS starts at an AUG; UniProt P09038
documents CUG-initiated isoforms that extend it:

```python
fgf2 = genome.transcript_by_id("ENST00000644866")
for start in (179, 221, 344):
    orf = TranscriptORF.at(fgf2, start, orf_id=f"FGF2@{start}", source="UniProt P09038")
    print(start, orf.start_codon, len(orf.protein_sequence), orf.end, orf.cds_relation)
```

```text
179 CTG 210 812 N_terminal_extension
221 CTG 196 812 N_terminal_extension
344 ATG 155 812 annotated_cds
```

A start outside the genetic code's initiation codons must be declared, as in
`TranscriptORF.at(rna, start, start_codon="GTG", ...)`; otherwise the read is
rejected. Catalogue models declare it with `ORFModel(start_codon=...)`.

## Recoded stop codons and frameshifts

`decoded_stop_codons` lists internal stop codons that the protein reads
through. GPX1 has a selenocysteine:

```python
gpx1 = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000419783"))
print(gpx1.decoded_stop_codons)
```

```text
((48, 'TGA', 'U'),)
```

Stop-codon readthrough can be described from its source. Eswarappa et al. (Cell
2014) report VEGFA's UGA stop decoded as serine, giving the 22-residue extension
of VEGF-Ax:

```python
vegfa = genome.transcript_by_id("ENST00000672860")
cds = TranscriptORF.from_annotation(vegfa)
ax = TranscriptORF.at(vegfa, cds.start, orf_id="VEGF-Ax", source="Eswarappa et al. 2014",
                      recoded={len(cds.protein_sequence): "S"})
print(ax.protein_sequence[len(cds.protein_sequence):], ax.cds_relation, ax.decoded_stop_codons)
```

```text
SAGQEEGASLRVSGTRSLTRKD C_terminal_extension ((395, 'TGA', 'S'),)
```

Ensembl models programmed ribosomal frameshifts as one- or two-base gaps in the
coding path. `subcodon_introns` reports them with the reading-frame shift:

```python
oaz1 = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000602676"))
print(oaz1.subcodon_introns, len(oaz1.protein_sequence))
```

```text
((204, 2271443, 2271443, 1),) 228
```

The skipped base is the T of the in-frame UGA after UCC, so antizyme 1 is read
in the +1 frame from there. Some such gaps instead correct reference-sequence
errors, so check the gene before interpreting one.

## Endogenous retrovirus proteins

Proteins from endogenous retroviruses that Ensembl annotates, such as
syncytin-1 (ERVW-1) and syncytin-2 (ERVFRD-1), are part of the annotated
proteome. Some catalogue ORFs have since become annotated proteins: the
catalogue's lncRNA ORF `c21norep46` is now ERVH48-1's CDS, which encodes
suppressyn:

```python
erv = orfs.for_orf("c21norep46")[0]
print(erv.orf_type, erv.gene_name, erv.cds_relation, len(erv.protein_sequence), erv.evidence.has_evidence("ms"))
```

```text
lncRNA ERVH48-1 annotated_cds 160 True
```

Translated ERV sequences without an annotated transcript need genomic ORF
models that orfdata does not yet support.

## Circular RNAs

`SourceRNA` represents an RNA that is not an annotated transcript, without
inventing an Ensembl ID for it. circ-ZNF609 is ZNF609's second exon
back-spliced into a circle (Legnini et al., Mol Cell 2017). Its ORF starts at
the linear mRNA's AUG and reads across the junction:

```python
znf609 = genome.transcript_by_id("ENST00000326648")
circ = SourceRNA.from_transcript(znf609, 261, 1135, id="circZNF609", source="Legnini et al. 2017",
                                 topology="circular")
orf = TranscriptORF.at(circ, 127, orf_id="circZNF609-ORF", source="Legnini et al. 2017")
print(len(circ.sequence), orf.segments, len(orf.coding_sequence), len(orf.protein_sequence), orf.coding_blocks)
```

```text
874 ((127, 874), (0, 6)) 753 250 ((64499293, 64499298), (64499420, 64500166))
```

The 753-nt ORF crosses the junction and stops 6 bases after it, as the paper
describes. The stop codon read after the junction is not in the linear mRNA, so
the circular protein differs from ZNF609 after residue 249.

## Processed precursors

`ProteinProduct` slices a precursor's protein. POMC is cleaved into several
hormones; positions here are from UniProt P01189 (entry version 239), whose
1-based ranges become Python slices. `products_from_uniprot(entry)` converts a
whole UniProt JSON entry.

```python
pomc = TranscriptORF.from_annotation(genome.transcript_by_id("ENST00000395826"))
products = [
    ProteinProduct("ACTH", "Corticotropin", "peptide", 137, 176, "UniProt P01189 v239"),
    ProteinProduct("alpha-MSH", "Melanotropin alpha", "peptide", 137, 150, "UniProt P01189 v239"),
    ProteinProduct("beta-endorphin", "Beta-endorphin", "peptide", 236, 267, "UniProt P01189 v239"),
]
for product in products:
    print(product.product_id, product.sequence(pomc)[:13], product.coding_slice(pomc), product.genomic_blocks(pomc))
```

```text
ACTH SYSMEHFRWGKPV (514, 631) ((25161357, 25161473),)
alpha-MSH SYSMEHFRWGKPV (514, 553) ((25161435, 25161473),)
beta-endorphin YGGFMTSEKSQTP (811, 904) ((25161084, 25161176),)
```

`coding_slice` gives the product's codons in the transcript sequence and
`genomic_blocks` their chromosome positions, so a variant in a precursor can be
placed in each mature product.

## IRES elements

`rfam_ires(rna)` finds the human cellular internal ribosome entry site (IRES)
sequences of Rfam 15.1 that occur exactly in an RNA:

```python
for gene, transcript_id in (("MYC", "ENST00000621592"), ("VEGFA", "ENST00000672860"),
                            ("NTRK2", "ENST00000277120"), ("ODC1", "ENST00000234111")):
    transcript = genome.transcript_by_id(transcript_id)
    cds = TranscriptORF.from_annotation(transcript)
    for element in rfam_ires(transcript):
        print(gene, element.feature_id, element.start, element.end, "CDS starts", cds.start, cds.start_codon)
```

```text
MYC RF00216:X00364.2/2502-2803 13 315 CDS starts 363 CTG
VEGFA RF00461:M63971.1/3107-3415 742 1051 CDS starts 496 CTG
NTRK2 RF00547:AF410902.1/3685-4037 179 532 CDS starts 532 ATG
ODC1 RF02535:AK292352.1/308-434 207 334 CDS starts 334 ATG
```

An IRES is an RNA feature, not an ORF category; an element can lie upstream of
a CDS or inside it. The NTRK2 and ODC1 elements end exactly at the AUG, and the
MYC element lies upstream of even the CUG-initiated CDS. The VEGFA element was
defined upstream of the AUG-initiated protein, so it lies inside the longer
CUG-initiated CDS that Ensembl annotates. The bundled families are those Rfam curates from published reports,
and several cellular IRESs are disputed: dicistronic reporter assays can be
confounded by cryptic promoters and splicing. A match locates a reported IRES
sequence; it does not establish activity in that transcript or condition.

## Variants

For annotating a variant, look up the ORFs and features at its position:

```python
print([(o.orf_id, o.transcript_id, o.cds_relation) for o in orfs.overlapping("7", 66205486, 66205486, genome_build="GRCh38")][:3])
myc_ires, = rfam_ires(genome.transcript_by_id("ENST00000621592"))
print(myc_ires.genomic_blocks, myc_ires.overlaps("8", 127736300, 127736300))
```

```text
[('c7riboseqorf56', 'ENST00000304842', 'upstream'), ('c7riboseqorf56', 'ENST00000904436', 'upstream'), ('c7riboseqorf56', 'ENST00000904438', 'upstream')]
((127736244, 127736545),) True
```

`overlapping` searches coding paths, optionally with initiation flanks; an RNA
feature reports overlap through its genomic blocks.
