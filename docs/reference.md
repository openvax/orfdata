# ORFs and evidence

The bundled human catalogue is the
[GENCODE Ribo-seq ORF](https://www.gencodegenes.org/pages/riboseq_orfs/)
Phase 1 set: 7,264 GRCh38 ORFs that ribosome-profiling studies found translated
but that the reference annotation (GENCODE v35) did not annotate as proteins.
2,186 also have supporting MS peptides. "Noncanonical" is relative to an
annotation and transcript: in Ensembl 115, 34 of the 7,108 ORFs that bind encode
a protein identical to an annotated ENSP, and an ORF can be annotated in one
transcript model but not another. The catalogue is source-backed; it does not
discover every possible ORF. The larger Phase 2 set is planned
([#16](https://github.com/openvax/orfdata/issues/16)). Most users start with a
[reference proteome](proteome.md) or `extra_orfs(genome)`; this page describes
the objects and evidence underneath.

## ORF categories

Each ORF keeps its source category verbatim as `orf_type`. The terms follow
[Mudge et al. 2022, Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757701/#T2):

| `orf_type` | Position |
| --- | --- |
| `uORF` | Within the 5′ UTR of a protein-coding gene's transcript |
| `uoORF` | Begins in the 5′ UTR and overlaps the annotated coding sequence in a different frame |
| `intORF` | Entirely within the annotated coding sequence, in a different frame |
| `dORF` | Within the 3′ UTR |
| `doORF` | Begins within the annotated coding sequence and ends in the 3′ UTR, in a different frame |
| `lncRNA` | On a transcript annotated as a long noncoding RNA |
| `processed_transcript` | On a transcript annotated as `processed_transcript`, with no annotated protein |

The category, the host transcript's biotype and the evidence are separate
properties. The category describes the source's annotation; a bound ORF's
`transcript.biotype` describes the annotation you loaded. In Ensembl 115, for
example, most `processed_transcript` ORFs are on `protein_coding_CDS_not_defined`
or `retained_intron` transcripts, and 48 uORFs are on nonsense-mediated-decay
transcripts. One RNA can carry several ORFs.

## Objects and queries

A `TranscriptORF` is one ORF on one PyEnsembl transcript. It covers both
annotated CDSs and the catalogue's extra ORFs, so all the ORFs on an RNA are
handled the same way. `ReferenceORF` is the same class under its original
name. Each has:

- `kind`: `"annotated"` for an annotated CDS, `"extra_orf"` for a catalogue ORF.
- `orf_id`: the source ORF ID, or the versioned protein ID of an annotated CDS.
- `start` and `end`: Python slice coordinates into the transcript's spliced
  sequence, so `transcript.sequence[orf.start:orf.end] == orf.coding_sequence`.
  The coding sequence includes the stop codon when there is one. The genomic
  path is in `coding_blocks`, 1-based and inclusive.
- `transcript`, `transcript_id`, `gene`, `gene_id` and `gene_name`: the
  transcript carrying the ORF and its gene in the annotation you loaded.
- `source_transcript_id`, `source_gene_id` and `source_gene_name`: what the ORF's
  source says, and `on_source_transcript`: whether this is that transcript.
- `protein_sequence`, `orf_type`, `model` and `evidence`. An annotated CDS has
  no `orf_type`, `model` or evidence of its own.
- `start_complete` and `stop_complete`, and for an annotated CDS
  `equivalent_ids`: other annotated proteins with the identical placement.

An extra ORF is checked strictly: the assembly, available species information,
the full spliced coding path, the protein sequence and the terminal stop must
all match. An annotated CDS is taken as the annotation defines it, with the
annotation's protein. Annotations include selenoproteins and CDSs without
marked start or stop codons, so translation is not re-checked;
`start_complete` and `stop_complete` report what the annotation marks.

```python
from orfdata import extra_orfs

orfs = extra_orfs(genome)  # Every catalogue ORF, on every transcript carrying it.
tpst1 = genome.transcript_by_id("ENST00000304842")
for orf in orfs.for_transcript(tpst1, include_annotated=True):
    print(orf.kind, orf.orf_id, orf.orf_type, orf.start, orf.end, len(orf.equivalent_ids))
```

```text
extra_orf c7riboseqorf56 uORF 149 269 0
annotated ENSP00000302413.5 None 290 1403 21
```

Results come in RNA order: the uORF ends 21 bases before the annotated CDS
starts, and 21 other proteins share TPST1's annotated placement exactly.
Another transcript's CDS is listed only when this RNA encodes the same ORF. A
CDS without a marked start or stop codon, such as a fragment on a truncated
transcript, is listed only when this RNA continues identically beyond that end;
otherwise translation would read on. Two such TPST1 fragments are not listed here. `annotated_orfs(transcript)` gives just the annotated ones;
`TranscriptORF.from_annotation(transcript)` gives one transcript's own CDS.

`orfs.overlapping(contig, start, end, genome_build=...)` finds ORFs whose coding
path overlaps a 1-based inclusive genomic interval; introns between coding
blocks don't count. Add `include_initiation_context=True` to include an ORF's
spliced start flanks, for example for a Kozak-position variant. Flanks come
from the transcript each placement is on (`orf.initiation_context`), since
transcripts sharing a coding path can differ upstream of the start codon.

### Transcripts carrying an ORF

An ORF is placed on every transcript that contains its exact spliced coding path
and translates to its protein, not only on the transcript its source names. In
Ensembl 115 the bundled ORFs have 67,529 such placements; only 1,314 ORFs are
carried by a single transcript. `on_source_transcript` marks the one the source
names. Another carrying transcript is not evidence that the ORF is translated
from it, and the ORF's `evidence` is the same on every placement.

`orfs.source_failures` lists ORFs whose source transcript is absent from this
annotation or no longer carries them, with the reason. Most are still placed
elsewhere: of 156 such ORFs in Ensembl 115, only the 25 in `orfs.excluded` are
carried by no transcript. `record.bind(transcript)` still binds strictly to
the source's transcript.

Annotated CDSs can also lie in other transcripts. MIURF's CDS sits in the
"5′ UTR" of MIEF1's mRNA, for example, and SNURF's in SNRPN's:

```python
mief1 = genome.transcript_by_id("ENST00000325301")
for orf in orfs.for_transcript(mief1, include_annotated=True):
    print(orf.kind, orf.orf_id, orf.source_gene_name, orf.start, orf.end, orf.on_source_transcript)
```

```text
extra_orf c22riboseqorf42 MIEF1 72 132 True
annotated ENSP00000490747.1 MIURF 176 389 False
annotated ENSP00000327124.2 MIEF1 487 1879 True
extra_orf c22norep69 MIEF1 527 677 True
```

MIURF's protein keeps its own ID and gene while sitting on MIEF1's RNA, so it is
not mistaken for a new protein. The [examples](examples.md) cover non-AUG starts,
recoded stop codons, frameshifts, circular RNAs, precursor products and IRES
elements.

### Catalogue records

`load_catalog()` loads the bundled human records offline using only the Python
standard library; like PyEnsembl, it defaults to `species="human"`.
`catalog.for_transcript(transcript)` binds the catalogue ORFs on one
transcript. To inspect source records without a genome:

```python
from orfdata import load_catalog

catalog = load_catalog()
records = catalog.query(gene="TPST1", evidence="ms")
```

## Identifiers

Bundled ORF IDs are the source's GENCODE Phase I IDs. They start with the
chromosome: `c7riboseqorf56` is on chromosome 7. The 3,085 `c…riboseqorf…`
ORFs were detected in two or more ribosome-profiling studies, and the 4,179
`c…norep…` ORFs in one. They are not Ensembl stable IDs.

Each record names one host transcript (`transcript_id`, an unversioned ENST)
and that transcript's gene (`gene_id`, an ENSG), as validated against Ensembl
101. For an ORF on an mRNA, this is the gene of the transcript's annotated
protein. Later releases can change both: of the 7,108 ORFs that bind to Ensembl
115, 87 host transcripts now belong to a different ENSG and 796 gene symbols
changed. A bound ORF's `gene_id`, `gene_name` and `transcript` describe the
genome you loaded; `orf.model` keeps the source's labels. In
`extra_orfs(genome)`, `for_gene` matches genes in the loaded annotation, or the
source's labels with `source=True`. `for_transcript` and `for_gene` reject a
transcript or gene from another reference, and a versioned Ensembl ID must match
the loaded version.

The host is one of possibly several transcripts that carry the same ORF; see
[transcripts carrying an ORF](#transcripts-carrying-an-orf).

## Evidence

Use `orf.evidence.has_evidence("ribosome")`, `"shotgun"` (protease-digest MS),
`"mhc"` (including human HLA), or `"ms"` (any supporting MS observation).
`"hla"` selects human HLA observations specifically. Source observations whose
MS assay type is unspecified remain labelled `"ms_unclassified"`.

```python
orf.evidence.has_evidence("ms", require_unique=True, require_reviewed=True)
```

This requires a supporting spectrum review and a single protein mapping
**within that source's search scope**. It does not imply uniqueness across every
possible protein. Inspect `orf.evidence.peptides` for peptide sequences, source
accessions, matching intervals, mapping counts and review quality. Missing counts
remain unknown; low or mixed spectrum reviews do not count as supporting MS.

`has_evidence` and `catalog.query(evidence=...)` also take rules that set
Ribo-seq and MS stringency separately, such as
`RiboSeq(min_studies=2) & MS(min_peptides=2)`; see
[Build a reference proteome](proteome.md#separate-ribo-seq-and-ms-stringency).

Some sources provide only aggregate support counts. `EvidenceSummary` keeps the
native metric, definition, value, source/release, source checksum and explicitly
scoped entity. For example, [OpenProt's native metrics](https://api.openprot.org/api/2.0/SC/downloads/yeast_s288c-openprot-2_2-altprots+isoforms-min_2_pep-R64-1-1.106.tsv.readme.html)
include an MS score that sums study-specific unique-peptide
counts; it is neither global peptide uniqueness nor protein abundance. Its
translation score counts studies with significant translation detection.
`None` means unknown and remains distinct from a reported zero.

Positive summaries in `orf.evidence.summaries` can satisfy general evidence
queries. To require actual peptide records or named ribosome studies, use:

```python
orf.evidence.has_evidence("ms", require_resolved=True)
catalog.query(evidence="ribosome", require_resolved=True)
```

The CLI uses `--evidence ms --require-resolved`. Uniqueness and review filters
always require resolved peptide observations. Summaries never create peptide
sequences, residue coverage, retained variant peptides or mutant-translation
evidence. Their protein checksum and species/reference checks prevent attaching
them to a different protein or incompatible source context.

A sequence-level summary may have several transcript assignments. Those records
share the summary's `key`; `catalog.evidence_summaries` returns each native metric
once and rejects conflicting facts for the same identity. Native accession
aliases must be resolved explicitly by the source adapter. Counts from different
metrics and sources are not summed, and a later peptide-resolved snapshot is not
an additional independent experiment merely because it has more detail.

This API includes source-backed OpenProt human, mouse and yeast test slices.
It does not yet publish those species' OpenProt ORF catalogues; coordinate and
sequence qualification remains tracked in [#18](https://github.com/openvax/orfdata/issues/18).

The distinction between translation, peptide detection and functional protein
evidence follows the [TransCODE study](https://www.nature.com/articles/s41586-026-10459-x).
Source evidence stays at its declared ORF or protein-sequence scope. It is not
automatically evidence for another isoform, a mutant protein, a fusion or
another species.

## Variants and RNA products

Varcode keeps the reference object alongside its modification:

```python
from pyensembl import EnsemblRelease
from orfdata import load_catalog
from varcode import Variant, annotate_orf_variants

genome = EnsemblRelease(95, species="human")  # already downloaded and indexed
transcript = genome.transcript_by_id("ENST00000304842")
orf, = load_catalog().for_transcript(transcript)
variant = Variant("7", 66205486, "C", "G", genome=genome)
annotation = annotate_orf_variants([variant], orf)
assert annotation.reference is orf
print(annotation.effect)  # P2A in the TPST1 upstream ORF
result = annotation.to_dict()  # reference, evidence and altered products
```

Isovar's `annotate_sv_orf_evidence` attaches the same reference annotations to
reconstructed RNA products and preserves each candidate's RNA support. For
example, the native TPST1 peptide KIILSKISV is absent from the Sid TPST1–CRCP
alternatives, so its MS evidence does not support those altered products.

The reference integration, including strict reference-context validation, is
available in Varcode 10.13.0 and Isovar 1.51.0.

## Data and limits

The bundled human GRCh38 catalogue, dataset `2026-10-11.1`, draws on the
[GENCODE Phase I reference](https://www.nature.com/articles/s41587-022-01369-0),
[TransCODE author exports](https://www.nature.com/articles/s41586-026-10459-x),
and the [community microprotein spectrum benchmark](https://www.nature.com/articles/s41467-025-68002-x).
All 7,264 ORF paths and proteins were checked against Ensembl 101 GTF, coding RNA
and noncoding RNA FASTAs. The 2,186 models with supporting MS observations use
the packaged source-specific filters; this is not a count of proven functional
proteins. Packaged manifests retain source accessions, checksums, licenses,
processing rules and dataset version.

Catalogue models bind as complete ORFs with a terminal stop under NCBI tables
1, 2 or 3 and their initiation codons, or a start codon the model declares with
`start_codon`; see [genetic codes](species.md#genetic-codes). All bundled models
use table 1 and start at ATG. Transcript-boundary contexts can be partial;
initiation efficiency and effects on the main coding region are not established
by these data.

ORFs that a source defines on an RNA rather than in genomic coordinates use
`TranscriptORF.at`: non-AUG starts, fragments without a known start, stop codons
read through, and ORFs on circular or other source RNAs (`SourceRNA`). See the
[examples](examples.md). The bundled catalogue contains none of these. Fusions,
edited RNAs and programmed frameshifts in catalogue models need further model
support, and ORF discovery is outside the current API.

`load_uorf_catalog()` preserves Isovar's upstream-only snapshot: 3,771 models,
dataset `2026-10-08.2`. New code should use the broader `load_catalog()` API.

## Other species

[ReferenceContext and ScopedTarget](species.md) keep taxa, strains, assemblies
and annotation versions separate. Custom genomes without species/taxonomy need
an explicit context when binding; source-pinned contexts require an exact match
for `record.bind(transcript)`. With `any_transcript=True`, as `extra_orfs` uses
to place ORFs by their exact coding path and protein, a record pinned to
another annotation release of the same species, strain and assembly is also
accepted. ORFs on another assembly can be
[lifted over](catalogs.md#lifting-a-catalogue-to-another-assembly).

`load_catalog(species="mouse")` returns a species' catalogue: the bundled one,
or else the one installed through [`CatalogRegistry`](catalogs.md).
`genome_build=` chooses between assemblies. It never downloads, and for a
species without a catalogue the error lists the catalogues that exist.

`load_catalog(path)` reads a custom versioned JSON catalogue. Each catalogue is
scoped to one species and assembly; each model records its NCBI taxon, scientific
name, assembly, translation table and source. Human evidence is never transferred
to a mouse ortholog. Mouse bindings are tested with synthetic ORFs on PyEnsembl
transcripts; there is no bundled mouse evidence.

To build a catalogue, create `ORFModel` objects from source coordinates and
sequences, attach source observations with `ORFRecord.from_model(model, evidence)`,
and construct `ORFCatalog(records, metadata)`. Its `to_dict()` produces the
`orfdata.catalog.v1` JSON representation. Genomic coding blocks are 1-based,
inclusive, sorted by genomic position, and include the terminal stop codon.

[OpenProt](https://archive.openprot.org/p/help) is a possible source for further species,
with organism/assembly-specific models, transcript mappings and experimental
evidence. A future importer must pin releases, preserve ambiguity and licensing,
and validate coordinate conventions. No OpenProt importer is included yet.

## Rebuild and verify

From a development checkout with the PyEnsembl dependency installed:

```sh
python -m pip install -e '.[dev]'
python -m orfdata.build --source-dir /path/to/sources --prepare-sources \
    --all-orfs --output orfdata/data/orf-evidence
./lint.sh
./test.sh
```

Source preparation downloads pinned public inputs and needs Rscript to extract
the author export. Ordinary queries need neither R nor spreadsheet libraries.
