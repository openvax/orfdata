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

`load_catalog()` loads the bundled human reference data offline. Catalogue
access uses only the Python standard library; the optional `pyensembl` extra
supplies PyEnsembl. Binding uses a transcript you have already loaded.

`catalog.for_transcript(transcript)` returns a `ReferenceORF` for each matching
ORF. Each has:

- `model`: source identity, ORF type, coordinates, species and protein sequence.
- `evidence`: ribosome observations, peptide-level MS and native aggregate support.
- `transcript`: the original PyEnsembl transcript, including its original biotype.
- `coding_sequence` and `protein_sequence`: checked against that transcript.

Binding checks the assembly, available species information, transcript identity,
full spliced coding path, protein sequence and terminal stop. A matching gene
name alone is insufficient. Incompatible annotations raise an error.

`extra_orfs(genome)` does the same for every transcript of a genome at once and
reports ORFs that do not bind. To inspect source records without loading a
transcript:

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

The host is one of possibly several transcripts that carry the same ORF. In
Ensembl 115, 5,830 of the 7,108 bound ORFs have an identical spliced coding
path on at least one other transcript. orfdata binds only the source's host.

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

The current binder accepts complete ORFs with a terminal stop under NCBI
tables 1, 2 or 3 and their initiation codons; see
[genetic codes](species.md#genetic-codes). All bundled models use table 1 and start
at ATG. Other genetic codes, initiation conventions and incomplete ORFs need
additional support. Transcript-boundary contexts can be partial; initiation efficiency and
effects on the main coding region are not established by these data.

Circular-RNA translation, stop-codon readthrough and translation from new
transcripts or fusions are not covered by this reference catalogue. Those need
additional source data and, where necessary, model support. ORF discovery is
outside the current API.

`load_uorf_catalog()` preserves Isovar's upstream-only snapshot: 3,771 models,
dataset `2026-10-08.2`. New code should use the broader `load_catalog()` API.

## Other species

[ReferenceContext and ScopedTarget](species.md) keep taxa, strains, assemblies
and annotation versions separate. Custom genomes without species/taxonomy need
an explicit context when binding; source-pinned contexts require an exact match.

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
See the [design notes](design.md) for the model and library boundaries.
