# Using orfdata

See the [README](../README.md) for the purpose, ORF types and basic example.

## Objects and queries

`load_catalog()` loads the bundled reference data offline. Catalogue access uses
only the Python standard library; the optional `pyensembl` extra supplies
PyEnsembl. Binding uses a transcript you have already loaded.

`catalog.for_transcript(transcript)` returns a `ReferenceORF` for each matching
ORF. Each has:

- `model`: source identity, ORF type, coordinates, species and protein sequence.
- `evidence`: ribosome observations and peptide-level MS observations.
- `transcript`: the original PyEnsembl transcript, including its original biotype.
- `coding_sequence` and `protein_sequence`: checked against that transcript.

Binding checks the assembly, available species information, transcript identity,
full spliced coding path, protein sequence and terminal stop. A matching gene
name alone is insufficient. Incompatible annotations raise an error.

To inspect source records without loading a transcript:

```python
from orfdata import load_catalog

catalog = load_catalog()
records = catalog.query(gene="TPST1", evidence="ms")
```

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

The distinction between translation, peptide detection and functional protein
evidence follows the [TransCODE study](https://www.nature.com/articles/s41586-026-10459-x).
Source evidence belongs to the exact reference ORF. It is not automatically
evidence for a different isoform, a mutant protein, a fusion or another species.

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

The current binder accepts complete ORFs with a terminal stop, NCBI translation
table 1, and ATG/CTG/TTG initiation codons. All bundled models start at ATG.
Other genetic codes, initiation conventions and incomplete ORFs need additional
support. Transcript-boundary contexts can be partial; initiation efficiency and
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
