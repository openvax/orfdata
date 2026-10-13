# Design notes

These notes record the library boundaries that the public API follows.

## Ownership

PyEnsembl remains responsible for reference genomes, annotations, transcripts,
exons and canonical proteins. A small separate library adds source-defined ORFs
on those transcripts and their experimental evidence. oncoref remains a consumer
for cancer-specific analysis. Neither Varcode nor Isovar owns reference datasets.
PepData's immune-epitope/property scope is not a fit.

The new library has a dependency-light data layer and a PyEnsembl-compatible
binding API. It must not import Varcode, Isovar or oncoref. No new PyEnsembl core
API is needed to bind a Transcript; avoid reverse dependencies or monkey-patching.

## Public objects

- ORFModel: immutable source ORF identity, assembly, species (NCBI taxon ID and
  scientific name), translation table, parent transcript, spliced coding blocks,
  protein sequence, ORF class. ORF class is distinct from transcript biotype.
- ORFEvidence: ribosome observations and peptide-level MS observations with their
  source, accession, mapping ambiguity, assay and quality. Native evidence does
  not become evidence for a modified product; orthology is not direct evidence.
- ORFRecord: model plus evidence, retaining legacy flat accessors for Isovar.
- ReferenceORF: validated binding of a record to a PyEnsembl Transcript, retaining
  the original transcript (including noncoding biotype), cDNA coding interval,
  coding sequence and actual resolved annotation/version. Multiple ORFs per RNA
  are normal. Annotation lookup never relabels the transcript protein-coding.
- ORFCatalog: explicitly one species/assembly per catalogue. Query by transcript,
  gene, locus, class or evidence. for_transcript(transcript) validates and returns
  ReferenceORFs. Load custom versioned catalogues offline. The initial bundled
  human data is identified as human-only; no implied mouse evidence.
- Proteome: annotated proteins of one PyEnsembl genome, with opt-in layers for
  evidence-selected catalogue ORFs and host-restricted viral proteins. Every
  layer records its source and counts; unbound ORFs are reported, not dropped.
- Varcode ORFAnnotation: reference ORF plus effect and realized products. Exports
  separately identify reference model/evidence and modifications, preserving
  phase/splice alternatives. A private transcript-engine adapter permits reuse
  of existing variant logic without exposing a fake transcript as the base ORF.
- Isovar: attach those same annotations to RNA products with candidate-specific
  RNA support, maintaining conditional translation and unresolved RNA polarity.

## Generality and compatibility

Support uORF, overlapping uORF, internal/out-of-frame, downstream and noncoding-
transcript ORFs with the same model and evidence schema. Preserve source classes
verbatim rather than equating them with Ensembl transcript biotypes. Complete
stopped ORFs are the initial binding contract; partial models are rejected
explicitly rather than given invented boundaries.

A generic v1 schema requires species, assembly, translation table and source
provenance. Legacy Isovar human data may be read only through its explicit known
human/GRCh38 provenance. Keep legacy Isovar uorf_* imports/call shapes working.
New generic APIs use orf terminology; defaults remain explicit about dataset scope.

Multi-species path: one versioned catalogue per organism/assembly, importer input
using the same model/evidence schema, optional organism-specific source adapters
(e.g. OpenProt). Test a mouse noncoding transcript with real PyEnsembl objects,
plus human coding transcripts, different frames, minus strand, splice junctions,
assembly/species mismatches and unknown translation tables. Synthetic examples
are identified as test fixtures, never scientific observations. OpenProt already
publishes multiple species with Ensembl/RefSeq mappings and evidence filters.

Primary sources:

- <https://www.nature.com/articles/s41587-022-01369-0>
- <https://www.nature.com/articles/s41586-026-10459-x>
- <https://archive.openprot.org/p/help>
- <https://github.com/openvax/pyensembl>
