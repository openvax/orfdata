# orfdata

**An evidence layer over [PyEnsembl](https://github.com/openvax/pyensembl) for
proteins beyond standard annotations.** orfdata tells you which open reading
frames (ORFs) have experimental support for producing peptides or proteins.

PyEnsembl supplies reference transcripts and annotated proteins. orfdata adds
other ORFs, their sequences, and evidence to those same transcript objects.
One RNA can have several ORFs, even when its annotation says it is noncoding.
The original transcript's biotype and annotated protein stay intact.

## Which ORFs?

The bundled catalogue includes:

- **Upstream ORFs (uORFs):** before the main protein-coding region.
- **Overlapping ORFs:** extending into or out of the main coding region.
- **Internal ORFs:** alternative reading frames within the main coding region.
- **Downstream ORFs (dORFs):** after the main coding region.
- **ORFs on noncoding RNAs:** including long noncoding RNAs (lncRNAs) and
  processed transcripts without an annotated protein.

These are additional possible protein sources that standard transcript protein
annotations may omit. The catalogue currently contains **7,264 human GRCh38
ORFs**, of which **2,186 have reported supporting mass-spectrometry observations**.
It is a source-backed reference catalogue; it does not discover every possible
ORF. [Data sources and current limits](docs/reference.md#data-and-limits).

## What does the evidence mean?

| Evidence | What it supports |
| --- | --- |
| Ribosome profiling | Translation of the ORF. |
| Protease-digest mass spectrometry (MS) | Detection of peptides matching its protein sequence. |
| MHC/HLA peptidomics | Detection of matching peptides presented by MHC molecules. |

Peptide matches can be ambiguous, so orfdata keeps their source, mapping and
quality information. These observations have different strengths: a translated
ORF or detected peptide does not by itself establish a stable, functional protein.
Missing evidence also does not establish that an ORF is never translated.
[Evidence interpretation](docs/reference.md#evidence).

## Use with PyEnsembl

With the selected PyEnsembl annotation already downloaded and indexed:

```python
from pyensembl import EnsemblRelease
from orfdata import load_catalog

genome = EnsemblRelease(95, species="human")
catalog = load_catalog()
transcript = genome.transcript_by_id("ENST00000304842")  # TPST1

for orf in catalog.for_transcript(transcript):
    print(orf.model.orf_id, orf.protein_sequence)
    print(orf.evidence.has_evidence("ms"))
    assert orf.transcript is transcript
```

The same call works for a noncoding parent, such as NUTM2B-AS1
(`ENST00000483080`). Each returned `ReferenceORF` keeps its model, evidence and
original PyEnsembl transcript together, after checking its coordinates and
sequence against that transcript.

Varcode describes changes to this reference ORF; Isovar connects it to observed
RNA products. Reference evidence stays separate from evidence for a variant or
fusion. See the [API and integration guide](docs/reference.md).

The catalogue API supports other species, with explicit species and assembly
checks. **The bundled experimental ORF catalogue is currently human-only**. Offline tests
bind real human, mouse, dog and yeast reference slices; those fixtures are not
new experimental ORF catalogues. [Species and reference identity](docs/species.md).

[Optional catalogue management](docs/catalogs.md) provides explicit, versioned
downloads and source builds while returning the same ORF objects and queries.

Install with `pip install orfdata[pyensembl]`. See the
[measurement API](docs/measurements.md) for native units, explicit samples and
assays, and DataCache-managed downloads. Downloadable
[human, mouse and dog protein abundance](docs/pfizer-protein.md) and
[yeast calibrated protein/RNA abundance](docs/copy-yeast.md) use the same
species-scoped API, as do [native dog and mouse RNA counts](docs/pfizer-rna.md)
and [matched human tissue RNA/protein measurements](docs/human-tissue.md).
[Yeast RNA and ribosome occupancy](docs/ingolia-yeast.md) retain their source
regions, including UTR and ncRNA labels.
Gene RNA and canonical protein measurements do not establish
noncanonical ORF abundance. The [native search references](docs/search-references.md)
preserve complete deposited sequences, including unreported entries. The
[quantitative, multispecies data roadmap](https://github.com/openvax/orfdata/issues/3)
tracks further RNA abundance, ribosome occupancy, absolute quantification and
additional protein data releases. See the [released coverage and gaps](docs/coverage.md).
