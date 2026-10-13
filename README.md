# orfdata

orfdata is a catalogue of open reading frames (ORFs) that standard protein
annotations leave out, with the experimental evidence that they are translated
or produce detectable peptides. It attaches each ORF, its protein sequence and
its evidence to the corresponding [PyEnsembl](https://github.com/openvax/pyensembl)
transcript, and provides downloadable RNA, ribosome and protein measurement
datasets through the same species-scoped API.

## Install

```sh
pip install 'orfdata[pyensembl]'
```

The bundled catalogue itself uses only the Python standard library; the
`pyensembl` extra is needed to attach ORFs to transcripts.

## Look up ORFs and their evidence

The bundled human GRCh38 catalogue loads offline, with no download:

```python
from orfdata import load_catalog

catalog = load_catalog()
for record in catalog.query(gene="TPST1", evidence="ms"):
    print(record.orf_id, record.transcript_id, record.protein_sequence)
```

To attach ORFs to a PyEnsembl transcript, first download and index that
annotation (`pyensembl install --release 95 --species human`), then:

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

Each returned `ReferenceORF` keeps its model, evidence and original PyEnsembl
transcript together, after checking its coordinates and sequence against that
transcript. The transcript's biotype and annotated protein stay intact, and one
RNA can have several ORFs. The same call works for a noncoding parent, such as
NUTM2B-AS1 (`ENST00000483080`).

## Which ORFs are included?

The bundled catalogue includes:

- Upstream ORFs (uORFs), before the main protein-coding region.
- Overlapping ORFs, extending into or out of the main coding region.
- Internal ORFs, in alternative reading frames within the main coding region.
- Downstream ORFs (dORFs), after the main coding region.
- ORFs on noncoding RNAs, including long noncoding RNAs (lncRNAs) and processed
  transcripts without an annotated protein.

It currently contains 7,264 human GRCh38 ORFs, of which 2,186 have reported
supporting mass-spectrometry observations. It is a source-backed reference
catalogue; it does not discover every possible ORF. See
[data sources and current limits](docs/reference.md#data-and-limits).

## What does the evidence mean?

| Evidence | What it supports |
| --- | --- |
| Ribosome profiling | Translation of the ORF. |
| Protease-digest mass spectrometry (MS) | Detection of peptides matching its protein sequence. |
| MHC/HLA peptidomics | Detection of matching peptides presented by MHC molecules. |

These observations have different strengths. A translated ORF or detected
peptide does not by itself establish a stable, functional protein, and missing
evidence does not establish that an ORF is never translated. Peptide matches can
be ambiguous, so orfdata keeps their source, mapping and quality information.
See [evidence interpretation](docs/reference.md#evidence).

Varcode describes changes to a reference ORF, and Isovar connects it to observed
RNA products. Reference evidence stays separate from evidence for a variant or
fusion; see [variants and RNA products](docs/reference.md#variants-and-rna-products).

## Measurement datasets

RNA, ribosome and protein measurements are separate, downloadable datasets.
Install the `data` extra, then download and query one:

```sh
pip install 'orfdata[data]'
orfdata datasets list
orfdata datasets download ingolia-yeast-rna-ribo
```

```python
from orfdata import DatasetRegistry, ScopedTarget

dataset = DatasetRegistry().open("ingolia-yeast-rna-ribo")
reference, = dataset.manifest.references
target = ScopedTarget(reference, "region", "HRA1", "Ingolia2009:feature")
for measurement in dataset.measurements(target=target, unit="RPKM"):
    print(measurement.sample_id, measurement.quantity, measurement.value)
```

| Guide | Species | Measurements |
| --- | --- | --- |
| [Protein abundance](docs/pfizer-protein.md) | Human, mouse, dog | TMT protein groups (Oziolor et al.) |
| [RNA abundance](docs/pfizer-rna.md) | Dog, mouse | Native gene counts (Oziolor et al.) |
| [Matched human tissues](docs/human-tissue.md) | Human | RNA and protein from adjacent cryosections (Wang et al.) |
| [Yeast protein and RNA](docs/copy-yeast.md) | Budding yeast | Calibrated protein and RNA abundance (CoPY) |
| [Yeast RNA and ribosome occupancy](docs/ingolia-yeast.md) | Budding yeast | Source regions, including UTR and ncRNA labels (Ingolia et al.) |

Gene RNA and canonical protein measurements do not establish noncanonical ORF
abundance. See [released coverage and gaps](docs/coverage.md); the
[data roadmap](https://github.com/openvax/orfdata/issues/3) tracks further RNA
abundance, ribosome occupancy, absolute quantification and protein releases.

## Other species

The catalogue API supports other species, with explicit species and assembly
checks, but the bundled experimental ORF catalogue is currently human-only.
Offline tests bind real human, mouse, dog and yeast reference slices; those
fixtures are not experimental ORF catalogues. See
[species and reference identity](docs/species.md).

## Documentation

- [Using orfdata](docs/reference.md): objects, evidence queries, Varcode and Isovar, data sources and limits.
- [Species and references](docs/species.md): `ReferenceContext`, `ScopedTarget` and binding custom genomes.
- [Optional ORF catalogues](docs/catalogs.md): versioned catalogue downloads and source builds.
- [Measurements](docs/measurements.md): native units, explicit samples and assays, DataCache-managed downloads and importing new sources.
- [Comparing RNA, ribosome and protein](docs/pairing.md): matched assays, derived ratios and orthology.
- [Source references](docs/source-references.md): source identities and explicit mappings to PyEnsembl.
- [Native search references](docs/search-references.md): complete deposited search FASTAs, including unreported entries.
- [Quantitative storage](docs/storage.md): SQLite storage versions and compatibility.

## License

Apache 2.0. Downloadable datasets keep their sources' own terms and
attribution, recorded in each release manifest.
