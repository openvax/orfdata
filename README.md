[![Tests](https://github.com/openvax/orfdata/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/orfdata/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/orfdata.svg)](https://pypi.org/project/orfdata/)

# orfdata

orfdata builds the protein reference for a species: the proteins annotated in
a [PyEnsembl](https://github.com/openvax/pyensembl) genome, plus noncanonical
open reading frames (ORFs) with experimental evidence of translation. These
include upstream and downstream ORFs, ORFs in alternative reading frames, and
ORFs on RNAs annotated as noncoding. Each added ORF stays bound to its
PyEnsembl transcript, with its source and evidence.

## Installation

Use Python 3.9 or later. Install orfdata and one PyEnsembl annotation:

```sh
python -m pip install 'orfdata[pyensembl]'
pyensembl install --release 115 --species human
```

## Example

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

genome = EnsemblRelease(115, species="human")
proteome = Proteome.from_genome(genome).with_extra_orfs(evidence="ms")
proteome.write_fasta("human-e115.fasta")
print(proteome.summary()["kinds"], len(proteome.excluded))
```

```text
{'annotated': 245535, 'extra_orf': 2140} 46
```

This takes every protein in the Ensembl 115 annotation and adds the 2,140
bundled human ORFs that have supporting mass-spectrometry peptides and bind
exactly to Ensembl 115 transcripts. The other 46 selected ORFs are listed in
`proteome.excluded` with the reason. Use `evidence="ribosome"` for all 7,108
ribosome-profiling ORFs that bind. From the command line:

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
```

A translated ORF or a detected peptide does not by itself establish a stable,
functional protein. Most MS support in the bundled catalogue comes from
peptides presented by HLA molecules.

## Species

The bundled noncanonical ORF catalogue is human GRCh38 only. Mouse, dog and
yeast proteomes contain their PyEnsembl-annotated proteins; `with_extra_orfs()`
raises for them until a qualified catalogue is released. Human evidence is never
transferred to another species. Viral proteomes can be added with
`with_viral(...)`, but only to a proteome for one of the virus's recorded hosts.
orfdata also provides [RNA, ribosome and protein measurements](https://openvax.github.io/orfdata/measurements/)
for human, mouse, dog and yeast.

## Documentation

The [documentation](https://openvax.github.io/orfdata/) covers:

| Task | Guide |
| --- | --- |
| Build a proteome: evidence levels, exclusions, contigs, viral proteins | [Build a reference proteome](https://openvax.github.io/orfdata/proteome/) |
| ORF classes, evidence definitions, Varcode and Isovar integration | [ORFs and evidence](https://openvax.github.io/orfdata/reference/) |
| Species and reference identity, mitochondrial genetic codes | [Species, references and genetic codes](https://openvax.github.io/orfdata/species/) |
| Versioned ORF catalogues managed through DataCache | [Optional ORF catalogues](https://openvax.github.io/orfdata/catalogs/) |
| Quantitative datasets and what has been released | [Measurements](https://openvax.github.io/orfdata/measurements/), [released data and gaps](https://openvax.github.io/orfdata/coverage/) |

The [data roadmap](https://github.com/openvax/orfdata/issues/3) tracks further
species and sources. Every PR includes a version bump, lint and tests;
maintainers publish with `./deploy.sh` from clean `main`.
