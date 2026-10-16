[![Tests](https://github.com/openvax/orfdata/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/orfdata/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/orfdata.svg)](https://pypi.org/project/orfdata/)

# orfdata

orfdata builds the protein reference for a species: the proteins annotated in
a [PyEnsembl](https://github.com/openvax/pyensembl) genome, plus open reading
frames (ORFs) that ribosome-profiling studies found translated outside those
annotations. These include upstream and downstream ORFs, ORFs in other reading
frames of a coding sequence, and ORFs on RNAs annotated as noncoding. Each added
ORF stays bound to its PyEnsembl transcript, with its source, category and
evidence.

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
{'annotated': 245535, 'extra_orf': 2181} 5
```

This takes every protein in the Ensembl 115 annotation and adds the 2,181
bundled human ORFs that have supporting mass-spectrometry peptides and are
encoded by an Ensembl 115 transcript. Each entry lists every transcript encoding
it; the other 5 selected ORFs are listed in `proteome.excluded` with the
reason. Use `evidence="ribosome"` for all 7,239 encoded ribosome-profiling ORFs,
or a rule such as `RiboSeq(min_studies=2) & MS()` to set Ribo-seq and MS
stringency separately (1,293 encoded ORFs).
`proteome.save(path)` writes a checksummed snapshot that `Proteome.load` reads
offline. From the command line:

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
```

A translated ORF or a detected peptide does not by itself establish a stable,
functional protein. Most MS support in the bundled human catalogue comes from
immunopeptidomics: peptides eluted from HLA molecules.

## Species

| Species | Noncanonical ORFs | Sources |
| --- | ---: | --- |
| Human | 7,264 bundled; 57,607 with `catalogs="all"` | GENCODE; optionally Chang et al. 2025 and Yang et al. 2024 Ribo-seq |
| Mouse | 48,805 | Chang et al. and Yang et al. Ribo-seq, OpenProt 2.2 MS |
| Zebrafish | 9,759 | Yang et al. Ribo-seq |
| C. elegans | 2,849 | Yang et al. Ribo-seq |
| Yeast | 1,128 | Yang et al. Ribo-seq |
| Dog | none: no public dog Ribo-seq or ncORF MS exists | |

Human's bundled catalogue needs no setup. For the others, `pip install
"orfdata[data,build-data]"`, then `orfdata setup --species mouse` (or
`zebrafish`, `worm`, `yeast`, or `human --catalogs all`) builds the catalogues
from their sources' own files ([details](https://openvax.github.io/orfdata/catalogs/)).
After `orfdata chains download` (add
`--accept-terms` for UCSC's non-commercial mouse chains), older assemblies
(GRCh37, NCBI36, GRCm38, NCBIM37) and T2T-CHM13 use these catalogues lifted
to them; `orfdata genomes install t2t-chm13` installs a CHM13 annotation. Dog proteomes contain their
PyEnsembl-annotated proteins. `load_catalog(species=...)` and
`with_extra_orfs()` pick the catalogue for a species and raise when none is
available or built.
Human evidence is never transferred to another species. Viral proteomes can be added with
`with_viral(...)`, but only to a proteome for one of the virus's recorded hosts.
orfdata also provides [RNA, ribosome and protein measurements](https://openvax.github.io/orfdata/measurements/)
for human, mouse, dog and yeast.

`viral_orfs(reference_dna, contig="HPV16")` loads source-defined viral CDSs as
existing RNA/ORF models, offline and qualified by exact reference sequence.
The first [coding reference](https://openvax.github.io/orfdata/viral/) contains
the nine RefSeq NC_001526.4 HPV16 CDSs, with native and GDC placements;
`viral_references()` lists the available identities and coverage. Broader
imports remain in [the viral roadmap](https://github.com/openvax/orfdata/issues/91).

## Measurement datasets

RNA, ribosome and protein measurements are separate, downloadable datasets with
the same species-scoped API. Install the `data` extra, then download and query
one:

```sh
pip install 'orfdata[data]'
orfdata datasets list
orfdata datasets download ingolia-yeast-rna-ribo
```

| Guide | Species | Measurements |
| --- | --- | --- |
| [Protein abundance](https://openvax.github.io/orfdata/pfizer-protein/) | Human, mouse, dog | TMT protein groups (Oziolor et al.) |
| [RNA abundance](https://openvax.github.io/orfdata/pfizer-rna/) | Dog, mouse | Native gene counts (Oziolor et al.) |
| [Matched human tissues](https://openvax.github.io/orfdata/human-tissue/) | Human | RNA and protein from adjacent cryosections (Wang et al.) |
| [Yeast protein and RNA](https://openvax.github.io/orfdata/copy-yeast/) | Budding yeast | Calibrated protein and RNA abundance (CoPY) |
| [Yeast RNA and ribosome occupancy](https://openvax.github.io/orfdata/ingolia-yeast/) | Budding yeast | Source regions, including UTR and ncRNA labels (Ingolia et al.) |

Gene RNA and canonical protein measurements do not establish noncanonical ORF
abundance. See [released data and gaps](https://openvax.github.io/orfdata/coverage/).

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

## License

Apache 2.0. Downloadable datasets keep their sources' own terms and
attribution, recorded in each release manifest.
