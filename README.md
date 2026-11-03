[![Tests](https://github.com/openvax/orfdata/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/orfdata/actions/workflows/tests.yml)
[![Documentation](https://github.com/openvax/orfdata/actions/workflows/docs.yml/badge.svg)](https://openvax.github.io/orfdata/)
[![PyPI](https://img.shields.io/pypi/v/orfdata.svg)](https://pypi.org/project/orfdata/)

# orfdata

orfdata builds protein references for a species from a
[PyEnsembl](https://github.com/openvax/pyensembl) annotation: the annotated
proteins, plus open reading frames (ORFs) that ribosome-profiling and
mass-spectrometry studies found translated outside those annotations, and viral
proteins when you ask for them. Each added ORF stays bound to its PyEnsembl
transcript, with its source, category and evidence. Its FASTA and GTF outputs
open as PyEnsembl genomes, so [varcode](https://github.com/openvax/varcode),
[isovar](https://github.com/openvax/isovar),
[topiary](https://github.com/openvax/topiary) and
[vaxrank](https://github.com/openvax/vaxrank) can use them. orfdata also
distributes RNA, ribosome and protein measurement datasets for human, mouse,
dog and yeast.

## Install

```sh
python -m pip install orfdata
pyensembl install --release 115 --species human
orfdata setup --species human --catalogs all
```

orfdata needs Python 3.11 or later. The second command downloads and indexes
Ensembl 115 for human GRCh38. The third builds human's optional Ribo-seq ORF
catalogues from their pinned sources (needed for evidence presets such as
`"ms"` and for `"all"`), and downloads the assembly chains and the GDC viral
collection (`--collections none` skips it). The bundled human ORF catalogue
needs no setup. orfdata never downloads data implicitly: every download is an
explicit command, checked against a pinned checksum.

The `canine-extraction` extra (`pip install 'orfdata[canine-extraction]'`)
adds APSW, needed only to extract the original canine vendor database.

## Quick start

Build the human Ensembl 115 proteome plus the noncanonical ORFs (ncORFs) that
have supporting mass-spectrometry peptides, and write it as FASTA:

```python
from pyensembl import EnsemblRelease
from orfdata import ReferenceTranscriptome

genome = EnsemblRelease(115, species="human")
transcriptome = ReferenceTranscriptome.from_genome(genome)
proteome = transcriptome.proteome(nc_orf_sources="ms")
proteome.write_fasta("human-e115.fasta")
print(proteome.summary()["protein_sources"], len(proteome.excluded))
```

```text
{'reference_cds': 245535, 'nc_orf': 2181} 5
```

`proteome()` on its own gives every reference CDS protein, with no ncORFs or
viruses. `nc_orf_sources="ms"` adds the 2,181 MS-supported ORFs whose exact
coding path an Ensembl 115 transcript carries; the other 5 are in
`proteome.excluded` with the reason. Other sources include `"ribosome"`,
`"all"` and named catalogues; `viral="HPV16"` adds a virus. A translated ORF or
a detected peptide does not by itself establish a stable, functional protein:
most MS support in the bundled human catalogue comes from peptides eluted from
HLA molecules.

## Command line

Run `orfdata` for the command guide, or `orfdata help COMMAND` for one
command. Tables are the default output; add `--json` for scripts.

```sh
orfdata data list                       # what can be downloaded or built
orfdata data status                     # what is installed locally
orfdata data fetch ingolia-yeast-rna-ribo
orfdata query --gene TPST1              # ORFs in the installed catalogue
orfdata proteome --species human --release 115 --nc-orf-sources ms --output human-e115.fasta
orfdata viral summary                   # installed viral reference layers
orfdata datasets query ingolia-yeast-rna-ribo --limit 5
```

`download` and `downloads` are aliases for `data`. See the
[command-line guide](https://openvax.github.io/orfdata/cli/) for every command.

## Documentation

The [documentation](https://openvax.github.io/orfdata/) covers, in reading
order:

1. [Getting started](https://openvax.github.io/orfdata/): install, a first
   proteome, ncORF evidence presets and supported species.
2. [Build a reference proteome](https://openvax.github.io/orfdata/proteome/):
   evidence filters, excluded ORFs, viral and sample-derived proteins, FASTA
   headers.
3. [Discover and download data](https://openvax.github.io/orfdata/downloads/):
   the `orfdata data` workflow and the cache.
4. [Downstream tools](https://openvax.github.io/orfdata/downstream/): ncORFs
   and viral references in varcode, isovar, topiary and vaxrank.
5. [ORFs and evidence](https://openvax.github.io/orfdata/reference/): ORF
   classes and what each kind of evidence supports.
6. [Measurements](https://openvax.github.io/orfdata/measurements/): the
   quantitative datasets and their query API.

Release notes are in [CHANGELOG.md](https://github.com/openvax/orfdata/blob/main/CHANGELOG.md);
the [data roadmap](https://github.com/openvax/orfdata/issues/3) tracks further
species and sources.

## Development

```sh
python -m pip install -e '.[dev,docs]'
./lint.sh    # ruff
./test.sh    # pytest
./docs.sh    # mkdocs build --strict
```

Every PR bumps `orfdata/version.py` and adds a `CHANGELOG.md` entry.
Maintainers publish with `./deploy.sh` from a clean, synchronized `main`.

## License

Apache 2.0. Downloadable datasets keep their sources' own terms and
attribution, recorded in each release manifest.
