[![Tests](https://github.com/openvax/orfdata/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/orfdata/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/orfdata.svg)](https://pypi.org/project/orfdata/)

# orfdata

orfdata builds protein references from transcript models in
a [PyEnsembl](https://github.com/openvax/pyensembl) annotation, plus open reading
frames (ORFs) that ribosome-profiling studies found translated outside those
annotations. These include upstream and downstream ORFs, ORFs in other reading
frames of a coding sequence, and ORFs on RNAs annotated as noncoding. Each added
ORF stays bound to its PyEnsembl transcript, with its source, category and
evidence.

## Installation

Use Python 3.11 or later. Install orfdata and one PyEnsembl annotation:

```sh
python -m pip install 'orfdata[pyensembl,build-data]'
pyensembl install --release 115 --species human
# For the ncORF evidence presets and "all":
orfdata setup --species human --catalogs all
```

List every bundled, downloadable or buildable resource and what is installed:

```sh
orfdata download list
orfdata download status
```

Viral download support and its Biopython/DataCache dependencies are included in
the regular installation. `orfdata setup --species human` also installs
`gdc-viruses`, the default 200-reference GDC viral screen for human
transcriptome analysis. To acquire only those source references:

```sh
orfdata download viral --species human
```

## Example

```python
from pyensembl import EnsemblRelease
from orfdata import ReferenceTranscriptome

genome = EnsemblRelease(115, species="human")
transcriptome = ReferenceTranscriptome.from_genome(genome)
reference = transcriptome.proteome()                    # Reference CDS proteins only; all isoforms.
one_per_gene = transcriptome.proteome(isoforms="longest")
proteome = transcriptome.proteome(nc_orf_sources="ms")          # Explicitly add MS-supported ncORFs.
proteome.write_fasta("human-e115.fasta")
print(proteome.summary()["protein_sources"], len(proteome.excluded))
```

```text
{'reference_cds': 245535, 'nc_orf': 2181} 5
```

`proteome()` defaults to reference CDS proteins, with no ncORFs or viruses.
`isoforms="longest"` selects one protein per gene by length, with accession
and transcript ID breaking ties. This does not assign Ensembl canonical or
MANE status. `nc_orf_sources="ms"` takes every protein in the Ensembl 115 annotation and adds the 2,181
bundled human ORFs that have supporting mass-spectrometry peptides and are
encoded by an Ensembl 115 transcript. Each entry lists every transcript encoding
it; the other 5 selected ORFs are listed in `proteome.excluded` with the
reason. `nc_orf_sources="ribosome"` includes Ribo-seq-supported ORFs from
all compatible registered datasets. Use a rule such as
`RiboSeq(min_studies=2) & MS()` to set Ribo-seq and MS stringency separately.
A dataset name such as `"bundled"` selects every ORF in that dataset;
a list such as `["bundled", "ncorf-mammals-human"]` combines datasets.
`nc_orf_sources="all"` includes every compatible registered dataset, including
optional ones. Missing datasets raise with installation instructions.
`proteome.save(path)` writes a checksummed snapshot that `Proteome.load` reads
offline. From the command line:

```sh
orfdata proteome --species human --release 115 --nc-orf-sources ms --output human-e115.fasta
```

A translated ORF or a detected peptide does not by itself establish a stable,
functional protein. Most MS support in the bundled human catalogue comes from
immunopeptidomics: peptides eluted from HLA molecules.

## Species

| Species | Noncanonical ORFs | Sources |
| --- | ---: | --- |
| Human | 7,264 bundled; 57,607 with `nc_orf_sources="all"` | GENCODE; optionally Chang et al. 2025 and Yang et al. 2024 Ribo-seq |
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
`with_nc_orfs("all")` select catalogues for a species and raise when none is
available or built.
Human evidence is never transferred to another species. Viral proteins can
be added only to a proteome for a recorded host. To build one dictionary and
FASTA with the reference CDS proteins, all binding catalogue ORFs and the
bundled viruses associated with human:

```python
combined = transcriptome.proteome(nc_orf_sources="all", viral="all")
proteins_by_id = combined.by_accession()  # Reference CDS, ncORF and viral entries.
combined.write_fasta("human-all-proteins.fasta")
combined.save("human-all-proteins.proteome.json.gz")  # Complete provenance.

nc_only = transcriptome.proteome(reference_cds=False, nc_orf_sources="all")
nc_only.write_fasta("human-reference-ncorfs.fasta")
```

For sample-specific altered sequences, `DerivedProtein` records the supplied
amino acids, parent protein accession, sample, change IDs and source before
`combined.with_derived_proteins(...)` adds them. These sequences are not
inferred from variants or included in the reference proteome by default.

```sh
orfdata proteome --species human --release 115 --nc-orf-sources all --viral all --output human-all-proteins.fasta
```

`viral="all"` selects one protein set per bundled virus and keeps its exact
reference ID in FASTA headers. To select one virus, use `viral="HPV16"` or
`with_viral("HPV16")`.
Viral proteins are included only when requested; a recorded host association
does not say which virus is present in a sample.
orfdata also provides [RNA, ribosome and protein measurements](https://openvax.github.io/orfdata/measurements/)
for human, mouse, dog and yeast.

`viral_references(host="human", group="HPV", virus_type="HPV16")` finds
source-qualified `ViralReference` objects. Choose a stable ID with
`viral_reference("NC_001526.4")`, then use `.dna` and `.orfs` to load its
packaged bases and coding paths. The functional loaders remain available.
When annotating an existing alignment, `viral_orfs(reference_dna,
contig="HPV16")` checks its actual reference sequence before applying the same
models.
The curated [coding reference](https://openvax.github.io/orfdata/viral/)
contains the nine RefSeq NC_001526.4 HPV16 CDSs, with native and GDC
placements; the host association is distinct from the virus's own taxon. The
separate `gdc-viruses` source bundle contains all 200 GDC source genomes and
their annotated CDS/protein FASTAs for human analysis. It does not assert a
human source host for every virus; source-host metadata remains separate.
Promoting further source records into the curated
`ViralReference` API remains in the [viral roadmap](https://github.com/openvax/orfdata/issues/91).

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
| [Dog peptides](https://openvax.github.io/orfdata/dog-peptides/) | Dog | DLA class I ligands from tumours (Kaabinejadian et al.); Canine PeptideAtlas peptides and protein calls, built locally |

Gene RNA and canonical protein measurements do not establish noncanonical ORF
abundance. See [released data and gaps](https://openvax.github.io/orfdata/coverage/).

## Documentation

The [documentation](https://openvax.github.io/orfdata/) covers:

| Task | Guide |
| --- | --- |
| Discover, inspect and download data | [Discover and download data](https://openvax.github.io/orfdata/downloads/) |
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
