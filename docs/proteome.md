# Build a reference proteome

A `Proteome` holds protein sequences for one species and genome. It is built
in explicit layers: annotated proteins first, then optional extra ORFs and
viral proteins. Each method returns a new proteome, and `proteome.layers`
records every addition with its source, selection and counts.

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

genome = EnsemblRelease(115, species="human")  # Installed with pyensembl install.
annotated = Proteome.from_genome(genome)
proteome = annotated.with_extra_orfs(evidence="ms")
```

The examples on this page continue that Python session.

## Annotated proteins

`Proteome.from_genome(genome)` takes every transcript with an annotated protein
sequence. For Ensembl this matches the annotation's protein FASTA: Ensembl 115
has 245,535 human proteins with 160,811 distinct sequences. They include
nonsense-mediated-decay and immunoglobulin/T-cell receptor segment translations,
and 11,540 proteins on patch, haplotype and unplaced contigs. Each entry
records its accession, gene, transcript, contig and transcript biotype, so you
can restrict them explicitly:

```python
primary = {str(n) for n in range(1, 23)} | {"X", "Y", "MT"}
chromosomes = annotated.filter(lambda e: e.contig in primary, description="primary chromosomes")
print(len(chromosomes))
```

```text
233995
```

`filter` records its description and the number of removed entries as a layer.
Mouse, dog, yeast and other PyEnsembl species work the same way; mouse Ensembl
115 has 66,668 annotated proteins. orfdata does not download or index
annotation data. An uninstalled genome raises an error naming the
`pyensembl install` command.

## Choose the evidence

`with_extra_orfs(evidence=...)` selects catalogue ORFs by evidence, binds each
one to this genome's transcript and adds those that match exactly. Counts below
are for the bundled human catalogue (dataset `2026-10-11.1`) and Ensembl 115:

| Selection | Human ORFs | Bound |
| --- | ---: | ---: |
| `evidence=None` or `"ribosome"`: translated in ribosome profiling | 7,264 | 7,108 |
| `evidence="ms"` (default): any supporting MS peptide | 2,186 | 2,140 |
| `evidence="hla"` or `"mhc"`: peptides presented by HLA | 1,881 | 1,843 |
| `evidence="ms_unclassified"`: MS with an unreported assay type | 537 | 525 |
| `evidence="shotgun"`: protease-digest MS | 73 | 71 |
| `"ms"` with `require_unique=True` | 1,623 | 1,587 |
| `"ms"` with `require_reviewed=True` | 44 | 43 |

All bundled ORFs come from ribosome-profiling catalogues, so `"ribosome"`
selects all of them. An MS-supported ORF can have peptides from several assay
types. `require_unique` keeps peptides that map to one protein *within the MS
source's search*, which is not uniqueness across every possible protein; in
this catalogue only the author HLA mappings report it. `require_reviewed`
keeps peptides with convincing community spectrum reviews. Because those two
properties come from different sources, no bundled peptide has both.
`require_resolved` excludes support that consists only of aggregate counts.
See [ORFs and evidence](reference.md#evidence) for the full definitions.

A translated ORF or a detected peptide does not by itself establish a stable,
functional protein. Missing evidence does not establish that an ORF is never
translated.

## Excluded and duplicate ORFs

Selected ORFs that do not bind are kept in `proteome.excluded`:

```python
from collections import Counter

print(Counter(e.reason for e in proteome.excluded))
print(dict(proteome.layers[-1]))
```

```text
Counter({'Transcript absent from this annotation': 25, 'ORF splice path mismatch': 21})
{'kind': 'extra_orf', 'source': 'orfdata-catalog:2026-10-11.1', 'evidence': 'ms', 'filters': {}, 'selected': 2186, 'added': 2140, 'excluded': 46, 'identical_to_annotated': 18}
```

The bundled models were validated against Ensembl 101; later releases retire
or change some of their transcripts. orfdata never moves an ORF to a different
transcript or sequence, even though 41 of these 46 ORFs have an identical
spliced coding path on another Ensembl 115 transcript. `identical_to_annotated` counts added ORFs whose
protein now matches an annotated protein. Both entries are kept because they
have different identities and evidence.

Without `catalog=`, the catalogue covering the genome's species and assembly is
used. Only human GRCh38 has one; other species raise. Pass a catalogue from
[`CatalogRegistry`](catalogs.md) or `load_catalog(path)` explicitly. A
catalogue for a different species or assembly is rejected.

## Viral proteomes

A `ViralProteome` holds one virus's proteins, its NCBI taxon, and the NCBI
species taxa recorded as its hosts, with the source of that host range. It can
only be added to a proteome for one of those hosts:

```python
from orfdata import ViralProteome

# Template: substitute your FASTA, the virus taxon and your host-range source.
hpv16 = ViralProteome.from_fasta(
    "hpv16-proteins.fasta", name="HPV16", taxon_id=333760,
    host_taxon_ids={9606}, host_source="Host range source and date",
)
with_virus = proteome.with_viral(hpv16)
```

Adding it to a mouse or dog proteome raises an error. orfdata does not infer
hosts from a virus's name or from related viruses. Viral entries keep the
virus's taxon in their FASTA headers.

## Lower-level ORF lookup

`extra_orfs(genome, evidence=None)` returns the bound ORFs that a proteome layer
uses, plus the same exclusions. Look them up with `for_transcript` and
`for_gene`. Each item is a [`ReferenceORF`](reference.md#objects-and-queries)
with the original PyEnsembl transcript.

## FASTA headers

`write_fasta` writes 60-column sequences. The first word of each header is the
accession: an Ensembl protein ID with its version, a source ORF ID or a viral
protein ID. Fields that apply follow as `key=value`:

```text
>ENSP00000269305.4 kind=annotated taxon=9606 source=ensembl:115 gene=ENSG00000141510 gene_name=TP53 transcript=ENST00000269305.9 contig=17 biotype=protein_coding
>c7riboseqorf56 kind=extra_orf taxon=9606 source=orfdata-catalog:2026-10-11.1 gene=ENSG00000169902 gene_name=TPST1 transcript=ENST00000304842.6 contig=7 biotype=uORF evidence=ribosome,mhc
```

For extra ORFs, the accession is the source ORF ID, `transcript` is its host
transcript in this genome, and `gene` and `gene_name` are that transcript's gene
in this genome. `biotype` is the source ORF class and `evidence` lists the
supporting evidence types; `mhc` includes HLA. See
[identifiers](reference.md#identifiers).

## Command line

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
orfdata proteome --species human --release 115 --extra-orf-evidence all
orfdata proteome --species mouse --release 115 --extra-orf-evidence none
```

`--extra-orf-evidence` accepts the evidence names above, `all` for every
catalogue ORF, or `none` for annotated proteins only. `--require-unique`,
`--require-reviewed`, `--require-resolved` and `--catalog FILE` match the Python
options. The FASTA goes to `--output` or standard output; a JSON summary of the
layers and excluded ORFs goes to standard error.

## Not yet included

Mouse and dog noncanonical ORFs, and translated endogenous retrovirus (ERV)
ORFs without annotated transcripts, need new source data; see the
[data roadmap](https://github.com/openvax/orfdata/issues/3). Translation
products from circular RNAs, stop-codon readthrough or fusions are outside this
reference.
