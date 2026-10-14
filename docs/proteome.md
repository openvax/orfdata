# Build a reference proteome

A `Proteome` holds protein sequences for one species and genome. It is built
in explicit layers: annotated proteins first, then optional extra ORFs and
viral proteins. Each method returns a new proteome. `proteome.layers` records
every addition with its role (`host`, `viral` or `filter`), source, selection
and counts, and each entry lists every transcript that encodes it.

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

`with_extra_orfs(evidence=...)` selects catalogue ORFs by evidence and adds each
one that a transcript of this genome encodes: the transcript contains the ORF's
exact spliced coding path and translates to its protein. Counts below are for
the bundled human catalogue (dataset `2026-10-11.1`) and Ensembl 115:

| Selection | Human ORFs | Encoded in Ensembl 115 |
| --- | ---: | ---: |
| `evidence=None` or `"ribosome"`: translated in ribosome profiling | 7,264 | 7,239 |
| `evidence="ms"` (default): any supporting MS peptide | 2,186 | 2,181 |
| `evidence="hla"` or `"mhc"`: peptides presented by HLA | 1,881 | 1,879 |
| `evidence="ms_unclassified"`: MS with an unreported assay type | 537 | 535 |
| `evidence="shotgun"`: protease-digest MS | 73 | 72 |
| `"ms"` with `require_unique=True` | 1,623 | 1,621 |
| `"ms"` with `require_reviewed=True` | 44 | 44 |

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

## Every transcript, and excluded ORFs

Each extra ORF is one entry, whatever the number of transcripts encoding it.
`entry.occurrences` lists them, the transcript the source names first. Each
`Occurrence` has the transcript and gene, `start` and `end` slicing the
transcript sequence, the genomic `coding_blocks` and `on_source_transcript`.
Annotated entries have one occurrence, on their own transcript.

```python
from collections import Counter

ms_orfs = [e for e in proteome if e.kind == "extra_orf"]
print(Counter(min(len(e.occurrences), 5) for e in ms_orfs))
print(Counter(e.reason for e in proteome.excluded))
```

```text
Counter({5: 1253, 1: 345, 2: 254, 3: 184, 4: 145})
Counter({'ORF splice path mismatch': 5})
```

More than half of the MS-supported ORFs, 1,253 of 2,181, are encoded by five
or more Ensembl 115 transcripts. The bundled models were validated against Ensembl 101, and later
releases retire or change some of their transcripts. When an ORF's source
transcript is gone or changed but another transcript still encodes it, the ORF
is added and `entry.source_binding` gives the reason: 41 of the MS-supported
ORFs, counted as `not_on_source_transcript` in the layer. The 5 in
`proteome.excluded` are encoded by no Ensembl 115 transcript.
`identical_to_annotated` in the layer counts added ORFs whose protein now
matches an annotated protein; both entries are kept because they have
different identities and evidence. `proteome.by_sequence()` groups entries with
identical sequences, for example before an exact peptide search.

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
virus's taxon, and their layer has the role `viral`; `proteome.host()` returns
the host's own proteins without them, for example as a self reference.

## Save and reload

```python
identity = proteome.save("human-e115.proteome.json.gz")
same = Proteome.load("human-e115.proteome.json.gz")  # Offline: no genome or catalogue needed.
assert same.identity == identity
```

A snapshot holds every entry with its occurrences and evidence, the layers with
their selections and catalogue identities, and the excluded ORFs. `identity` is
the SHA-256 of that content, so it changes when any sequence, occurrence,
selection, catalogue or exclusion changes, but not with file compression or
formatting. Loading rejects a snapshot whose content does not match its
identity, host proteins from another species, viral proteins without a
matching host layer and layers declaring another reference. The human Ensembl
115 proteome with MS-supported ORFs is a 57 MB snapshot.

## Lower-level ORF lookup

`extra_orfs(genome, evidence=None)` returns every catalogue ORF on every
transcript carrying it, as [`TranscriptORF`](reference.md#objects-and-queries)
objects, with `source_failures` and `excluded`. Its default selects no evidence
level: use it for annotation, such as variant effects on uORFs, where the
proteome's MS default would be the wrong scope. Look ORFs up with
`for_transcript` and `for_gene`; see [identifiers](reference.md#identifiers)
for how versions and references are checked.
`for_transcript(transcript, include_annotated=True)` adds the annotated CDSs
on that transcript, and `overlapping(...)` finds ORFs at a genomic interval.

## FASTA headers

`write_fasta` writes 60-column sequences. The first word of each header is the
accession: an Ensembl protein ID with its version, a source ORF ID or a viral
protein ID. Fields that apply follow as `key=value`:

```text
>ENSP00000269305.4 kind=annotated taxon=9606 source=ensembl:115 gene=ENSG00000141510 gene_name=TP53 transcript=ENST00000269305.9 contig=17 transcript_biotype=protein_coding
>c7riboseqorf56 kind=extra_orf taxon=9606 source=orfdata-catalog:2026-10-11.1 gene=ENSG00000169902 gene_name=TPST1 transcript=ENST00000304842.6 contig=7 transcript_biotype=protein_coding orf_type=uORF evidence=ribosome,mhc transcripts=18
>c10riboseqorf115 kind=extra_orf taxon=9606 source=orfdata-catalog:2026-10-11.1 gene=ENSG00000099204 gene_name=ABLIM1 transcript=ENST00000369253.6 contig=10 transcript_biotype=protein_coding orf_type=intORF evidence=ribosome,mhc source_transcript=ENST00000277895 transcripts=10
```

For extra ORFs, the accession is the source ORF ID, `transcript` is its host
transcript in this genome, and `gene`, `gene_name` and `transcript_biotype`
describe that transcript in this genome. `orf_type` is the source's ORF
category, kept verbatim; see [ORF categories](reference.md#orf-categories).
`evidence` lists the supporting evidence types; `mhc` includes HLA.
`transcripts` counts the encoding transcripts when there are several, and
`source_transcript` names the source's transcript when the header's transcript
is a different one, as for ABLIM1 above, whose source transcript is gone. See
also [identifiers](reference.md#identifiers). The snapshot keeps every
occurrence; the FASTA header shows the first.

## Command line

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
orfdata proteome --species human --release 115 --extra-orf-evidence all
orfdata proteome --species mouse --release 115 --extra-orf-evidence none
```

`--extra-orf-evidence` accepts the evidence names above, `all` for every
catalogue ORF, or `none` for annotated proteins only. `--require-unique`,
`--require-reviewed`, `--require-resolved` and `--catalog FILE` match the Python
options. The FASTA goes to `--output` or standard output; `--snapshot FILE`
also saves a snapshot. A JSON summary of the layers, identity and excluded ORFs
goes to standard error.

## Not yet included

Mouse and dog noncanonical ORFs, and translated endogenous retrovirus (ERV)
ORFs without annotated transcripts, need new source data; see the
[data roadmap](https://github.com/openvax/orfdata/issues/3). Products from
circular RNAs, stop-codon readthrough or non-AUG starts can be described with
`TranscriptORF.at` and `SourceRNA` (see the [examples](examples.md)), but no
bundled catalogue supplies them yet, so proteomes do not include them.
