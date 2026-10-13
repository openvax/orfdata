# orfdata

orfdata builds the protein reference for a species: the proteins annotated in
a [PyEnsembl](https://github.com/openvax/pyensembl) genome, plus open reading
frames (ORFs) that ribosome-profiling studies found translated outside those
annotations. Each added ORF stays bound to its PyEnsembl transcript, with its
source, category and evidence. Viral proteomes can be added explicitly, but only for a species that
is a recorded host of the virus.

## Install

Use Python 3.9 or later. Install orfdata and one PyEnsembl annotation:

```sh
python -m pip install 'orfdata[pyensembl]'
pyensembl install --release 115 --species human
```

The second command downloads and indexes Ensembl 115 for human GRCh38. orfdata
never downloads or indexes annotation data itself; it reports what is missing.

## Build a human proteome

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

`Proteome.from_genome` takes every protein in the Ensembl 115 annotation.
`with_extra_orfs(evidence="ms")` adds the 2,186 bundled human ORFs that have
supporting mass-spectrometry (MS) peptides. The 2,140 that bind exactly to
Ensembl 115 transcripts are added. The other 46 are listed in
`proteome.excluded` with the reason: their transcript is absent from this
release, or its splice path changed. Building this human proteome takes about a
minute.

The same proteome from the command line:

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
```

## Choose the evidence

| `evidence=` | Selects ORFs with | Human ORFs | Bound to Ensembl 115 |
| --- | --- | ---: | ---: |
| `"ribosome"` | Ribosome profiling (all bundled ORFs) | 7,264 | 7,108 |
| `"ms"` (default) | Any supporting MS peptide | 2,186 | 2,140 |
| `"hla"` | Peptides presented by HLA molecules | 1,881 | 1,843 |
| `"shotgun"` | Peptides from protease-digest MS | 73 | 71 |

Most MS support in the bundled human catalogue comes from HLA
immunopeptidomics. A presented peptide shows that the ORF was translated; it
does not establish a stable protein. See
[building a reference proteome](proteome.md#choose-the-evidence) for
stricter filters and what each level means.

## Look up the ORFs on a transcript

Continue in the same Python session:

```python
from orfdata import extra_orfs

orfs = extra_orfs(genome)
for orf in orfs.for_transcript(genome.transcript_by_id("ENST00000304842")):  # TPST1
    print(orf.model.orf_id, orf.model.orf_type, orf.protein_sequence)
    print(orf.evidence.has_evidence("ms"))
```

```text
c7riboseqorf56 uORF MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV
True
```

`extra_orfs(genome)` binds the catalogue to every matching transcript of this
genome; without `evidence=`, it includes every catalogue ORF. Each result is a
`ReferenceORF`: the source ORF, its evidence and the original PyEnsembl
transcript, after checking the spliced coding path and protein sequence. One RNA
can have several ORFs, even when Ensembl annotates it as noncoding; the
transcript's biotype and annotated protein are unchanged.

## Species

| Species | Annotated proteins | Noncanonical ORFs |
| --- | --- | --- |
| Human | From PyEnsembl | 7,264 GRCh38 ORFs, bundled |
| Mouse | From PyEnsembl | None released; [OpenProt adapter](https://github.com/openvax/orfdata/issues/18) in progress |
| Dog | From PyEnsembl | None released; no qualified source yet |
| Budding yeast | From PyEnsembl | None released |

`with_extra_orfs()` raises for a species without a catalogue rather than
returning an annotated-only proteome. Human ORF evidence is never transferred
to another species. The same API also serves
[RNA, ribosome and protein measurements](measurements.md) for human, mouse, dog
and yeast; see [released data and gaps](coverage.md).

## Next steps

- [Build a reference proteome](proteome.md): evidence filters, excluded ORFs,
  restricting contigs, viral proteomes and FASTA headers.
- [ORFs and evidence](reference.md): the ORF classes, what each kind of
  evidence supports, and the per-transcript API used by Varcode and Isovar.
- [Species, references and genetic codes](species.md): reference identity and
  mitochondrial translation tables.
- [Measurements](measurements.md): quantitative datasets with native units.
