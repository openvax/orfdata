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
{'annotated': 245535, 'extra_orf': 2181} 5
```

`Proteome.from_genome` takes every protein in the Ensembl 115 annotation.
`with_extra_orfs(evidence="ms")` adds the bundled human ORFs that have
supporting mass-spectrometry (MS) peptides: 2,181 of the 2,186 are encoded by
an Ensembl 115 transcript, which must contain the ORF's exact spliced coding
path. The other 5 are listed in `proteome.excluded` with the reason. Each entry
lists every transcript encoding it. Building this human proteome takes under a
minute, and `proteome.save(path)` keeps it as a checksummed snapshot.

The same proteome from the command line:

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
```

## Choose the evidence

| `evidence=` | Selects ORFs with | Human ORFs | Encoded in Ensembl 115 |
| --- | --- | ---: | ---: |
| `"ribosome"` | Ribosome profiling (all bundled ORFs) | 7,264 | 7,239 |
| `"ms"` (default) | Any supporting MS peptide | 2,186 | 2,181 |
| `"hla"` | Peptides presented by HLA molecules | 1,881 | 1,879 |
| `"shotgun"` | Peptides from protease-digest MS | 73 | 72 |

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

`extra_orfs(genome)` places every catalogue ORF on every transcript of this
genome that carries its exact coding path; without `evidence=`, it includes
every catalogue ORF. Each result is a `TranscriptORF`: the source ORF, its
evidence and the original PyEnsembl transcript, with `start` and `end` slicing
the transcript sequence. One RNA can have several ORFs, even when Ensembl
annotates it as noncoding; the transcript's biotype and annotated protein are
unchanged. Add `include_annotated=True` to list the annotated CDSs on that
transcript alongside them.

## Species

| Species | Annotated proteins | Noncanonical ORFs |
| --- | --- | --- |
| Human | From PyEnsembl | 7,264 GRCh38 ORFs, bundled; lifted to GRCh37 and NCBI36 |
| Mouse | From PyEnsembl | 2,976 GRCm39 ORFs from OpenProt 2.2, [built locally](catalogs.md#mouse-orfs-from-openprot); lifted to GRCm38 and NCBIM37 |
| Dog | From PyEnsembl | None released; no qualified source yet |
| Budding yeast | From PyEnsembl | None released |

Other assemblies use a catalogue lifted to them once `orfdata chains download`
has fetched the chain files; see [every assembly of a species](catalogs.md#every-assembly-of-a-species).
`with_extra_orfs()` raises for a species without a catalogue, or whose
catalogue or chains aren't installed yet, rather than returning an
annotated-only proteome. Human ORF evidence is never transferred
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
