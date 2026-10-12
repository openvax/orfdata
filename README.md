# orfdata

Reference open reading frames and experimental translation evidence, as an
auxiliary layer on [PyEnsembl](https://github.com/openvax/pyensembl) transcripts.
This is a development draft of a new package; version 0.1.0 is not published.

```python
from pyensembl import EnsemblRelease
from orfdata import load_catalog

reference = EnsemblRelease(95, species="human")  # install/cache explicitly
transcript = reference.transcript_by_id("ENST00000304842")
orf, = load_catalog().for_transcript(transcript)

orf.transcript             # original PyEnsembl Transcript
orf.model                  # immutable ORF identity, class, coordinates and sequence
orf.evidence               # ribosome and peptide-level observations
orf.protein_sequence       # validated against the selected transcript sequence
orf.evidence.has_evidence("ms")
```

The reference transcript retains its canonical protein and annotated biotype.
Several ORFs may bind one transcript, including transcripts annotated as lncRNA
or processed transcripts. Binding validates the entire spliced coding path and
protein, including the stop, against the selected annotation. A matching gene
name or unversioned transcript ID alone cannot transfer evidence.

Varcode describes modifications to this same base object:

```python
from varcode import Variant, annotate_orf_variants

variant = Variant("7", 66205486, "C", "G", genome=reference)
annotation = annotate_orf_variants([variant], orf)
assert annotation.reference is orf
print(annotation.effect)  # P2A in the TPST1 upstream frame
result = annotation.to_dict()  # reference model/evidence + realized altered products
```

Isovar's `annotate_sv_orf_evidence` connects the same API to reconstructed RNA
paths and retains each candidate's full-interval RNA support. The normal TPST1
peptide KIILSKISV is absent from the reconstructed Sid TPST1–CRCP alternatives;
its reference MS evidence does not establish translation of those products.

## Data and evidence

`load_catalog()` includes 7,264 human GRCh38 reference ncORFs: upstream,
overlapping upstream, internal, downstream, overlapping downstream, lncRNA and
processed-transcript ORFs. There are 2,186 models with reported supporting MS
observations under the source-specific filters. This is **not** a count of
validated functional proteins. The library distinguishes ribosome evidence,
protease-digest MS, MHC/HLA peptidomics, ambiguous mappings and spectrum reviews.

Data are drawn from the [GENCODE Phase I reference](https://www.nature.com/articles/s41587-022-01369-0),
[TransCODE author exports](https://www.nature.com/articles/s41586-026-10459-x),
and the community microprotein spectrum benchmark. The packaged manifests retain
source accessions, checksums, licenses, processing rules and dataset version.
Ensembl 101 coding **and noncoding** RNA FASTAs validate all 7,264 full ORF paths.
Contexts can be partial at transcript boundaries; biological cap sites and
initiation efficiency remain unmeasured. Dataset version is `2026-10-11.1`.

`load_uorf_catalog()` retains the released Isovar upstream-only snapshot
(3,771 models, version `2026-10-08.2`). These compatibility imports are kept
separate from the generic API. Core catalogue access is offline and uses the
Python standard library. PyEnsembl binding accepts already-loaded transcripts;
install the `pyensembl` extra when needed.

## Other species

The API is species-independent; the bundled observations are human-only. Load a
versioned JSON catalogue for another species with `load_catalog(path)`. Each
model explicitly identifies the NCBI taxon, scientific name, assembly, source
and translation table. Species and assembly must agree with the supplied
PyEnsembl transcript. Human observations are never assigned to a mouse ortholog.

```python
from orfdata import ORFModel, ORFRecord, ORFEvidence, ORFCatalog

# Obtain real model coordinates and evidence from the selected source release.
model = ORFModel(
    orf_id=source_orf_id, transcript_id=transcript_id, gene_id=gene_id,
    gene_name=gene_name, contig=contig, strand=strand,
    coding_blocks=coding_blocks, protein_sequence=protein_sequence,
    orf_type=source_orf_class, genome_build="GRCm39", taxon_id=10090,
    species="mus_musculus", translation_table=1, source="source-release-id",
)
record = ORFRecord.from_model(model, ORFEvidence(ribosome_studies=(), peptides=()))
catalog = ORFCatalog([record], metadata)  # species/assembly/source/version metadata
# JSON serialization: json.dump(catalog.to_dict(), stream)
```

[OpenProt](https://openprot.org/p/help) is a credible source for additional species:
it provides organism/assembly-specific models, Ensembl/RefSeq transcript mappings,
and separate MS/translation evidence. An OpenProt source adapter must pin its
release, preserve mapping ambiguity and licensing, and validate its coordinate
conventions before importing observations. No OpenProt importer or mouse evidence
bundle is claimed here. Tests exercise two synthetic mouse ORFs on one real
PyEnsembl lncRNA transcript; fixture sequences are explicitly not biological data.

The initial sequence binder supports complete stopped ORFs using NCBI table 1
and its ATG/CTG/TTG initiation convention (all bundled models start at ATG).
Other near-cognate initiation conventions require explicit future support.
Other genetic codes and incomplete ORFs fail explicitly. ORF discovery, inferred
initiation efficiency and downstream CDS regulation are separate future work.

## Rebuild and verify

```sh
python -m pip install -e '.[dev]'
python -m orfdata.build --source-dir /path/to/sources --prepare-sources \
    --all-orfs --output orfdata/data/orf-evidence
./lint.sh
./test.sh
```

Source preparation downloads pinned public inputs; extracting the author export
requires Rscript. Ordinary catalogue queries require neither R nor spreadsheet
libraries. See [design](docs/design.md) and the bundled notices for provenance.
