# Use ORFdata with varcode, isovar, topiary and vaxrank

Those tools read a [PyEnsembl](https://github.com/openvax/pyensembl) genome: a
GTF of transcript models, plus transcript and protein FASTAs whose identifiers
match the GTF. ORFdata writes exactly that for the two things a reference
annotation leaves out — noncanonical ORFs and viral references — so a variant
in either gets a coding consequence from the same tools that annotate reference
transcripts.

## Three layers, three jobs

ORFdata keeps the molecular layers apart and names them the same way in every
command and API:

| Layer | Molecule | What it answers | Who reads it |
| --- | --- | --- | --- |
| `genome` | DNA | Where a read or a variant call sits | Aligners; the contigs a VCF names |
| `transcriptome` | RNA | What a variant does to a reading frame | varcode, isovar |
| `proteome` | protein | Which peptide a patient could present | topiary, vaxrank, MS search |

A DNA reference such as GDC's `GRCh38.d1.vd1` is the genome layer: human
chromosomes plus the viral sequences as extra contigs. It carries no transcript
models and no proteins, which is what the other two layers supply.

## Viral transcripts in a tumour sample

Install the collection, or one cancer subset of it, and open it as a genome:

```sh
orfdata download viral --species human --subset cancer-causing
orfdata viral summary --subset cancer-causing
```

```python
from varcode import Variant

from orfdata import index_viral_references, viral_genome

index = index_viral_references(
    subsets=["cancer-causing"], species="human", download_if_missing=True
)
genome = viral_genome(index)  # A PyEnsembl genome over the viral contigs.

variant = Variant(contig="AY446894.2", start=1500, ref="C", alt="T", genome=genome)
print(variant.effects().top_priority_effect().short_description)
```

`viral_genome()` indexes the set's `annotation.gtf`, `transcripts.fasta` and
`proteins.fasta` into a PyEnsembl database kept beside the cache, never inside
the checksummed reference directory. The same object goes to isovar with the
tumour RNA alignment, and to topiary or vaxrank with an MHC predictor, wherever
those take a genome.

Four things are worth knowing before the first run:

- **Contig names are the source `accession.version`**, such as `AY446894.2`,
  because that is the record the sequences and the annotation come from.
  `genomes.tsv` keeps the inventory's `source_genbank` and `source_refseq`
  spellings beside it, so you can map the names your alignment used.
- **A non-standard genetic code is carried**, as `transl_table` on the
  transcript and CDS rows, so PyEnsembl translates with the code the source
  record declares.
- **A partial or frameshifted source CDS gets no start or stop codon rows**, so
  effect predictors will not translate it. `annotation_complete` says which
  transcripts are whole and `partial` says which end the source marks as
  running past the record: 78 of 1,703, mostly HPV E4 and E5 ORFs that are
  really expressed from spliced transcripts, plus HIV-1's `pol`, which is a
  Gag-Pol frameshift fusion. HPV16's E6, E7, E1, E2, L1 and L2 are complete,
  and HBV's polymerase is joined back together across the circular origin, so
  variants on either side of it get a coding consequence.
- **A coding transcript here is the CDS span**, not a UTR-containing mRNA. The
  source records annotate 1,678 CDS features but only 39 `mRNA` features, so
  the transcriptome layer is built from the CDS features and says so rather
  than inventing transcript bounds.
- **Three source records annotate no CDS at all** — HPV2, HPV8 and HPV51 — so
  their coding features come from PaVE instead, on the identical sequence, and
  every transcript records which source annotated it. Use `--no-annotation-fill`
  to keep strictly to the GenBank submissions; `orfdata viral summary` then
  names the genomes that are left without transcripts. For the opposite
  trade — one source's uniform, spliced gene models everywhere it publishes
  the same molecule, in place of each submission's own naming — use
  `--alternate-annotations`, which leaves the genome FASTA byte-identical and
  adds 235 coding transcripts across the GDC collection. `orfdata viral
  completeness` shows what either choice leaves on the table.
- A transcript's `initiation` says how the source starts it (`atg`,
  `alternative`, `declared`, `non_canonical`, `partial`, `unknown`) and
  `recoding` names any frameshift the source claims. Effect predictors handle
  a non-canonical start correctly, because the protein a source publishes is
  the one in `proteins.fasta`; they do not handle a variant at the base a
  ribosomal-frameshift join reads twice, which is two positions in HTLV-1.

## Noncanonical ORFs

Build a proteome as usual and ask for the matching annotation:

```sh
orfdata proteome --species human --release 115 --nc-orf-sources ms \
  --output proteins.fasta --annotation ncorfs/
```

That writes `ncorfs/annotation.gtf`, `ncorfs/transcripts.fasta` and
`ncorfs/proteins.fasta`: one coding transcript per ncORF, placed on its host
gene, with the spliced coding sequence and the translation. In Python the same
is `proteome.write_annotation("ncorfs")`.

```python
from pyensembl import Genome

genome = Genome(
    reference_name="GRCh38",
    annotation_name="orfdata-ncorfs",
    annotation_version="1",
    gtf_path_or_url="ncorfs/annotation.gtf",
    transcript_fasta_paths_or_urls=["ncorfs/transcripts.fasta"],
    protein_fasta_paths_or_urls=["ncorfs/proteins.fasta"],
    cache_directory_path="ncorfs/index",
)
genome.index()
```

A variant inside an ncORF now has a protein consequence — a stop gained in a
lncRNA ORF reads as `p.K2*` rather than disappearing into "non-coding
transcript".

To keep the reference transcripts as well, concatenate the files instead of
replacing them. PyEnsembl takes several FASTAs, and the ncORF GTF writes no
`gene` rows precisely so it can be appended to an annotation that has them:

```sh
cat Homo_sapiens.GRCh38.115.gtf ncorfs/annotation.gtf > combined.gtf
```

```python
genome = Genome(
    reference_name="GRCh38", annotation_name="ensembl-plus-ncorfs", annotation_version="115",
    gtf_path_or_url="combined.gtf",
    transcript_fasta_paths_or_urls=["Homo_sapiens.GRCh38.cdna.all.fa", "ncorfs/transcripts.fasta"],
    protein_fasta_paths_or_urls=["Homo_sapiens.GRCh38.pep.all.fa", "ncorfs/proteins.fasta"],
    cache_directory_path="combined/index",
)
```

For the models alone, without building a proteome, the query command writes the
same GTF rows for any selection:

```sh
orfdata --evidence ms --format gtf > ncorfs.gtf
orfdata --gene LARP4B --format gtf
```

Two deliberate choices in that file:

- ncORF transcripts are marked `transcript_biotype "protein_coding"`, because
  each has a complete CDS with start and stop codons and effect predictors only
  translate transcripts marked that way. The source's own category stays in
  `orf_biotype`, and the host transcript in `host_transcript_id`.
- The transcript is the ORF's coding path, so its identifier is the catalogue's
  `orf_id`. The protein is `<orf_id>.p`. Nothing is renamed to look like an
  Ensembl accession.

## What this does not claim

An ORF model here is a catalogue ORF placed on an assembly, with the evidence
the source reported; a predicted coding consequence is not evidence that a
stable protein exists or that a peptide is presented. Viral collection and
context membership record where a reference was selected, not that a virus
drives a particular tumour. The layers stay separate so each of those questions
can be asked of the right one.
