# Build a reference proteome

A `ReferenceTranscriptome` names the transcript models in one installed
PyEnsembl annotation. A `Proteome` selects protein products from those models:
translations of the reference CDSs, catalogue ncORFs bound to the transcripts,
and optionally viral proteins associated with the host. All three sources have
annotations; `reference_cds` identifies the proteins induced by the reference
annotation's CDS features. Each method returns a new proteome. `proteome.layers` records
every addition with its role (`host`, `viral`, `component`, `external`,
`derived` or `filter`), source, selection
and counts, and each protein lists every transcript that encodes it.

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome, ReferenceTranscriptome

genome = EnsemblRelease(115, species="human")  # Installed with pyensembl install.
transcriptome = ReferenceTranscriptome.from_genome(genome)
reference_cds = transcriptome.proteome()                  # Default: all reference CDS isoforms.
one_per_gene = transcriptome.proteome(isoforms="longest")
proteome = transcriptome.proteome(nc_orf_sources="ms")            # Add MS-supported ncORFs.
```

The examples on this page continue that Python session. `nc_orf_sources=()` and
`viral=()` are the defaults; a reference-only build does not load ncORF datasets.
The ncORF argument accepts one string or a list:

| Selection | Included ncORFs |
| --- | --- |
| `[]` or `None` | None |
| `"all"` | Every compatible registered dataset, including optional ones |
| `"bundled"` | Every ORF in the bundled human dataset |
| `["bundled", "ncorf-mammals-human"]` | The union of these datasets |
| `"ms"` | MS-supported ORFs across all compatible datasets |
| `["proteomics", "ribosome"]` | ORFs with either kind of evidence across all compatible datasets |
| `["bundled", "ms"]` | Every bundled ORF, plus MS-supported ORFs from the other datasets |

Evidence rules and loaded `ORFCatalog` objects also work as sources. A list
always takes the union; combine evidence rules with `&` when both are required.
Species and assembly are taken from the transcriptome. Required datasets and
liftover chains must already be installed. For all human datasets, run
`orfdata setup --species human --catalogs all` once. Use `"bundled"` for the
bundled data alone. The dataset loader/setup commands retain their own catalogue
selectors; `proteome()` and `Proteome.build()` have only `nc_orf_sources`.

## Reference CDS proteins

`transcriptome.proteome()` and `Proteome.from_genome(genome)` take every
transcript with a reference CDS protein sequence. For Ensembl this matches the annotation's protein FASTA: Ensembl 115
has 245,535 human proteins with 160,811 distinct sequences. They include
nonsense-mediated-decay and immunoglobulin/T-cell receptor segment translations,
and 11,540 proteins on patch, haplotype and unplaced contigs. Each protein
records its accession, gene, transcript, contig and transcript biotype, so you
can restrict them explicitly:

```python
primary = {str(n) for n in range(1, 23)} | {"X", "Y", "MT"}
chromosomes = reference_cds.filter(lambda e: e.contig in primary, description="primary chromosomes")
print(len(chromosomes))
```

```text
233995
```

`filter` records its description and the number of removed entries as a layer.
Set `isoforms="longest"` on either constructor to choose one reference protein
per gene. Selection uses amino-acid length without a terminal stop, then
accession and transcript ID in lexical order to break ties. Every selected
protein retains its original identifiers and occurrence. This is a length-based
representative, not an Ensembl canonical or MANE designation. A missing gene
ID raises an error. `layers[0]` records `isoforms`, the number `available`
before selection and the number `added`.

Isoform selection applies only to reference CDS proteins. Added ncORFs retain
all matching transcripts, and viral and derived proteins keep their own
identities. Identical sequences from different genes remain separate entries.

Mouse, dog, yeast and other PyEnsembl species work the same way; mouse Ensembl
115 has 66,668 annotated proteins. orfdata does not download or index
annotation data. An uninstalled genome raises an error naming the
`pyensembl install` command.

## RefSeq proteins (experimental)

Build from NCBI's protein FASTA and genome-download feature table, with pinned
byte checksums and an explicit `ReferenceContext`:

```python
from orfdata import ReferenceContext

context = ReferenceContext(9615, "dog", "UU_Cfam_GSD_1.0", "NCBI RefSeq", "106",
                           assembly_accession="GCF_011100685.1")
reference = Proteome.from_refseq(
    "protein.faa.gz", "feature_table.txt.gz", reference=context,
    expected_sha256={"protein": protein_sha256, "feature_table": feature_table_sha256},
)
reference.save("dog.proteome.json.gz")
```

This preserves NCBI Gene IDs, versioned RNA/protein accessions, every mapped
occurrence and all FASTA proteins. Wrong checksums, assemblies, missing proteins
and conflicting assignments fail. CDS rows without protein products are counted
in `layers[0]["cds_without_protein"]`.

The annotated layer retains both asset hashes and identifier namespaces. Its
`genomic_spans` maps each protein accession to spans in occurrence order. These
are genomic envelopes, so occurrence CDS offsets/blocks and transcript biotypes
remain unknown. No annotation indexing or downloads occur.

## Pinned Ensembl proteins

`Proteome.from_ensembl()` provides the same reference inventory and provenance
contract for Ensembl's full GTF, peptide, cDNA and ncRNA files. It reads local
plain or gzip files without building annotation indexes:

`ReferenceContext.from_genome()` from a matching PyEnsembl release is also
accepted. Provider names `Ensembl` and `ensembl` are compared without case
differences; the supplied context spelling is preserved in the snapshot.

```python
context = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115",
                           assembly_accession="GCA_000001405.29")
reference = Proteome.from_ensembl(
    "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz",
    "Homo_sapiens.GRCh38.cdna.all.fa.gz",
    "Homo_sapiens.GRCh38.ncrna.fa.gz",
    reference=context,
    expected_sha256={"protein": protein_sha256, "gtf": gtf_sha256,
                     "rna": rna_sha256, "ncrna": ncrna_sha256},
)
reference.save("human.proteome.json.gz")
```

Obtain and pin the matching files from the same [Ensembl release](https://www.ensembl.org/info/data/ftp/index.html).
The reader hashes the exact source bytes it parses. The GTF assembly accession,
versioned transcript/protein accessions and versioned gene identities must agree
with the supplied reference and sequence headers. Every peptide FASTA accession
must have a GTF assignment and its exact RNA record; a missing coding sequence or
conflicting duplicate fails the build. Identical protein sequences at different
accessions and genes remain separate entries. Identical duplicate FASTA records
are counted in the layer manifest.
Protein residues, including Ensembl's internal stop symbols and ambiguous
residues, are preserved exactly. Downstream prediction must assess whether it
supports those sequences; the reference reader does not repair them.

The annotated layer carries all four source hashes, identifier namespaces and
the complete protein inventory checksum. Its `transcript_spans` preserves each
original genomic transcript envelope in occurrence order. These are not spliced
CDS coordinates: occurrence `start`, `end` and `coding_blocks` remain unknown.
Noncoding annotation transcripts without RNA bytes are listed explicitly in
`untranslated_transcripts_without_rna`; they never become protein entries.

The same method accepts other species' matching Ensembl files and an explicit
`ReferenceContext`. Reference membership does not establish RNA abundance,
empirical translation, protein abundance or tissue restriction. The reader
does not compare translations or transfer evidence between assemblies.

## Add ncORFs with selected evidence

Pass `nc_orf_sources="ms"`, `"ribosome"`, `"all"` or an evidence rule to `proteome()`,
or add to an existing proteome with `with_nc_orfs("ms")`. Evidence presets query
all compatible registered datasets, including optional ones.
This selects catalogue ORFs by evidence and adds each
one that a transcript of this genome encodes: the transcript contains the ORF's
exact spliced coding path and translates to its protein. Counts below are for
the bundled human catalogue (dataset `2026-10-11.1`) and Ensembl 115.
The table uses the lower-level `load_catalog().query(evidence=...)` notation;
`nc_orf_sources=None` instead means no ncORFs:

| Selection | Human ORFs | Encoded in Ensembl 115 |
| --- | ---: | ---: |
| `evidence=None` or `"all"`: every catalogue ORF | 7,264 | 7,239 |
| `evidence="ribosome"` or `"riboseq"`: translated in ribosome profiling | 7,264 | 7,239 |
| `evidence="ms"`: any supporting MS peptide | 2,186 | 2,181 |
| `evidence="immunopeptidomics"`: peptides eluted from MHC (in human, HLA) molecules | 1,881 | 1,879 |
| `evidence="ms_unclassified"`: MS with an unreported assay type | 537 | 535 |
| `evidence="proteomics"`: peptides from digested cellular proteins | 73 | 72 |
| `"ms"` with `require_unique=True` | 1,623 | 1,621 |
| `"ms"` with `require_reviewed=True` | 44 | 44 |

All bundled ORFs come from ribosome-profiling catalogues, so `"ribosome"`
selects all of them. An MS-supported ORF can have peptides from several assay
types. `require_unique` keeps peptides that map to one protein *within the MS
source's search*, which is not uniqueness across every possible protein; in
this catalogue only the author HLA ligand mappings report it. `require_reviewed`
keeps peptides with convincing community spectrum reviews. Because those two
properties come from different sources, no bundled peptide has both.
`require_resolved` excludes support that consists only of aggregate counts.
See [ORFs and evidence](reference.md#evidence) for the full definitions.

A list of evidence names selects ORFs with any of them:
`nc_orf_sources=["proteomics", "immunopeptidomics"]` keeps ORFs with either kind of
MS peptide, unlike `"ms"`, which also counts unclassified MS. For ORFs that
need two kinds of evidence at once, combine rules with `&`, below.

Length is a separate consumer policy, never an import rule. RibORF catalogues
retain valid short source ORFs because they can matter for variant annotation
even when they are unsuitable for a peptide search. Use
`catalog.query(min_protein_length=5)`, `extra_orfs(...,
min_protein_length=5)`, or:

```python
search_reference = transcriptome.proteome(
    nc_orf_sources="ribosome", min_protein_length=5
)
```

The corresponding CLI option is `--min-protein-length 5`. A proteome snapshot's
ncORF layer records both the threshold and `length_filtered`, so the omission is
auditable and changes the cache identity. Reference CDS proteins are unaffected.

### Separate Ribo-seq and MS stringency

For more control, pass a rule. `RiboSeq(min_studies=n)` needs ribosome
profiling in at least `n` independent studies the source names, and
`MS(min_peptides=n, unique=..., reviewed=..., mhc=...)` needs at least `n`
distinct supporting peptide sequences. `MS(mhc=True)` counts only peptides
eluted from MHC molecules (the immunopeptidome) and `MS(mhc=False)` only
peptides from digested cellular proteins (the whole proteome). `MS()` counts
both, and the MS peptides whose source doesn't say which. `methods=` names the
kinds directly, for combinations such as `("proteomics", "ms_unclassified")`.
`&` requires both rules and `|` accepts either:

```python
from orfdata import MS, RiboSeq

reference_cds.with_nc_orfs(RiboSeq(min_studies=2))                        # replicated Ribo-seq
reference_cds.with_nc_orfs(RiboSeq(min_studies=2) & MS())                 # Ribo-seq and MS
reference_cds.with_nc_orfs(RiboSeq(min_studies=2) | MS())                 # either
reference_cds.with_nc_orfs(RiboSeq() & MS(mhc=True))                     # Ribo-seq and an MHC ligand
reference_cds.with_nc_orfs(MS(methods=("proteomics", "ms_unclassified")))  # not immunopeptidomics-only
```

| Rule | Human ORFs | Encoded in Ensembl 115 |
| --- | ---: | ---: |
| `RiboSeq(min_studies=2)` | 3,085 | 3,072 |
| `RiboSeq(min_studies=3)` | 1,675 | 1,671 |
| `MS(min_peptides=2)` | 682 | 681 |
| `MS(methods=("proteomics", "ms_unclassified"))` | 610 | 607 |
| `RiboSeq(min_studies=2) & MS()` | 1,296 | 1,293 |
| `RiboSeq(min_studies=2) & MS(min_peptides=2)` | 457 | 456 |
| `RiboSeq(min_studies=2) \| MS()` | 3,975 | 3,960 |

GENCODE's `c…riboseqorf…` IDs are the 3,085 ORFs seen in two or more studies.
Each method name above is the simplest rule: `"ms"` is `MS()`,
`"immunopeptidomics"` is `MS(mhc=True)`, `"proteomics"` is `MS(mhc=False)` and
`"ribosome"` (or `"riboseq"`) is `RiboSeq()`. The names before 0.28.0 (`"shotgun"`, and `"mhc"` or `"hla"`)
still work, with a deprecation warning, and catalogues that use them still
load. A peptide counts once
however many assay types report it, and so do peptides that differ only by I
and L (equal in mass). Nested peptides, such as MHC ligand length variants, count as
one: the count is the largest set of supporting peptides none of which lies
inside another, so a long peptide never hides two shorter, separate ones, and
a stricter rule always selects a subset of a looser one. Of the 916 ORFs with
two or more distinct peptide sequences, 234 have fewer than two that aren't
nested. Each rule in `&` is checked on its own:
`MS(mhc=False) & MS(unique=True)` needs a proteomics peptide and a
unique peptide, not one unique proteomics peptide. A source's aggregate count stands for one observation:
it satisfies `RiboSeq()` and `MS()` but never `min_studies` or `min_peptides`
above 1, `unique`, `reviewed` or `resolved=True`. The layer stores the rule,
so a snapshot records exactly how its ORFs were chosen; `EvidenceRule.from_dict`
reads it back. Equal rules are stored the same way, whatever their order or
repeats, so they give one proteome identity.

A translated ORF or a detected peptide does not by itself establish a stable,
functional protein. Missing evidence does not establish that an ORF is never
translated.

## Every transcript, and excluded ORFs

Each extra ORF produces one protein, whatever the number of transcripts encoding it.
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

Source names resolve against the genome's species and assembly. `"all"` and
evidence presets query every compatible registered dataset, using installed
liftover chains where needed. `"bundled"` selects only the bundled human dataset.
A missing dataset or chain raises with installation instructions; a species
without a dataset cannot supply ncORFs. A loaded `ORFCatalog`, or a list of
loaded catalogues, can also be passed directly as `nc_orf_sources`.

### Merging catalogues

Pass several catalogues for one species and assembly to use their ORFs together:

```python
from orfdata import load_catalog, merge_catalogs

catalogs = [load_catalog(), load_catalog("study-orfs.json.gz")]  # The bundled one and a local file.
proteome = reference_cds.with_nc_orfs(catalogs)  # Or pass merge_catalogs(catalogs).
```

Records from different catalogues are one ORF when they have the same spliced
coding path, protein, genetic code and start codon; declaring a start the code
already uses, such as ATG, changes nothing. Records within one catalogue stay
as it gave them. A record joins an ORF found on its own transcript before any
other.

- **IDs.** The merged ORF keeps the first catalogue's ID. Every source's record
  stays in `record.members`, and the others' IDs are in `equivalent_ids` as
  `source:id`. `catalog.get` and `ExtraORFs.for_orf` accept any of them.
- **ID collisions.** When two sources use one ID for different ORFs, both
  become `source:id`, and the bare ID raises rather than picking one.
  Direct `ORFCatalog(records, metadata)` construction follows the same rule,
  retaining each native record in `members` without merging its evidence.
  IDs must still be unique within a source. If a source or native ID contains
  punctuation that would collide with a generated ID, generated components
  are percent-escaped; native IDs remain available through the source records.
  Use `catalog.get("native-id", source="study")` or
  `orfs.for_orf("native-id", source="study")` for explicit source selection.
  A combined string such as `study:native-id` raises if it also identifies
  another native record; explicit source selection resolves that ambiguity.
- **Evidence.** Each source's evidence is kept, and evidence rules see the union:
  a peptide counts once per sequence. A Ribo-seq study two catalogues label
  differently counts once when their `ribosome_study_pmids` metadata give it
  the same PMID, as the bundled catalogue's does; without PMIDs it counts
  twice toward `RiboSeq(min_studies=...)`.
- **Lookups.** Queries by transcript, gene or biotype, `for_transcript` and
  `ExtraORFs.for_gene(..., source=True)` match every source's labels, and an
  ORF is on its source transcript when it is on any source's.
- **Annotation pins.** Sources pinned to different annotation releases of the
  assembly leave the merged ORF unpinned; each member keeps its own pin.
- **Rejected merges.** Catalogues must share a strain, where an unknown strain
  counts as its own value, and must not repeat a source, so two releases of
  one catalogue can't both be merged. A merged catalogue can be merged again.
- **Conflicts.** `catalog.metadata["conflicts"]` lists ORFs that share a coding
  path but differ in protein, genetic code or declared start codon.

The layer's `source` joins the catalogues' releases, such as
`orfdata-catalog:2026-10-11.1+orfdata-catalog:study-v1`. On the command line, repeat
`--nc-orf-sources FILE`.

## Viral proteomes

Build one proteome with reference CDS proteins, every catalogue ORF that binds to
the genome, and every distinct bundled viral protein set recorded for this host:

```python
combined = transcriptome.proteome(nc_orf_sources="all", viral="all")
proteins_by_id = combined.by_accession()
combined.write_fasta("human-all-proteins.fasta")
combined.save("human-all-proteins.proteome.json.gz")
```

`nc_orf_sources="all"` includes every compatible registered dataset with no
evidence filter. ORFs that do not bind to this genome are reported in
`combined.excluded`. Reference CDS and viral proteins are unaffected.

`proteins_by_id` maps each unique protein accession to a `ProteinSequence` with
its sequence, `protein_source`, taxon and source. `by_sequence()` groups entries with
identical amino-acid sequences while retaining their separate identities.
`viral="all"` chooses the native RefSeq placement when bundled references
contain the same protein set at different DNA coordinate placements, as with
HPV16. To request only one virus, use `viral="HPV16"` or an exact reference ID
such as `viral="NC_001526.4"`. `with_viral("HPV16")` accepts the same choices,
as well as a `ViralReference` returned by `viral_reference(...)`. Viral proteins
are opt-in: a recorded host association does not mean every sample contains
that virus. A virus without a recorded association to the host is rejected.
Requesting `viral="all"` for a host with no bundled viral references raises
instead of silently returning a host-only proteome.

`ProteinSet` represents named protein sequences supplied outside the bundled
references. It has no virus or host behavior of its own. When the sequences
are viral, pass the host association to `with_viral` so the same host check
applies:

```python
from orfdata import ProteinSet

# Template: substitute your FASTA, the virus taxon and your host-range source.
hpv16: ProteinSet = ProteinSet.from_fasta(
    "hpv16-proteins.fasta",
    name="HPV16",
    taxon_id=333760,
    source="Source database and release",
)
with_virus: Proteome = proteome.with_viral(
    hpv16,
    host_taxon_ids={9606},
    host_source="Host range source and date",
)
```

`viral=[first, second]` accepts several specific inputs. The tests exercise
this with two synthetic viruses; HPV16 is the only released biological viral
reference at present.

All protein kinds use the same `ProteinSequence` and appear in the same
`Proteome`. `ProteinSet` is an input, not a second output collection. Viral
proteins keep the virus's taxon and a
source label such as `viral:HPV16@NC_001526.4`; their layer records the source
accession, checksum and host-range source. `proteome.host()` returns the host's
own proteins without viral entries, for example as a self reference.

The same input type handles a selected donor, transgene or custom component.
That is a separate operation because a component's source reference and the
proteome it is selected for are different facts:

```python
from orfdata import ProteinSet, Proteome, ReferenceContext

human: ReferenceContext = ReferenceContext(
    9606, "human", "GRCh38", "Ensembl", "115"
)
mouse_reference: ReferenceContext = ReferenceContext(
    10090, "mouse", "GRCm39", "Ensembl", "115"
)
mouse: Proteome = Proteome(mouse_reference, ())
insert: ProteinSet = ProteinSet.from_fasta(
    "selected-human-insert.fasta",
    name="selected-human-insert",
    reference=human,
    source="Construct manifest 7",
)
chimera: Proteome = mouse.with_protein_set(
    insert,
    role="donor",
    target_taxon_ids={10090},
    target_source="Mouse line manifest 3",
)
```

This adds only the listed sequences, not the human proteome. The component
layer retains the human source reference and the separate mouse-target
assertion. No genomic placement, expression or sample membership is inferred.
If the `ProteinSet.reference` exactly matches the proteome reference, the
target assertion is unnecessary. Readding the same set is idempotent; changed
content under the same name is rejected.

`ProteinSet` is deliberately sequence-only. It is suitable when the source has
accessions, amino acids and provenance but no coordinate model. An ncORF source
that has a transcript, coding path, initiation context or evidence belongs in
an `ORFCatalog`/`ORFRecord` instead. Binding those ORF records to a reference
produces the same `ProteinSequence` output in the same `Proteome`; the richer
ORF facts are not flattened into a `ProteinSet`. Conversely, importing a bare
ncORF protein FASTA as a `ProteinSet` does not invent a canonical ORF,
alternative start, transcript placement or evidence.

### ncORFs with or without an installed transcriptome

If you have a PyEnsembl reference transcriptome, build a proteome containing
only the ncORFs that bind to its transcripts:

```python
nc_only = transcriptome.proteome(reference_cds=False, nc_orf_sources="all")
nc_only.write_fasta("reference-ncorfs.fasta")
```

```sh
orfdata proteome --species human --release 115 --no-reference-cds --nc-orf-sources all --output reference-ncorfs.fasta
```

This uses the installed annotation and transcript sequences to check exact
ORF paths. It does not read or include reference CDS proteins, and it does not
need a whole-genome DNA FASTA. ORFs that do not bind are reported in
`nc_only.excluded`.

If no transcriptome is installed, use the source ncORF catalogue directly;
its source transcript and assembly identifiers remain attached, but no
binding to a particular reference transcriptome is claimed:

```python
from orfdata import load_catalog

catalog = load_catalog(species="human")
nc_orfs_by_id = {record.orf_id: record for record in catalog.query(evidence="all")}
```

```sh
orfdata --species human --evidence all --format fasta > nc-orfs.fasta
```

For viral source records alone, use `viral_references(host="human")` and each
reference's `.products` and `.orfs`; this also needs no PyEnsembl genome:

```python
from orfdata import viral_reference

ref = viral_reference("NC_001526.4")
viral_by_id = {orf.orf_id: orf for orf in ref.orfs}
```

The released collection currently has HPV16 in two DNA coordinate placements.
Specifying `viral=` in `Proteome.build` only applies when building a host
proteome. With an installed genome, filter entries by `kind` to obtain only
`extra_orf` or only `viral` proteins; the result retains `write_fasta`,
`by_accession`, and `save`.

```python
nc_only = combined.filter(lambda entry: entry.kind == "extra_orf")
virus_only = combined.filter(lambda entry: entry.kind == "viral")
```

## Sample-derived proteins

Protein origin and sample-specific sequence changes are separate. A changed
reference CDS protein or ncORF retains its parent protein source, while its
derivation names the parent accession, sample and change IDs. Supply the amino
acid sequence yourself; orfdata does not infer it from a VCF or claim that the
sample expresses it:

```python
from orfdata import DerivedProtein

altered = DerivedProtein(
    accession="sample-A:protein-1",
    sequence=altered_amino_acids,
    parent_accessions=("ENSP00000269305.4",),
    sample_id="sample-A",
    change_ids=("source-variant-id",),
    source="variant-caller-and-translation-release",
)
sample_proteome = combined.with_derived_proteins(altered)
sample_proteome.write_fasta("sample-A-proteins.fasta")
```

The parent must already be in the proteome. A derived protein has `kind="derived"`
and `protein_source` inherited from its host parent (`mixed` for parents from
different sources); `entry.derivation` stores
the parent, its gene/transcript identifiers, sample, changes and optional
source asset SHA-256. Its FASTA header includes `sample`, `derived_from`,
`parent_genes`, `parent_transcripts`, `changes` and the supplied `source`.
When a parent maps to several transcripts, pass `parent_transcript_ids=(...,)`
to select the transcript for each parent; omitting an ambiguous choice raises.
Reference transcript coordinates and reference protein measurements are not
copied onto the altered sequence. Viral parents are not supported for derived
proteins. The snapshot retains the complete derivation record.

## Proteins supplied by other tools

An adapter for an RNA assembly or DNA prediction can submit its amino-acid
sequence without naming an included parent protein. It supplies source claims
as `ProducerInfo`; these do not create genomic `Occurrence` mappings:

```python
from orfdata import ProducerInfo, ProteinSequence

reference = transcriptome.proteome()
rna_product: ProteinSequence = ProteinSequence.from_producer(
    "MPEPTIDE",
    producer=ProducerInfo(
        origin="nc_orf", method="rna_assembly", tool="isovar",
        tool_version="2.0", native_id="ASM42", sample_id="sample-A",
        gene_names=("GENE1",), reference=reference.reference,
        start_complete=False, stop_complete=None,
    ),
)
joint = reference.with_proteins([rna_product])
joint.save("joint.proteome.json.gz")
```

`origin` describes the source protein or ORF category; `method` and `tool`
describe how the supplied sequence was produced. The protein keeps provider,
export paths, sample/run, native record, gene/transcript/ORF IDs, source hash and
known reference when supplied. Its `accession` qualifies the native ID by
producer and sample, so two samples can reuse the same native ID. Equal amino
acid sequences stay distinct entries. Readding the same record is idempotent;
changing content under the same identity raises.

The wrapper that serialized a producer record is provenance, not part of its
source identity. Supply `exporters=("isovar",)` or
`exporters=("vaxrank",)`. Merging direct and wrapped copies produces one
protein and one source record with
`exporters == ("isovar", "vaxrank")`. The singular `exporter=` constructor
argument remains readable for snapshots and callers from 0.35, but new code
should use `exporters`.

Build a standalone producer collection with
`Proteome.from_proteins(proteins, reference=reference.reference)`, then use
`reference.merge(producer_collection)`. If a producer does not know its exact
reference, the caller must explicitly choose
`compatibility="allow_unknown"` when adding it to a known host reference.
A standalone collection with mixed or unknown references uses
`Proteome.from_proteins(proteins, compatibility="mixed")`; merging it into a
host reference also requires an explicit compatibility choice. Known wrong
species or reference declarations always fail. These collection operations
preserve source entries through offline save/load.

An externally supplied fragment may use `taxon_id=None` when its source did
not state a taxon. Its snapshot keeps `None` and its FASTA token says
`taxon=unknown`. Adding it to a known host requires
`compatibility="allow_unknown"`; that policy is stored separately and does not
rewrite the source taxon. A known conflicting taxon still fails.

### Producer records and evidence

`ProteinSequence` carries amino acids. A producer can also attach versioned records
for translations, evidence, events and outcomes that have no sequence. ORFdata
stores these records in the snapshot without importing the producing package or
interpreting its namespaced `payload`:

```python
from orfdata import ProducerLink, ProducerRecord

translation = ProducerRecord(
    namespace="isovar", schema_version="export.v2", kind="translation",
    native_id="T1", sample_id="sample-A", payload={"cdna_id": "C1"},
)
evidence = ProducerRecord(
    namespace="isovar", schema_version="export.v2", kind="evidence",
    native_id="reads-1", sample_id="sample-A", support="rna_sequence",
    scope={"sample": "sample-A", "read_set": "R1"},
    payload={"read_count": 7},
)
unresolved = ProducerRecord(
    namespace="isovar", schema_version="export.v2", kind="event",
    native_id="E2", sample_id="sample-A", status="capped",
    complete=False, limits={"max_hypotheses": 4},
)
joint = joint.with_records(
    [translation, evidence, unresolved],
    links=[
        ProducerLink(translation.key, rna_product.accession, "protein", "translates_to"),
        ProducerLink(evidence.key, translation.key, "record", "supports_translation"),
    ],
)
```

The evidence record says it supports an RNA sequence. Its link to a translation
does not claim that the protein was observed. A producer reporting a protein
observation should use a separate evidence record with an explicit support
type and observation scope. Two translations may link to one protein,
and one evidence record may link to several translations; merging the same
record again does not add its read count again.

An event may have status `zero_product`, `unknown`, `filtered`, `incomplete` or
`capped` without a protein link. `complete` and `limits` retain the producer's
search boundary. These status strings and the versioned payload belong to the
producer; ORFdata retains unknown namespaces and versions exactly. Each record
keeps its native ID separately from its qualified `key`. Links must resolve
to records or proteins in the same proteome. Saving and loading preserve all
records and links, and changing either changes the proteome identity.

## Save and reload

```python
identity = proteome.save("human-e115.proteome.json.gz")
same = Proteome.load("human-e115.proteome.json.gz")  # Offline: no genome or catalogue needed.
assert same.identity == identity
```

A new save uses snapshot schema v5. `load` also reads existing v1–v4
snapshots. Loading v4 producer data normalizes exporter-specific protein and
record keys to the shared source identity and retains the export path.

A snapshot holds every protein with its occurrences and evidence, the layers with
their selections and catalogue identities, and the excluded ORFs. `identity` is
the SHA-256 of that content, so it changes when any sequence, occurrence,
selection, catalogue or exclusion changes, but not with file compression or
formatting. Loading rejects a snapshot whose content does not match its
identity, host proteins from another species, viral proteins without a
matching host layer and layers declaring another reference. The human Ensembl
115 proteome with MS-supported ORFs is a 57 MB snapshot that saves in about 8 s
and loads in about 11 s. Saving is deterministic: the same proteome gives the
same bytes.

### Cache a built proteome

`ReferenceTranscriptome.proteome` delegates to `Proteome.build`. Both default
to reference CDS proteins, then add only the requested ncORFs and viruses.
The equivalent incremental methods are `from_genome(isoforms=...)`,
`with_nc_orfs(nc_orf_sources)` and `with_viral(viral)`.
`viral=` accepts `"all"`/`True`, a virus type, exact bundled reference ID,
`ViralReference`, or a sequence of those specific inputs. `with_viral` also
accepts an external `ProteinSet` when it is accompanied by
`host_taxon_ids=` and `host_source=`; bundled references already carry those
facts.
With `cache=True` it keeps the result in the OpenVax DataCache root
(`ORFDATA_CACHE`, `OPENVAX_DATA_CACHE`, or `cache_root=`). `cache_root=` is
also where the species' registered catalogue is looked for, with or without
`cache=True`:

```python
proteome = transcriptome.proteome(nc_orf_sources=RiboSeq(min_studies=2) & MS(), cache=True)
```

The cache key includes the reference, isoform selection, evidence rule, catalogue and viral
identities, snapshot schema and orfdata version. A later identical build,
from any program, loads that verified snapshot, so tools share one identity.
Earlier bundled-only measurements on Ensembl 115 took about 30 s for the
first build, including saving, and about 7 s for later builds. Loading the
optional datasets adds time and memory before the cache lookup. A cached snapshot is checked
by its DataCache checksum, so it skips the re-serialization `load` does. Local
annotation files of a custom genome are fingerprinted by size and SHA-256, so
genomes that share names never share a protein. Proteins of a cached or loaded
proteome don't carry `orf`; use `extra_orfs(genome)` for ORF objects. Each cached build sits next
to a `recipe.json` recording its key and identity. Nothing is downloaded.
`reference_cds=False` omits reference CDS proteins; `nc_orf_sources=()` is the
reference-only default. Use `proteome(reference_cds=False, nc_orf_sources="all")`
for only ncORFs bound to this transcriptome. `viral="all"` adds every distinct
bundled viral proteome for this host. `nc_orfs`, `catalog`, `catalogs`,
`evidence`, `nc_orf_evidence`, `extra_orf_evidence`, `annotated`, `extra_orfs`
and `host_sources` are no longer constructor arguments. The lower-level
`with_extra_orfs(evidence, catalog=...)` method remains available for existing
callers that query a supplied catalogue with a separate evidence rule.

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

## An annotation for the ncORFs

`write_annotation` writes the proteome's ncORFs as a GTF with the matching
transcript and protein FASTAs, so PyEnsembl-based tools can place a variant
inside an ncORF and translate it:

```python
written = combined.write_annotation("ncorfs")
print(written["annotation_gtf"], written["transcript_count"])
```

```sh
orfdata proteome --species human --release 115 --nc-orf-sources ms \
  --output proteins.fasta --annotation ncorfs/
```

One coding transcript per ncORF is placed on its host gene, named by the
catalogue's `orf_id`, with the ORF's spliced coding sequence and its
translation. Reference CDS entries are left out, because the reference
annotation already describes them, and viral entries are too, because they are
not placed on this genome. [Downstream tools](downstream.md) shows how to read
the result, alone or concatenated with the assembly's own annotation.

## FASTA headers

`write_fasta` writes 60-column sequences and, by default, one FASTA record per
source protein. The first whitespace-delimited token has the same policy for
reference CDS, ncORF, viral and supplied proteins:

```text
>origin=nc_orf|orf=uORF|method=rna_assembly|tool=isovar|gene=GENE1|aa=8|taxon=9606|source=producer%3Aisovar|sample=T1|id=ASM42|key=... kind=external protein_source=nc_orf ...
```

The `key` suffix is a 16-character digest of the source-qualified accession;
it disambiguates records without replacing the readable origin, ORF kind,
method, tool, gene, amino-acid length, source, sample and native ID. Producer
tokens also carry available tool version, provider, export paths and exact
reference, plus `start` and `stop` as `true`, `false` or `unknown`. Optional
`run` and `gene_count` fields appear when needed. `gene=unknown` means no gene
was supplied. Multiple gene names appear as a comma-separated value with
`gene_count`; punctuation and whitespace inside values are percent-escaped
with UTF-8, including `|`, `=`, commas and spaces. `orf=unknown` means the
producer did not state an ORF category. Supplied products may set
`ProducerInfo.orf_kind` (for example `uORF`).

The complete source-qualified accession and original native ID stay in the
snapshot. After `Proteome.load`, `proteome.fasta_index()` maps every FASTA
first token to a tuple containing its source protein. A token collision raises
before writing. `write_fasta(..., deduplicate=True)` groups *exactly equal*
amino-acid sequences; its neutral token gives the record count, length and a
group key without borrowing one producer's metadata, and
`fasta_index(deduplicate=True)` maps it to every contributing protein. This
option never groups a partial sequence with a longer sequence that contains
it. Its description lists every source accession and the distinct protein
sources and taxa. Keep the snapshot beside the FASTA to recover all producer
records, evidence and no-sequence outcomes.

Additional `key=value` fields follow the first token in the description.
All values are percent-escaped, so whitespace, `|`, `=` and newlines in source
metadata cannot add fields or FASTA records. Producer descriptions include
tool version, provider, export paths, native ID, exact reference, source
checksum and explicit start/stop completeness.
`protein_source` uses `reference_cds`, `nc_orf`, `protein_set`, `viral` or the
producer's declared origin. The older `kind` field remains in FASTA and
snapshots for compatibility; `annotated` there means a reference CDS protein
and `extra_orf` means a catalogue ncORF.
For extra ORFs, the native accession is the source ORF ID, `transcript` is its host
transcript in this genome, and `gene`, `gene_name` and `transcript_biotype`
describe that transcript in this genome. `orf_type` is the source's ORF
category, kept verbatim; see [ORF categories](reference.md#orf-categories).
`evidence` lists the supporting evidence types.
`transcripts` counts the encoding transcripts when there are several, and
`source_transcript` names the source's transcript when the header's transcript
is a different one, as for ABLIM1 above, whose source transcript is gone. See
also [identifiers](reference.md#identifiers). The snapshot keeps every
occurrence; the FASTA header shows the first. For bundled viral proteins,
`source` names the exact viral reference; `gene_name` and `product` are the
source's protein labels. A supplied protein FASTA may lack those annotations.
For supplied derived proteins, `derived_from` links to their parent accessions,
and `parent_genes` and `parent_transcripts` label those parents without
asserting that their reference coordinates still describe the changed sequence.
The FASTA header shows one host transcript and a count when several encode an
entry; it does not list every occurrence, evidence record, source checksum or
host-range attribution. Save a `.proteome.json.gz` snapshot alongside FASTA
when those details must be recovered exactly.

## Command line

```sh
orfdata proteome --species human --release 115 --output human-e115.fasta
orfdata proteome --species human --release 115 --isoforms longest --output human-one-per-gene.fasta
orfdata proteome --species human --release 115 --nc-orf-sources all
orfdata proteome --species human --release 115 --nc-orf-sources all --viral all --output human-all-proteins.fasta --snapshot human-all-proteins.proteome.json.gz
orfdata proteome --species human --release 115 --viral HPV16 --output human-hpv16.fasta
orfdata proteome --species mouse --release 115
orfdata proteome --species human --release 115 --no-reference-cds --nc-orf-sources all
orfdata proteome --species human --release 115 --no-reference-cds --viral HPV16
```

The default is reference CDS proteins with all isoforms. `--isoforms longest`
selects one per gene. `--nc-orf-sources` accepts the same string selections as
Python, or a local JSON/JSON.gz dataset file. Repeat the option or use commas
to combine sources:

```sh
orfdata proteome --species human --release 115 --nc-orf-sources bundled,ncorf-mammals-human
orfdata proteome --species human --release 115 --nc-orf-sources ms --nc-orf-sources ribosome
orfdata proteome --species human --release 115 --nc-orf-sources my-orfs.json.gz
```

The separate `--catalog`, `--catalogs` and evidence-selection flags were removed
from this command. For composed evidence rules use the Python argument, e.g.
`nc_orf_sources=RiboSeq(min_studies=2) & MS()`.

`--viral all` adds each distinct bundled viral protein set associated with
the host, choosing its native placement when equivalent references exist.
Use `--viral` with an exact reference ID to choose a placement; repeat the
option for several specific viruses.

The FASTA goes to `--output` or
standard output; `--snapshot FILE` also saves a snapshot, and `--cache`
reuses or keeps the build in the cache; `--cache-root DIR` chooses the cache,
for catalogues as well as builds. A
JSON summary of the layers, identity and excluded ORFs goes to standard error.

## Not yet included

Mouse noncanonical ORFs come from Chang et al.'s and Yang et al.'s Ribo-seq and
OpenProt's MS catalogues once you install them; see [Mouse](catalogs.md#mouse). Dog
noncanonical ORFs, and translated endogenous retrovirus (ERV) ORFs without
annotated transcripts, need new source data; see the
[data roadmap](https://github.com/openvax/orfdata/issues/3). Products from
circular RNAs, stop-codon readthrough or non-AUG starts can be described with
`TranscriptORF.at` and `SourceRNA` (see the [examples](examples.md)), but no
bundled catalogue supplies them yet, so proteomes do not include them.
