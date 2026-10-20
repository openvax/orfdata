# Build a reference proteome

A `ReferenceTranscriptome` names the transcript models in one installed
PyEnsembl annotation. A `Proteome` selects protein products from those models:
translations of the reference CDSs, catalogue ncORFs bound to the transcripts,
and optionally viral proteins associated with the host. All three sources have
annotations; `reference_cds` identifies the proteins induced by the reference
annotation's CDS features. Each method returns a new proteome. `proteome.layers` records
every addition with its role (`host`, `viral` or `filter`), source, selection
and counts, and each entry lists every transcript that encodes it.

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome, ReferenceTranscriptome

genome = EnsemblRelease(115, species="human")  # Installed with pyensembl install.
transcriptome = ReferenceTranscriptome.from_genome(genome)
proteome = transcriptome.proteome(nc_orf_evidence="ms")
reference_cds = transcriptome.proteome(nc_orfs=False)
```

The examples on this page continue that Python session.

## Reference CDS proteins

`Proteome.from_genome(genome)` and `nc_orfs=False` take every
transcript with a reference CDS protein sequence. For Ensembl this matches the annotation's protein FASTA: Ensembl 115
has 245,535 human proteins with 160,811 distinct sequences. They include
nonsense-mediated-decay and immunoglobulin/T-cell receptor segment translations,
and 11,540 proteins on patch, haplotype and unplaced contigs. Each entry
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

## Choose the evidence

`with_extra_orfs(evidence=...)` selects catalogue ORFs by evidence and adds each
one that a transcript of this genome encodes: the transcript contains the ORF's
exact spliced coding path and translates to its protein. Counts below are for
the bundled human catalogue (dataset `2026-10-11.1`) and Ensembl 115:

| Selection | Human ORFs | Encoded in Ensembl 115 |
| --- | ---: | ---: |
| `evidence=None` or `"all"`: every catalogue ORF | 7,264 | 7,239 |
| `evidence="ribosome"` or `"riboseq"`: translated in ribosome profiling | 7,264 | 7,239 |
| `evidence="ms"` (default): any supporting MS peptide | 2,186 | 2,181 |
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

A list (or tuple or set) of names selects ORFs with any of them:
`evidence=["proteomics", "immunopeptidomics"]` keeps ORFs with either kind of
MS peptide, unlike `"ms"`, which also counts unclassified MS. For ORFs that
need two kinds of evidence at once, combine rules with `&`, below.

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

reference_cds.with_extra_orfs(RiboSeq(min_studies=2))                        # replicated Ribo-seq
reference_cds.with_extra_orfs(RiboSeq(min_studies=2) & MS())                 # Ribo-seq and MS
reference_cds.with_extra_orfs(RiboSeq(min_studies=2) | MS())                 # either
reference_cds.with_extra_orfs(RiboSeq() & MS(mhc=True))                     # Ribo-seq and an MHC ligand
reference_cds.with_extra_orfs(MS(methods=("proteomics", "ms_unclassified")))  # not immunopeptidomics-only
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

Without `catalog=`, `load_catalog` picks the catalogue for the genome's species
and assembly: a bundled one, or else one installed through
[`CatalogRegistry`](catalogs.md), several merged. Human GRCh38 has the bundled
one, and mouse GRCm39 three once [installed](catalogs.md#mouse); for
other species the error names the catalogues that exist. Pass `catalog=` to use
another, such as `load_catalog(path)`. A catalogue for a different species or
assembly is rejected.

### Merging catalogues

Pass several catalogues for one species and assembly to use their ORFs together:

```python
from orfdata import load_catalog, merge_catalogs

catalogs = [load_catalog(), load_catalog("study-orfs.json.gz")]  # The bundled one and a local file.
proteome = reference_cds.with_extra_orfs("ms", catalog=catalogs)  # Or merge_catalogs(catalogs).
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
`--catalog FILE`.

## Viral proteomes

Build one proteome with reference CDS proteins, every catalogue ORF that binds to
the genome, and every distinct bundled viral protein set recorded for this host:

```python
combined = transcriptome.proteome(nc_orf_evidence="all", viral="all")
proteins_by_id = combined.by_accession()
combined.write_fasta("human-all-proteins.fasta")
combined.save("human-all-proteins.proteome.json.gz")
```

`nc_orf_evidence` filters only the ncORF catalogue added to the reference CDS
proteome. `"all"` selects every catalogue ORF; those that do not bind to this
genome are reported in `combined.excluded`. Reference CDS and viral proteins are
unaffected by that option. The older `extra_orf_evidence=` and `evidence=`
spellings still work on `Proteome.build`.

`proteins_by_id` maps each unique protein accession to a `ProteomeEntry` with
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

`ViralProteome` is the input adapter for a protein FASTA supplied outside the
bundled references. It holds that virus's NCBI taxon, recorded host taxa and
host-range source, so the same host check applies:

```python
from orfdata import ViralProteome

# Template: substitute your FASTA, the virus taxon and your host-range source.
hpv16 = ViralProteome.from_fasta(
    "hpv16-proteins.fasta", name="HPV16", taxon_id=333760,
    host_taxon_ids={9606}, host_source="Host range source and date",
)
with_virus = proteome.with_viral(hpv16)
```

`viral=[first, second]` accepts several specific inputs. The tests exercise
this with two synthetic viruses; HPV16 is the only released biological viral
reference at present.

All three protein kinds use the same `ProteomeEntry` and appear in the same
`Proteome`. The distinct `ViralProteome` type describes an external input,
not a second output collection. Viral entries keep the virus's taxon and a
source label such as `viral:HPV16@NC_001526.4`; their layer records the source
accession, checksum and host-range source. `proteome.host()` returns the host's
own proteins without viral entries, for example as a self reference.

### ncORFs with or without an installed transcriptome

If you have a PyEnsembl reference transcriptome, build a proteome containing
only the ncORFs that bind to its transcripts:

```python
nc_only = transcriptome.nc_orf_proteome(evidence="all")
nc_only.write_fasta("reference-ncorfs.fasta")
```

```sh
orfdata proteome --species human --release 115 --no-reference-cds --nc-orf-evidence all --output reference-ncorfs.fasta
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

The parent must already be in the proteome. A derived entry has `kind="derived"`
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

## Save and reload

```python
identity = proteome.save("human-e115.proteome.json.gz")
same = Proteome.load("human-e115.proteome.json.gz")  # Offline: no genome or catalogue needed.
assert same.identity == identity
```

A new save uses snapshot schema v2 for derivation metadata; `load` also reads
existing v1 snapshots.

A snapshot holds every entry with its occurrences and evidence, the layers with
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

`ReferenceTranscriptome.proteome` delegates to `Proteome.build`, which makes
the usual proteome in one call: `from_genome`, then
`with_extra_orfs(nc_orf_evidence, catalog=...)` and `with_viral` for each of `viral=`.
`viral=` accepts `"all"`/`True`, a virus type, exact bundled reference ID,
`ViralReference`, external `ViralProteome`, or a sequence of specific inputs.
With `cache=True` it keeps the result in the OpenVax DataCache root
(`ORFDATA_CACHE`, `OPENVAX_DATA_CACHE`, or `cache_root=`). `cache_root=` is
also where the species' registered catalogue is looked for, with or without
`cache=True`:

```python
proteome = transcriptome.proteome(nc_orf_evidence=RiboSeq(min_studies=2) & MS(), cache=True)
```

The cache key is the reference, evidence rule, catalogue and viral
identities, snapshot schema and orfdata version. A later identical build,
from any program, loads that verified snapshot, so tools share one identity.
For human Ensembl 115 with MS-supported ORFs, the first build takes about
30 s, including saving, and later ones about 7 s. A cached snapshot is checked
by its DataCache checksum, so it skips the re-serialization `load` does. Local
annotation files of a custom genome are fingerprinted by size and SHA-256, so
genomes that share names never share an entry. Entries of a cached or loaded
proteome don't carry `orf`; use `extra_orfs(genome)` for ORF objects. Each cached build sits next
to a `recipe.json` recording its key and identity. Nothing is downloaded.
`reference_cds=False` omits reference CDS proteins; `nc_orfs=False` omits
catalogue ncORFs. Setting both to False with `viral=` builds a viral-only
proteome for this host. `transcriptome.nc_orf_proteome(evidence="all")` is the
short form for the ncORF-only case. The older `annotated=`, `extra_orfs=` and
`host_sources=` selectors still work for existing callers.

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
>ENSP00000269305.4 kind=annotated protein_source=reference_cds taxon=9606 source=ensembl:115 gene=ENSG00000141510 gene_name=TP53 transcript=ENST00000269305.9 contig=17 transcript_biotype=protein_coding
>c7riboseqorf56 kind=extra_orf protein_source=nc_orf taxon=9606 source=orfdata-catalog:2026-10-11.1 gene=ENSG00000169902 gene_name=TPST1 transcript=ENST00000304842.6 contig=7 transcript_biotype=protein_coding orf_type=uORF evidence=ribosome,immunopeptidomics transcripts=18
>c10riboseqorf115 kind=extra_orf protein_source=nc_orf taxon=9606 source=orfdata-catalog:2026-10-11.1 gene=ENSG00000099204 gene_name=ABLIM1 transcript=ENST00000369253.6 contig=10 transcript_biotype=protein_coding orf_type=intORF evidence=ribosome,immunopeptidomics source_transcript=ENST00000277895 transcripts=10
>NP_041327.2 kind=viral protein_source=viral taxon=333760 source=viral:HPV16@NC_001526.4 gene_name=E1 product=replication_protein_E1
```

`protein_source` uses `reference_cds`, `nc_orf` or `viral`. The older `kind`
field remains in FASTA and snapshots for compatibility; `annotated` there
means a reference CDS protein and `extra_orf` means a catalogue ncORF.
For extra ORFs, the accession is the source ORF ID, `transcript` is its host
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
orfdata proteome --species human --release 115 --nc-orf-evidence all
orfdata proteome --species human --release 115 --nc-orf-evidence all --viral all --output human-all-proteins.fasta --snapshot human-all-proteins.proteome.json.gz
orfdata proteome --species human --release 115 --viral HPV16 --output human-hpv16.fasta
orfdata proteome --species mouse --release 115 --no-nc-orfs
orfdata proteome --species human --release 115 --no-reference-cds --nc-orf-evidence all
orfdata proteome --species human --release 115 --no-reference-cds --no-nc-orfs --viral HPV16
```

`--nc-orf-evidence` accepts the evidence names above, `all` for every
catalogue ORF, or `none` for reference CDS proteins only. `--require-unique`,
`--require-reviewed`, `--require-resolved` and `--catalog FILE` match the Python
options. The earlier `--extra-orf-evidence` spelling remains accepted.

`--viral all` adds each distinct bundled viral protein set associated with
the host, choosing its native placement when equivalent references exist.
Use `--viral` with an exact reference ID to choose a placement; repeat the
option for several specific viruses.

Evidence rules have their own options:

```sh
orfdata proteome --species human --release 115 --min-ribo-seq-studies 2 --min-ms-peptides 1
orfdata proteome --species human --release 115 --min-ribo-seq-studies 3 --min-ms-peptides 2 --either
orfdata proteome --species human --release 115 --min-ms-peptides 1 --ms-methods proteomics,ms_unclassified
```

Given both a Ribo-seq and an MS requirement, an ORF needs both unless
`--either` is set. `--require-unique` and `--require-reviewed` refine the MS
requirement, and `--require-resolved` both. The FASTA goes to `--output` or
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
