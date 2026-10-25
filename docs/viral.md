# Viral coding references for human analysis

## GDC panel coverage

`gdc-viruses` is ORFdata's default 200-reference GDC viral collection for
human transcriptome analysis. Its published source inventory is named
`GRCh38.d1.vd1`. The [published viral
inventory](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt)
contains 200 rows: 189 HPV entries and 11 other viral entries. The published
file name spells `GRCh83`; ORFdata retains that exact source URL and checksum.
Use the bundled, offline inventory to see both available and missing
coding references:

```python
from orfdata import viral_panel_coverage

coverage = viral_panel_coverage()            # one row per screened reference
available = viral_panel_coverage(available=True)
missing = viral_panel_coverage(available=False)
print(len(coverage), len(available), len(missing))  # 200 1 199
print(available.iloc[0]["abbreviation"], available.iloc[0]["bundled_reference_id"])
# HPV16 GRCh38.d1.vd1:HPV16
```

An inventory is a `pandas.DataFrame`, so the arguments below are a convenience
rather than the only way to narrow one — anything pandas does to a frame works
here, and `coverage.attrs` carries the packaged table's own header: which panel
it follows and where it came from.

Use `viral_panel_coverage(available=True)` for bundled references,
`available=False` for rows still missing from ORFdata, `virus="HPV18"` for one
exact GDC abbreviation, or `search="NC_001357"` to search labels, names and
source accessions. The equivalent CLI works offline:

```sh
orfdata viral coverage --available
orfdata viral coverage --missing --virus HPV18
orfdata viral coverage --json
orfdata viral coverage --format tsv > gdc-viral-coverage.tsv
```

To download the **published inventory metadata** into the OpenVax DataCache
for audit purposes, run:

```sh
orfdata viral inventory download
```

The Python equivalent is `download_viral_inventory(cache_root=...)`. Repeated
calls reuse a verified cache entry. This metadata download is optional; the
inventory is bundled and needs no network access.

Each row preserves the GDC inventory's names, abbreviations, source accession
fields and PaVE notes as published. Those accession fields may lack version
numbers; a versioned, sequence-qualified coding reference is reported only
after it is bundled. `source_cds_bundled` means that source CDS models are
available; it does not claim that all reported ORFs or full transcripts are
covered. `not_bundled` describes ORFdata coverage, not the availability of
annotations elsewhere. The inventory contains metadata only and does not make
unavailable entries selectable through `viral_references()` or `viral="all"`.

Beside the inventory, ORFdata packages one **source ledger** per source the
collection draws on, all in the same format: what that source publishes for
every entry in the inventory, with the identifier it resolved to, the
molecule's length and checksums, how many coding features it annotates, its
NCBI taxon and its own description.

| ledger | rows | records |
| --- | --- | --- |
| [`gdc-grch38-d1-vd1-genbank.tsv`](https://github.com/openvax/orfdata/blob/main/orfdata/data/viral/gdc-grch38-d1-vd1-genbank.tsv) | 200 | the versioned NCBI record behind each inventory accession |
| [`gdc-grch38-d1-vd1-pave.tsv`](https://github.com/openvax/orfdata/blob/main/orfdata/data/viral/gdc-grch38-d1-vd1-pave.tsv) | 189 | the PaVE reference the inventory names for each papillomavirus row |

The ledgers hold checksums and counts, not sequence. They are what makes the
decisions below offline ones: which molecule a row should carry, whether a
downloaded record is still the molecule that was pinned, and how much more a
second source annotates. `collection.ledger(source)` reads either one as a
frame keyed by entry abbreviation, and `tools/build_reference_ledger.py
--source NAME` rebuilds it.

```python
from orfdata import viral_collection

ledger = viral_collection("gdc-viruses").ledger("pave")
print(ledger.at["HPV16", "sequence_md5"], int(ledger.at["HPV16", "coding_features"]))
# f2463721818c29998ef5ab917e6df74a 10
```

### What is in each table

Every table ORFdata hands back is a data frame with a fixed, documented column
set. `orfdata.viral_coverage.COVERAGE_COLUMNS`,
`orfdata.viral_coverage.CONTEXT_COLUMNS` and
`orfdata.reference_ledger.FRAME_COLUMNS` are those sets as code, and a test
holds each frame to them.

`viral_panel_coverage()` — one row per screened reference, default integer index:

| column | dtype | meaning |
| --- | --- | --- |
| `index` | Int64 | the row's 1-based position in the published inventory |
| `name` | str | the virus name as the panel publishes it |
| `abbreviation` | str | the panel's short label, unique within a collection |
| `refseq` | str | the RefSeq accession the panel lists, or `""` |
| `genbank` | str | the GenBank accession the panel lists, or `""` |
| `pave_id` | str | the PaVE reference the panel lists, or `""` |
| `pave_notes` | str | the panel's own note about that PaVE reference |
| `bundled_reference_id` | str | the ORFdata coding reference for this row, or `""` |
| `coding_models` | Int64 | how many source CDS models that reference carries |
| `status` | str | `source_cds_bundled` or `not_bundled` |
| `available` | bool | whether `bundled_reference_id` is set |
| `contexts` | object | tuple of context names this row belongs to |

`collection.ledger(source)` — one row per entry the source names, indexed by
`abbreviation`:

| column | dtype | meaning |
| --- | --- | --- |
| `index` | Int64 | that entry's position in the inventory |
| `source` | str | which source this row speaks for, e.g. `genbank` |
| `identifier` | str | what the source publishes the entry under, versioned |
| `state` | str | `published`, `unavailable` or `inconsistent` |
| `length_nt` | Int64 | bases in the molecule the source serves |
| `sequence_md5` | str | MD5 of the uppercased sequence |
| `sequence_sha256` | str | SHA-256 of the same bases |
| `coding_features` | Int64 | how many CDS features the source annotates on it |
| `taxon_id` | Int64 | the NCBI taxon the record cross-references, if any |
| `title` | str | the source's own one-line description of the record |

Only a `published` row pins anything: everything measured off the molecule —
`length_nt`, `sequence_md5`, `sequence_sha256`, `coding_features`, `taxon_id`
— is blank on any other state, and `taxon_id` is also blank where the source
names no taxon.

`viral_contexts()` carries `name`, `title`, `description` and `source` as str,
`references` as the number of inventory rows in that context, and `aliases` as
a tuple of the shorter names it also answers to.

## Download viral genomes, coding transcripts and proteins

To install all 200 GenBank source records for annotation, run:

```sh
orfdata download viral --species human
```

Viral support is part of a regular `pip install orfdata`; only the data is
downloaded. `orfdata setup --species human` installs this default human collection
alongside the human catalogues and chains, and `--collections none` skips it.
The command above installs only this checked bundle under `ORFDATA_CACHE` (or
the normal OpenVax cache), rather than asking for an export directory. Use
`--cache-root DIRECTORY` to select another cache and `--force` to refresh or
repair the installed files.

Here `--species human` is the target analysis species. The current collection
is for human transcriptome analysis; a future collection may target another
species. It does not assert that every source record declared a human host.
The default collection is `gdc-viruses`; pass that value with `--collection`
when an explicit configuration is useful. It names the 200-reference GDC viral
screen, whose published source inventory is `GRCh38.d1.vd1`. The collection's
`cancer-transcriptome-screening` context (with the convenient alias `cancer`)
is descriptive, nonexclusive metadata. It does not classify each virus as
cancer-causing. Future collections can carry several context tags and share
exact accession/version records. This is the only viral download collection
currently published; unsupported species or collections fail with the
available choice.

## How complete is the annotation?

A reference inventory is only as good as the submissions behind it, and those
differ in how much of a genome they annotate. `orfdata viral completeness`
walks the installed collection and says, genome by genome:

```sh
orfdata viral completeness                 # the whole collection
orfdata viral completeness --incomplete    # only genomes with a gap
orfdata viral completeness --subset cancer-causing --json
```

```text
gdc-viruses (GRCh38.d1.vd1): 200 genomes, 1758 coding transcripts
fully annotated: 33
with an unnamed gene: 40
with an incomplete transcript: 49
with a spliced transcript: 40
behind the alternate source: 155 genomes, 235 coding features
initiation: alternative 3, atg 1686, declared 1, non_canonical 2, partial 66
```

`behind the alternate source` is read from the packaged ledger, which records
how many coding features the alternate publishes for each entry, so the
comparison needs no network request — and it is only counted where the
alternate publishes the identical molecule, so it is always a difference of
annotation and never of sequence. No gene list is built in, which is what
keeps the survey true of any collection, species and source;
`viral_annotation_survey()` returns the same thing as data.

What the survey finds in the GDC papillomavirus rows is mostly submission
noise rather than biology. 42 distinct gene labels describe an eight-gene
virus family — `E1^E4`, `E1-E4` and `E1/E4 fusion protein` for one product,
`E5_ALPHA`, `E5A`, `E5B`, `E5_BETA` and `E5` for another, `X`, `Y`, `Z` for
unassigned frames, and `HpV126gp1`…`gp7` for a whole genome. 216 of the 1,367
papillomavirus CDS features carry no `/gene` qualifier at all, leaving the
name in prose: `unnamed protein product; putative E6 ORF`, or
`5' end of the coding region for ORF E6 not determined` on all eight features
of HPV42. The same source annotated at a different time disagrees with itself.

The alternate source does not have this problem, because it is curated rather
than deposited: across the 187 papillomavirus types it publishes, PaVE uses 15
gene labels in all, annotates `E1^E4` and `E8^E2` for every one of them, and
every one of its 1,582 coding features begins on `ATG` and reproduces the
protein it publishes. 184 of the 187 carry the full `E6 E7 E1 E2 E1^E4 L2 L1`
set; the three that do not — HPV101, HPV103 and HPV108 — have no `E6` and an
`E10` instead, which is their genus, not a gap. See
[Whose annotation each molecule carries](#whose-annotation-each-molecule-carries)
for taking it.

## Reading the source records

ORFdata reads GenBank flat files itself, in `orfdata.genbank`, rather than
depending on a parser. The format has many edges, so the reader is compared
against an independent implementation — `gb-io`, a development dependency —
over every record in a downloaded collection:

```sh
pip install 'orfdata[dev]'
ORFDATA_GENBANK_CORPUS=<cache>/viral/references/homo_sapiens/gdc-viruses/<version>/genbank \
  pytest tests/test_genbank_equivalence.py
```

The two agree on all 200 records and 1,678 CDS features. Where they differ,
ORFdata keeps the more faithful reading: `order(...)` stays distinct from
`join(...)`, and `complement(order(>104408,...))` on an EBV gene, which the
other reader returns unparsed, is read properly. Locations ORFdata does not
support — a `123^124` site, an uncertain `(1.5)` position, another record's
coordinates — raise rather than being guessed at.

## Cancer subsets

Three named subsets cover the intended cancer selections:

```sh
orfdata viral subsets
orfdata download viral --subset cancer-integrated-driver  # 15 references
orfdata download viral --subset cancer-causing            # 20 references
orfdata download viral --subset cancer-associated         # 37 references
```

`cancer-causing` is the 20-reference [IARC Group 1](https://monographs.iarc.who.int/agents-classified-by-the-iarc/)
selection. `cancer-associated` contains it, then adds HPV68 (Group 2A), CMV,
HIV-2, and 14 Group 2B HPV references. `cancer-integrated-driver` selects the
12 Group 1 HPV types plus HBV, HTLV-1 and MCV because host-genome integration
and a direct driving role are established in relevant tumour settings.

The remaining 163 GDC references are 162 HPV references outside those IARC
tiers plus SV40. Their exclusion from `cancer-associated` is not a finding that
they lack a cancer connection. Group 3 means IARC cannot classify an agent's
carcinogenicity in humans; it does not mean the agent is safe or noncarcinogenic.

The unfiltered download remains the whole 200-reference GDC screen. Contexts
remain available for future overlapping analysis categories, such as
autoimmune or respiratory, without changing source-reference identity.

The bundle keeps three layers, and names them the same way everywhere:

| Layer | Molecule | Files |
| --- | --- | --- |
| `genome` | DNA | `genomes.fasta`, `genomes.tsv` |
| `transcriptome` | RNA | `transcripts.fasta`, `transcripts.tsv`, `annotation.gtf` |
| `proteome` | protein | `proteins.fasta` |

`genbank/` holds each source record's text exactly as NCBI sent it, gzipped, so
`SHA256SUMS` attests to the source record rather than to a re-serialization of
it. `metadata.json` carries the schema, collection, counts and retrieval time.

A viral genome is not a transcript and a coding transcript is not a protein.
Genome FASTA IDs are the source `accession.version`, which is also the GTF
sequence name. Transcript IDs are `<accession.version>_cds<N>`, one per source
CDS feature, including joined or reverse-strand features. Protein IDs are the
GenBank protein accession when it names exactly one protein in the collection,
and `<transcript_id>.p` otherwise. Each identifier starts from the exact
accession/version rather than from the collection, so importing the same record
through another collection or context keeps its identity.

The transcriptome layer is built from CDS features, which is all these sources
provide: the 200 records carry 1,678 CDS features but only 39 explicit `mRNA`
features across five records. So a "transcript" here is the coding span of one
CDS feature, not a UTR-containing mRNA model, and `transcripts.tsv` records the
source `codon_start`, the translation table, the exon count and whether the
reading frame closes (`annotation_complete`), with `partial` naming the end the
source marks as running past the record.

78 of the 1,703 transcripts are partial, and the reason differs by virus:

- **HPV E4 and E5** (most of them) carry no start codon of their own. E4 is
  translated from the spliced **E1^E4** transcript, whose ATG sits in E1, so
  GenBank writes the E4 ORF as `<3332..3619`. The complete molecule is a
  different transcript, which PaVE annotates for every HPV type.
- **HIV-1 `pol`** is `<1631..4642` because it is translated as a Gag-Pol fusion
  through a −1 ribosomal frameshift, not from its own start codon.
- **EBV's repeat-region genes** (BWRF1 and others) are truncated at the
  boundaries of the collapsed repeats in `AJ507799.2`.
- A few, such as HPV56's E6 `102..>566`, are simply 3'-truncated in the record.

None of these get start or stop codon rows, so effect predictors will not
translate them.

One class is repairable without leaving the record. **HBV's polymerase** is
written as two partial halves, `2307..>3215` and `<1..1623`, because the gene
crosses the origin of the circular genome and a linear record cannot write it
as one span. ORFdata joins them into `join(2307..3215,1..1623)`: 2,532 bases,
843 amino acids, and `genbank:origin_join` as its `annotation_source`. The join
happens only when the result proves itself — whole codons, a start codon, a
terminal stop, no stop in between, and the two published translations meeting
end to end. PyEnsembl reads the joined transcript in order, so varcode
annotates variants on either side of the origin.

Three of the 200 source records annotate no CDS at all: HPV2 (`X55964.1`),
HPV8 (`M12737.1`) and HPV51 (`M62877.1`), the last of them an IARC Group 1
agent in the `cancer-causing` subset. This is not a gap in the GDC inventory,
which names the sequence and its PaVE reference; it is a gap in those GenBank
submissions. NCBI has no RefSeq for HPV8 or HPV51, and no other GenBank record
carries the same sequence annotated.

The [Papillomavirus Episteme](https://pave.niaid.nih.gov/) does annotate them,
on the identical sequence, so ORFdata fills the gap from the link GDC itself
provides — 26 transcripts and proteins across the three — and records
`annotation_source` on every transcript and genome:

```text
annotation sources: genbank 1678, pave 26
```

The fill is used only where the source record annotates nothing, only when
PaVE publishes the identical sequence, and only when each feature's
coordinates reproduce the protein PaVE publishes for it. It adds to the source
records rather than standing in for them, so a PaVE that is unreachable, or
that no longer publishes the same sequence, leaves that genome without
transcripts instead of discarding the download; the summary then names the
genome and the reason. Pass `--no-annotation-fill` (or
`fill_annotation_gaps=False`) to keep strictly to the GenBank submissions.

### Which molecule a collection carries

A reference inventory can name an entry twice — here a GenBank accession and a
PaVE reference — and those are not always the same molecule. The one to carry
is the one the reference the collection follows was built from, and for GDC's
viral decoys that is PaVE's: PaVE's `HPV16REF` reproduces the GDC HPV16 DNA
MD5 pinned in this repository (`f2463721818c29998ef5ab917e6df74a`) exactly,
while `K02718.1`, the accession the inventory lists first, is a different
sequence two bases shorter.

24 of the 189 papillomavirus rows differ this way, among them HPV16, HPV18,
HPV56 and HPV59 — four IARC Group 1 agents. The inventory's own `PaVE Notes`
column marks 18 of them `Revised`. For the other 163 the two are the same
molecule and the question does not arise.

The alternate source's ledger records what each entry looks like there, so the
choice is made offline and the fetched molecule is verified against pinned
checksums before anything is written. Each genome
records `genome_source`, `alternate_source`, `alternate_id` and
`alternate_sequence`, and `orfdata viral summary` reports the mix:

```text
annotation sources: genbank 1520, genbank:origin_join 1, pave 237
genome sources: genbank 176, pave 24
```

The difference this makes is not cosmetic. On its inventory accession, HPV16
carried 7,904 bases and eight coding transcripts with E4 and E5 5'-partial; on
the molecule GDC uses it carries 7,906 bases and ten complete ones, including
the spliced **E1^E4**, **E8^E2** and **E6\*** and an E5 of the canonical 83
residues. Variants called against GRCh38.d1.vd1 are positioned on that
molecule, so this is also what makes their coordinates land.

`--inventory-molecules` keeps every row on the accession its inventory names
first; the summary then names the genomes that are not the alternate
reference. An alternate that is unreachable, or no longer the molecule the
ledger pins, leaves the row on its own record with the reason recorded rather
than discarding the download.

### Whose annotation each molecule carries

A second source can annotate a molecule the collection already has. PaVE
annotates the spliced **E1^E4** and **E8^E2** transcripts for every HPV type,
and the GDC GenBank submissions behind the same types generally do not: of the
189 papillomavirus rows, only three GenBank records annotate a spliced product
at all, and they call it `E1^E4`, `E1-E4` and `E1/E4 fusion protein`.

By default ORFdata changes nothing it does not have to. A row whose own record
annotates no coding feature takes its features from the alternate, when that
is the identical molecule — three rows do. Every other row keeps the
annotation of the submission the inventory names, because that submission is
the record being referenced.

`--alternate-annotations` widens that: wherever the alternate publishes the
identical molecule and annotates more of it, its coding features are used
instead. This never changes a base — the genome FASTA is identical either way
— only which reading of it the transcriptome and proteome record. Across the
GDC collection it moves 155 genomes and adds 235 coding features, nearly all
of them the spliced papillomavirus transcripts, at the cost of replacing each
submission's own gene naming with one source's. `orfdata viral completeness`
reports what it would change before you ask for it.

Spliced models for a papillomavirus type whose genome also came from
elsewhere remain a different question: that is a different collection, which
is a new entry in the registry rather than a new code path.

### How translation starts, and where it is not straightforward

The transcriptome layer records how each coding sequence begins, because a
source's own claim is the only reference for it. `transcripts.tsv` carries
`codon_start` — how many bases of the feature precede the first codon, which
the GTF frame column repeats so readers stay in frame — alongside
`start_codon`, the three bases actually read, and `initiation`, which is one
of:

| `initiation` | meaning |
| --- | --- |
| `atg` | the canonical initiation codon |
| `alternative` | another initiation codon of the same genetic code, such as `CTG` |
| `declared` | neither, and the source declares the initiator with `transl_except` |
| `non_canonical` | not an initiation codon, and nothing declares it |
| `partial` | the 5' end runs past the record, so there is no start codon to read |
| `unknown` | the source names a genetic code ORFdata does not carry |

Across the GDC collection: `atg` 1,686, `partial` 66, `alternative` 3,
`non_canonical` 2, `declared` 1. The declared one is CMV's `UL30A`, an `ACG`
start carrying `/transl_except=(pos:complement(38131..38133),aa:Met)` and
ribosome-profiling evidence. The two `non_canonical` ones are HPV63 and HPV65,
whose submitters say in a note that the protein is in conflict with the
conceptual translation — the classification finds them without being told.

**The first residue belongs to the proteome, not the transcriptome.** A source
writes the first residue of a coding sequence as methionine whatever codon
initiates it, so `proteins.fasta` is the reference for what the protein starts
with and the transcript table never second-guesses it. `agrees_with` applies
that convention when checking a source's own translation, and nothing else
relies on it.

An **IRES** is orthogonal to all of this. It governs how the ribosome is
recruited to a message, not which codon is read first or what it codes for, so
a CDS translated from one is annotated like any other and comes through
unchanged. No record in the GDC collection annotates an IRES feature, so none
is lost; were one annotated, it would belong to the untranslated region of the
transcript rather than to its coding sequence.

Where the ribosome does something a codon-by-codon reading of the genome
cannot reproduce, the `recoding` column names it. Two transcripts in the
collection carry `ribosomal_slippage`: HTLV-1's `gag-pro` and `gag-pro-pol`,
written as joins that read the slippage base twice. Their sequences and
proteins are correct, and varcode reads them per product — but a variant at
the reused base itself is reported as an exon loss, because a GTF exon is a
genomic interval and effect predictors assume exons do not overlap.

Because the GTF, transcript FASTA and protein FASTA share identifiers, the set
opens as a PyEnsembl genome — see [Downstream tools](downstream.md):

```python
from orfdata import viral_genome

genome = viral_genome(bundle)
print(genome.transcript_by_id("AY446894.2_cds1").protein_sequence)
```

The Python API returns the generated paths and record counts:

```python
from orfdata import (
    download_viral_references,
    index_viral_references,
    viral_reference_bundle,
)

bundle = viral_reference_bundle(
    download_if_missing=True,
    species="human",
)
print(bundle.genome_fasta)
print(bundle.transcript_fasta, bundle.protein_fasta, bundle.annotation_gtf)
print(bundle.genome_count, bundle.protein_count)  # 200, 1678 for the pinned records
print(bundle.collection, bundle.contexts)

for layer in bundle.layers:  # genome, transcriptome, proteome
    print(layer.name, layer.molecule, layer.records, layer.fasta)

# Later, reopen the verified cached files without a network request.
bundle = viral_reference_bundle(species="human")

# The explicit downloader is equivalent to download_if_missing=True.
bundle = download_viral_references(species="human")

# A named subset index reuses the source bundle and can fetch it on demand.
cancer_index = index_viral_references(
    subsets=["cancer-causing"],
    species="human",
    download_if_missing=True,
)
print(cancer_index.protein_fasta)
```

Summarize what is installed, or what a selection of it holds, without writing
anything:

```sh
orfdata viral summary
orfdata viral summary --subset cancer-causing --json
```

```python
from orfdata import viral_reference_summary

summary = viral_reference_summary(subsets=["cancer-causing"])
print(summary["genome_count"], summary["transcript_count"])
print(summary["complete_transcript_count"])
```

Discover this collection alongside every catalogue, quantitative dataset,
pinned genome and chain file, without downloading anything:

```sh
orfdata download list
orfdata download list --kind viral --species human
orfdata download status
orfdata download list --missing
```

Use `--json` for machine-readable output and `--verify` to hash installed files.
The Python equivalent is `list_downloads()`.

This is a source-record download for alignment and annotation workflows. It
does not change the curated, sequence-qualified references exposed through
`viral_references()` or selected by `Proteome.build(..., viral=...)`.

The [viral panel roadmap](https://github.com/openvax/orfdata/issues/91) tracks
the remaining imports. GDC panel membership and a documented host association
are separate facts; `viral="all"` still selects bundled references associated
with the proteome's host.

Browse the packaged viral references by host, convenience group and virus type,
then load their source-annotated coding paths without fetching the genome:

```python
from orfdata import viral_orfs, viral_reference, viral_references

references = viral_references(host="human", group="HPV", virus_type="HPV16")
print([(r.id, r.coding_models) for r in references])
# [('NC_001526.4', 9), ('GRCh38.d1.vd1:HPV16', 9)]

ref = viral_reference("NC_001526.4")  # choose the native RefSeq placement
reference_dna = ref.dna
for orf in ref.orfs:
    print(orf.orf_id, orf.gene_name, orf.coding_blocks, orf.protein_sequence)
```

The small bundled catalogue loads offline without PyEnsembl, a GenBank
parser or a human genome. It supplies the complete DNA for the native RefSeq NC_001526.4
and GDC HPV16 placements; both map to the same nine source CDSs. The two DNA
strings differ by a verified circular rotation. `viral_reference_sequence(id)`
returns the exact forward FASTA bases, with T rather than U. A missing ID raises
`ValueError`. `viral_references()` returns immutable `ViralReference` data
classes; `.dna` and `.orfs` load only on access, so browsing does not copy
genome bases or construct ORF models for every result. Choose by reference ID,
not by result position: the native and GDC sequences have different coordinates.
Use `ref.to_dict()` when serializing reference metadata to JSON; it does not
load DNA or ORF models.

For a combined protein FASTA, callers can pass the virus type or exact
reference ID directly to `Proteome.build` without constructing a separate
viral input object:

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

genome = EnsemblRelease(115, species="human")  # Already installed with pyensembl install.
combined = Proteome.build(genome, nc_orf_sources="all", viral="HPV16")
combined.write_fasta("human-hpv16.fasta")
```

`viral="all"` includes each distinct bundled viral protein set associated
with the host. When equivalent references differ only in DNA coordinate
placement, it chooses the native source ID. The viral FASTA header records
the chosen exact reference ID, source gene and product label; the snapshot
keeps its source checksum and host-range attribution. See
[Build a reference proteome](proteome.md#viral-proteomes).

`host` accepts a known species name such as `"human"` or an NCBI taxon ID such
as `9606`. The released HPV16 reference objects have viral `taxon_id=333760`
and `host_taxon_ids=(9606,)`; these are different identities. Their documented
human host association is sourced to PaVE. `ref.group == "HPV"`
is a convenient selection label, not a formal taxonomic rank. The
`ref.virus_type == "HPV16"` label is the virus
type; ICTV places it in the species *Alphapapillomavirus 9*. Each reference also
retains the association and taxonomy source URLs. Host associations can cover
several host taxa as new references are curated; an unlisted host is unknown,
not proof that infection is impossible. NCBI Virus's isolation-host metadata is
not an exhaustive host-range list. `ref.host_taxon_ids` and `ref.host_source` have the
same meaning as in [ViralProteome](proteome.md#viral-proteomes).

When annotating an existing BAM, continue to bind against its **actual** full
alignment-reference DNA:

```python
orfs = viral_orfs(alignment_reference_dna, contig="HPV16")
```

This checks the sequence hash and its qualified coordinates. A contig name or
related strain cannot transfer CDS annotations to mismatching DNA. It raises
`ValueError` for an unavailable reference rather than returning an empty ORF
set. `viral_references()` reports coverage and source checksums, and
`viral_orfs()` returns the existing `SourceRNA`/`TranscriptORF` models.

Each model uses a `SourceRNA` containing only the available CDS, including
ordered splice blocks. Full transcript ends are not asserted. The existing
`TranscriptORF.at()` validates its source protein and coding boundaries.
These models have `orf_type="reference_CDS"`. Their internal kind `extra_orf`
describes a source-defined ORF outside a loaded PyEnsembl annotation, not a
claim that a RefSeq CDS is biologically noncanonical. No Ribo-seq or MS
support is inferred from the reference annotation.

The nine annotations are not an exhaustive HPV16 ORF catalogue. Other HPV
types/lineages, EBV/KSHV, CMV, hepatitis viruses including HDV, HTLV-1/HIV,
and polyomaviruses remain explicit coverage work in
[issues #91–98](https://github.com/openvax/orfdata/issues/91).

Isovar owns reconstructed sample RNA and read evidence. It can use these
coding paths to anchor a frame on observed sequence, preserving variants
and incomplete starts/stops. ORFdata does not assemble sample RNA, replace
missing bases, infer protein abundance or establish viral tumor causality.
RNA and coding-reference matches are not direct protein evidence.

Source attribution and rebuilding instructions are in
[`orfdata/data/viral/NOTICE.md`](https://github.com/openvax/orfdata/blob/main/orfdata/data/viral/NOTICE.md).
