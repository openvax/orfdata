# Viral coding references

## GDC panel coverage

The GDC reference is named `GRCh38.d1.vd1`. Its [published viral
inventory](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt)
contains 200 rows: 189 HPV entries and 11 other viral entries. The published
file name spells `GRCh83`; ORFdata retains that exact source URL and checksum.
Use the bundled, offline coverage ledger to see both available and missing
coding references:

```python
from orfdata import viral_panel_coverage

coverage = viral_panel_coverage()
available = viral_panel_coverage(available=True)
missing = viral_panel_coverage(available=False)
print(len(coverage), len(available), len(missing))  # 200 1 199
print(available[0].abbreviation, available[0].bundled_reference_id)
# HPV16 GRCh38.d1.vd1:HPV16
```

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
coverage ledger is bundled and needs no network access.

Each row preserves the GDC inventory's names, abbreviations, source accession
fields and PaVE notes as published. Those accession fields may lack version
numbers; a versioned, sequence-qualified coding reference is reported only
after it is bundled. `source_cds_bundled` means that source CDS models are
available; it does not claim that all reported ORFs or full transcripts are
covered. `not_bundled` describes ORFdata coverage, not the availability of
annotations elsewhere. The ledger contains metadata only and does not make
unavailable entries selectable through `viral_references()` or `viral="all"`.

ORFdata also records NCBI ESummary metadata for the 200 GenBank accession
leads in
[`gdc-grch38-d1-vd1-genbank-summaries.tsv`](../orfdata/data/viral/gdc-grch38-d1-vd1-genbank-summaries.tsv):
current accession versions, titles, sequence lengths, taxon IDs and record
dates. This bundled metadata lookup matched all 200 inventory accessions; it
does not itself contain sequence payloads or CDS annotations.

## Download source genomes, coding sequences and proteins

To install all 200 GenBank source records for annotation, run:

```sh
orfdata download viral --species human
```

The downloader and its Biopython/DataCache dependencies are part of a regular
`pip install orfdata`; the sequence payload stays on demand. The command
installs a checked bundle under `ORFDATA_CACHE` (or the normal OpenVax cache),
rather than asking for an export directory. Use `--cache-root DIRECTORY` to
select another cache and `--force` to refresh or repair the installed files.
`orfdata viral download` remains a compatibility alias.

Here `--species human` is the target analysis species. It does not assert that
every source record declared a human host. The default collection is
`gdc-grch38-d1-vd1`; pass that value with `--collection` when an explicit
configuration is useful. The collection's
`cancer-transcriptome-screening` context (with the convenient alias `cancer`)
is descriptive, nonexclusive metadata.
It does not classify each virus as cancer-causing. A user may request CMV or EBV
for both cancer and autoimmune workflows, so a disease context is not a mutually
exclusive download selector or part of cache identity. Future collections can
carry several context tags and share exact accession/version records. This is
the only viral download collection currently published; unsupported species or
collections fail with the available choice.

`cancer-causing` selects the 20 viral-agent references represented
here that [IARC classifies as Group 1, carcinogenic to humans](https://monographs.iarc.who.int/agents-classified-by-the-iarc/).
`cancer-associated` selects the 37 references classified by IARC in Groups 1,
2A or 2B; it includes the Group 1 set, HPV68, CMV, HIV-2, and HPV types with
weaker or conditional IARC classifications. These selectors describe evidence
tiers, rather than declaring that the remaining GDC screening entries lack a
cancer connection.

## Tumour mechanism selectors

Contexts describe why a collection may be useful. Mechanism selectors describe
an evidence-backed virus--tumour relationship and can be combined to make a
derived FASTA index:

```sh
orfdata viral mechanisms
# 12 Group 1 HPV types, HBV, HTLV-1 and MCV: 15 references
orfdata download viral \
  --mechanism host-genome-integrating \
  --mechanism direct-tumor-cell-driver \
  --mechanism-match all
```

`host-genome-integrating` also includes HIV-1 and HIV-2, because retroviral
integration is part of their lifecycle. They are removed by the
`direct-tumor-cell-driver` intersection: HIV-associated cancer risk commonly
involves immune dysfunction and other oncogenic viruses rather than HIV as a
driver in the tumour cell. EBV and KSHV are direct drivers in relevant tumours
but normally persist as episomes, so the same intersection excludes them.

`viral-oncogene-dependence-evidence` is deliberately narrower: it currently
contains HPV16, HPV18 and MCV, for which direct experimental evidence supports
continued viral-oncogene expression in named virus-positive tumour settings.
It is not a claim that every tumour carrying the virus is dependent, nor a
prediction that antiviral treatment will eradicate the cancer. The exact
tumour type, viral integration, viral expression and clonal distribution must
be resolved in the sample. The mechanism table links the supporting primary
studies for HPV integration and E6/E7 dependence, HBV integration, HTLV-1
proviral integration, and MCV integration and T-antigen dependence.

List registered contexts with `orfdata viral contexts`. Supplying one or more
contexts creates a derived index from the shared collection: repeated
`--context` values select a union and `--context-match all` selects an
intersection. The index has its own checked FASTAs but reuses the same cached
source records, so a reference shared by cancer and autoimmune contexts keeps
one accession-based identity.

```sh
orfdata download viral --context cancer  # all 200 GDC screening references
orfdata download viral --context cancer-associated  # IARC Groups 1, 2A and 2B
orfdata download viral --context cancer-causing  # `iarc-group-1` is an alias
# Excludes the Group 1 subset from the GDC screen; it does not assert the rest are safe.
orfdata download viral --exclude-context cancer-causing
# Once another listed context is published, select its union or intersection:
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context>
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context> --context-match all
```

The bundle contains `genomes.fasta`, `coding_sequences.fasta`,
`proteins.fasta`, `cds.tsv`, `manifest.tsv`, individual gzipped GenBank records,
and `SHA256SUMS`. The version-pinned accession summary selects the exact NCBI
records. Each genome FASTA record starts with `viral:<accession.version>`. The
coding-sequence and protein FASTAs contain one paired entry per source CDS
feature and use the same FASTA identifier, including for joined or
reverse-strand CDS features. The identifier starts with the exact
accession/version rather than the collection, so importing the same record
through another collection or context retains its sequence and product
identity.

PyEnsembl calls these protein FASTAs rather than proteome FASTAs. Its transcript
FASTA contains spliced cDNA transcripts, including UTRs. The GDC source records
have 1,678 CDS features but only 39 explicit `mRNA` features across five of the
200 records. ORFdata therefore calls the nucleotide export
`coding_sequences.fasta`; calling it a transcriptome would imply transcript
models that the sources do not provide.

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
print(bundle.coding_sequence_fasta, bundle.protein_fasta)
print(bundle.genome_count, bundle.protein_count)  # 200, 1678 for the pinned records
print(bundle.collection, bundle.contexts)

# Later, reopen the verified cached files without a network request.
bundle = viral_reference_bundle(species="human")

# The explicit downloader is equivalent to download_if_missing=True.
bundle = download_viral_references(species="human")

# A context index reuses the source bundle and can fetch it on demand.
cancer_index = index_viral_references(
    contexts=["cancer-causing"],
    species="human",
    download_if_missing=True,
)
print(cancer_index.protein_fasta)
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

The small bundled catalogue loads offline without PyEnsembl, Biopython or a
human genome. It supplies the complete DNA for the native RefSeq NC_001526.4
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
