# Discover and download data

ORFdata has one common entry point for data acquisition:

```sh
orfdata download list
orfdata download status
orfdata download viral --species human
orfdata download catalog NAME
orfdata download dataset NAME
orfdata download reference DATASET
orfdata download genome NAME
orfdata download chains --species human
```

The operation comes first so the same discovery and installation model applies
to every resource kind. `orfdata download viral` is therefore the primary
command. Focused command groups such as `orfdata viral coverage` and
`orfdata viral summary` remain useful for domain-specific inspection.

`list` shows every bundled, downloadable or locally buildable primary resource
and its current state, including native protein search references declared by
quantitative datasets. `status` shows resources already usable from local
files. Filter either command with `--kind` or `--species`; use `--installed` or
`--missing` with `list`, `--json` for programs, and `--verify` to hash installed
files. Discovery and status commands never acquire or repair data.

```sh
orfdata download list --kind viral --species human
orfdata download list --missing
orfdata download status --json
```

Catalogue and dataset releases that cannot be redistributed are built from
their pinned sources when selected through `orfdata download`. The listing's
`ACQUIRE` column says `download` or `build`. Existing focused commands such as
`orfdata catalogs download`, `orfdata datasets build`, `orfdata genomes
install`, and `orfdata chains download` remain available.

The same inventory is available from Python:

```python
from orfdata import list_downloads

all_resources = list_downloads()
human_missing = list_downloads(species="human", installed=False)
installed = list_downloads(installed=True)
machine_record = all_resources[0].to_dict()
```

Each `DownloadStatus` reports its kind, stable name and version, species,
nonexclusive analysis contexts, local status, acquisition method and equivalent
command, and resolved cache path. Status values follow DataCache: `available`,
`missing`, `corrupt`, or `inaccessible`. An unknown kind, species or context is
an error rather than an empty listing, so a typo never reads as "nothing is
available".

## Viral collections and contexts

The published viral collection currently has these properties:

| Property | Value | Meaning |
| --- | --- | --- |
| Species | `human` | Current target analysis species, NCBI taxon 9606 |
| Collection | `gdc-viruses` | Default 200-reference GDC viral screen for human analysis |
| Contexts | `cancer-transcriptome-screening` (`cancer`) | GDC human cancer-transcriptome screening context |

Run the default collection with:

```sh
orfdata download viral --species human
```

Viral support is part of a regular `pip install orfdata`; what needs acquiring
is the data. `orfdata setup --species human` installs this collection along
with the human catalogues and chains; `--collections none` sets up catalogues
and chains alone.

`--collection gdc-viruses` is the default and is available for explicit
configuration in scripts that will later choose among several collections. The
current collection is for human analysis. A future collection may target another
species. Target species does not assert a human source host for every virus.
Context is metadata rather than a sequence identity or a mutually exclusive
selector. CMV and EBV, for example, may be requested in both cancer and
autoimmune workflows. A future collection may carry several context tags, and
the same exact accession/version can belong to several collections without
becoming several biological references.

Pass `--context` to materialize a derived FASTA index from the shared source
bundle. Repeated contexts select their union; `--context-match all` selects
their intersection. It reuses source records already downloaded for the
collection, so overlapping contexts do not create conflicting copies of a
reference.

```sh
orfdata viral contexts
orfdata download viral --context cancer  # all 200 GDC screening references
# Once another listed context is published, select its union or intersection:
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context>
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context> --context-match all
```

The target species is an admission and organization choice for the reference
set. It does not replace a source record's viral taxon or assert that the source
declared a host. Unknown source-host metadata remains unknown.

The four facts stay separate in the downloaded tables and FASTA headers:
the versioned viral accession and CDS feature identify the source product;
collection and context record where it was selected; the target species records
the analysis it supports; and a source-host taxon is populated only when the
source establishes one. This keeps the same reference idempotent when it is
obtained through several collections or exporter paths.

The published GDC source inventory is named `GRCh38.d1.vd1`; `gdc-viruses` is
the ORFdata collection name for its 200 references for human analysis.

Use the full GDC collection for broad screening. For its three cancer
selections, list named subsets with `orfdata viral subsets`:

```sh
orfdata download viral --subset cancer-integrated-driver
orfdata download viral --subset cancer-causing
orfdata download viral --subset cancer-associated
```

`cancer-causing` is the 20-reference IARC Group 1 selection.
`cancer-associated` is its 37-reference Group 1, 2A and 2B superset.
`cancer-integrated-driver` contains the 12 Group 1 HPV types, HBV, HTLV-1 and
MCV. The remaining 163 GDC entries are 162 HPV references outside those IARC
tiers plus SV40; their absence from the cancer-associated set does not establish
that they are biologically irrelevant. IARC Group 3 means an agent is not
classifiable as to carcinogenicity in humans; it is not a safety finding. Future
respiratory and autoimmune contexts can overlap the GDC screening context.

## Three layers, named separately

A downloaded collection keeps the molecular layers apart, and ORFdata uses the
same three words everywhere for them:

| Layer | Molecule | Files | What it is |
| --- | --- | --- | --- |
| `genome` | DNA | `genomes.fasta`, `genomes.tsv` | Viral genomic DNA: the same molecules a DNA reference such as `GRCh38.d1.vd1` carries as extra contigs |
| `transcriptome` | RNA | `transcripts.fasta`, `transcripts.tsv`, `annotation.gtf` | One coding transcript per source-annotated CDS feature, spliced, placed on the genome |
| `proteome` | protein | `proteins.fasta` | The source translation of each coding transcript |

A viral genome is not a transcript and a coding transcript is not a protein, so
nothing in the cache mixes them: each layer has its own FASTA, its own
identifiers and its own record count. Every transcript also records which
source annotated it, and which end, if either, that source marks as running
past the record.

```python
from orfdata import viral_reference_bundle

bundle = viral_reference_bundle(species="human")
for layer in bundle.layers:
    print(layer.name, layer.molecule, layer.records, layer.fasta)
transcriptome = bundle.layer("transcriptome")
print(transcriptome.annotation)  # annotation.gtf
```

Summarize an installed collection, or one selection of it, without writing
anything:

```sh
orfdata viral summary
orfdata viral summary --subset cancer-causing --json
orfdata viral completeness --incomplete
```

Because the GTF, transcript FASTA and protein FASTA share identifiers, the set
opens as a PyEnsembl genome for downstream variant tools — see
[Downstream tools](downstream.md).

## Cache output

Downloads have a canonical, versioned path under `ORFDATA_CACHE`,
`OPENVAX_DATA_CACHE`, or the normal OpenVax cache. This is why the command does
not require `--output`: an arbitrary export path would also need to define
resource identity, update and verification behavior. Use `--cache-root` to put
the managed cache elsewhere. The command prints the installed paths; add
`--json` for structured output.

## Python acquisition convention

ORFdata follows the OpenVax split between explicit acquisition and offline
lookup:

```python
from orfdata import download_viral_references, viral_reference_bundle

download_viral_references(species="human")
bundle = viral_reference_bundle(species="human")
```

The convenience form is explicit and opt-in:

```python
bundle = viral_reference_bundle(species="human", download_if_missing=True)
```

Context selection is a separate, type-stable API. It returns a
`ViralReferenceIndex`; it uses the base bundle above and fetches it only when
requested:

```python
from orfdata import index_viral_references

cancer_index = index_viral_references(
    subsets=["cancer-causing"],
    species="human",
    download_if_missing=True,
)

integrated_drivers = index_viral_references(
    subsets=["cancer-integrated-driver"],
    download_if_missing=True,
)
```

The default is `download_if_missing=False`, so opening or inspecting a resource
does not unexpectedly use the network. The name matches the lower-level
PyEnsembl/DataCache convention; high-level PyEnsembl and MHCflurry workflows
also expose explicit download or fetch commands.

Across ORFdata download APIs, repeated acquisition reuses a valid installation
and `force=True` refreshes or repairs it, status objects provide `to_dict()`,
human CLI output is the default, and `--json` provides machine-readable output.

A derived index is local: `index_viral_references(force=True)` rebuilds it from
the cached bundle without downloading, and a damaged index is rebuilt rather
than reported, since nothing in it is unrecoverable. Adding
`download_if_missing=True` to that call lets it repair the source bundle too,
which does download.
