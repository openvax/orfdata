# Discover and download data

Run `orfdata` for the command guide, or `orfdata data` to see its subcommands.
ORFdata has one common entry point for data acquisition:

```sh
orfdata data list
orfdata data list --species dog
orfdata data info pfizer-atlas-dog-peptide-evidence
orfdata data fetch pfizer-atlas-dog-peptide-evidence
orfdata data fetch peptideatlas-dog-peptides openprot-mouse
orfdata data path pfizer-atlas-dog-peptide-evidence
orfdata data status
```

`data`, `downloads` and `download` are aliases for the same command group.
The `list/info/fetch/path` vocabulary follows other OpenVax download tools.
`fetch` accepts one or more explicit resource names and chooses the registered
acquisition route: download an artifact when available, otherwise build locally
from checksum-pinned source files. All names and required terms are checked
before the first acquisition. There is no implicit download-all operation.

`info NAME` reports the version, source/coverage manifest and distribution
metadata, including sizes and redistribution terms where declared. `path NAME`
prints an installed file or directory for scripts; missing data is an error,
and never triggers a download. Use `--verify` on reads to validate local bytes,
`--cache-root DIR` to select a cache, and `--json` for machine output.
`fetch --json` returns one result per distinct requested name. Human tables
preview 25 rows; use `--limit 0` for all rows. `--json` always returns complete
results. See the [CLI guide](cli.md) for output formats and migration notes.

Use `--kind` to disambiguate colliding resource names. For one resource,
`--version VERSION` selects a historical release and `--source-path ASSET=PATH`
provides already downloaded catalogue/dataset sources or genome files. Source
overrides remain checksum-checked. `fetch` reuses a valid installation and
replaces one whose files fail verification, reporting it as `repaired`, as
`setup` does; `--force` installs a fresh generation even over valid files, and is
needed to replace files that can't be read.
Restricted chain downloads require `--accept-terms` after reviewing `info`.

The same operations are available in Python:

```python
from orfdata import data_info, data_path, fetch_data

details = data_info("openprot-mouse")  # Offline; shows the build route and terms.
result = fetch_data("openprot-mouse")  # Explicit source download and processing.
path = data_path("openprot-mouse")    # Offline file/directory lookup.
```

The focused forms remain available for additional resource-specific options:

```sh
orfdata download viral --species human
orfdata download catalog NAME
orfdata download dataset NAME
orfdata download reference DATASET
orfdata download genome NAME
orfdata download chains --species human
```

The operation comes first so the same discovery and installation model applies
to every resource kind. Focused command groups such as `orfdata viral coverage` and
`orfdata viral summary` remain useful for domain-specific inspection.

`list` shows every bundled, downloadable or locally buildable primary resource
and its current state, including native protein search references declared by
quantitative datasets. `status` shows resources with local files: usable ones
and damaged or unreadable installs. `data fetch` replaces damaged files;
unreadable ones are left alone until their permissions are fixed or `--force` is
given. For catalogues, datasets, search references, genomes and chains, a repair
installs a fresh generation and leaves the damaged one in place until the cache
is pruned, as `--force` does; a viral collection is rebuilt in place. Filter either command with `--kind` or `--species`; use `--installed` or
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
from orfdata import DataStatus, list_data

all_resources: tuple[DataStatus, ...] = list_data()
human_missing: tuple[DataStatus, ...] = list_data(species="human", installed=False)
installed: tuple[DataStatus, ...] = list_data(installed=True)
machine_record: dict[str, object] = all_resources[0].to_dict()
```

Each `DataStatus` reports its kind, stable name and version, species,
nonexclusive analysis contexts, local status, acquisition method and equivalent
command, and resolved cache path. Status values follow DataCache: `available`,
`missing`, `corrupt`, or `inaccessible`. An unknown kind, species or context is
an error rather than an empty listing, so a typo never reads as "nothing is
available".

## Interrupted downloads

Downloads restart from the beginning by default. On an unreliable connection,
add `--resume` to `orfdata data fetch` or `orfdata setup` (`resume=True` in
`fetch_data` and `install_species`, or `resume="auto"` in a registry's download
options). Rerunning the same command then requests only the bytes not yet
received for published catalogues, genome and chain files, search references,
the source files of a local build, and split archives such as the four parts
of `pfizer-atlas-dog-peptide-evidence`.

Resuming costs disk space, which is why it is off by default:

- Partial downloads and their staging stay in the cache after a failure, in
  owner-only `.staging-*`, `.inputs-*` and `.datacache-resume-*` directories
  beside the resource. A successful install of the same resource removes them.
  Delete them yourself to reclaim the space from a download you abandon.
- An archive keeps every downloaded part until the dataset is installed. For
  `pfizer-atlas-dog-peptide-evidence` that adds the 6.9 GB of parts to the peak
  space the install needs, which is otherwise about 32 GB.
- An archive's staging directory isn't named per user, so in a shared cache
  one user's interrupted download stops another's until it is removed.

Datasets published as one compressed SQLite file can't resume: that is every
published dataset except `pfizer-atlas-dog-peptide-evidence`, including the
1.27 GB `pfizer-atlas-mouse-rna`. DataCache can't resume a transfer it
decompresses while downloading, so with `--resume` these still restart, and
publishing them in a resumable form is tracked in
[#156](https://github.com/openvax/orfdata/issues/156). Sources that `setup`
shares between catalogues are staged for one run and restart too. A viral
collection is fetched record by record and doesn't resume.

A resumed file is checked against the same pinned SHA-256 and size as a
complete one. If the finished file doesn't match, the install fails; rerun
without `--resume`, or delete the partial, to download it again in full.

## Viral collections and contexts

The Python API uses the same collection object for discovery, inventory reads,
source summaries, downloads and offline opens. The target species comes from the
collection; pass `species` only when you want an explicit consistency check.

```python
from pandas import DataFrame

from orfdata import ViralCollection, viral_collection, viral_collections

human_collections: tuple[ViralCollection, ...] = viral_collections(species="human")
collection: ViralCollection = viral_collection("gdc-viruses")
assert viral_collection(collection.panel) is collection

inventory: DataFrame = collection.inventory()
genbank: DataFrame = collection.source_summary()       # primary source
pave: DataFrame = collection.source_summary("pave")    # another declared source
```

Names, panel names and aliases resolve through `viral_collection()` everywhere.
An unknown species, collection or source is an error rather than an empty
result. A second species is another `ViralCollection` declaration with its own
checked inventory and source summaries; download code does not branch on the species.

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

Pass `--context` to materialize a derived FASTA set from the shared source
collection. Repeated contexts select their union; `--context-match all` selects
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
from orfdata import open_viral_references

references = open_viral_references(species="human")
for layer in references.layers:
    print(layer.name, layer.molecule, layer.records, layer.fasta)
transcriptome = references.layer("transcriptome")
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
from orfdata import ViralReferenceFiles, download_viral_references, open_viral_references

download_viral_references(collection="gdc-viruses")
references: ViralReferenceFiles = open_viral_references(collection="gdc-viruses")
```

The convenience form is explicit and opt-in:

```python
references = open_viral_references(collection="gdc-viruses", download_if_missing=True)
```

Context selection is a separate, type-stable API. It returns paths and counts
for a derived reference set; it uses the installed collection above and fetches it only
when requested:

```python
from orfdata import select_viral_references

cancer = select_viral_references(
    subsets=["cancer-causing"],
    species="human",
    download_if_missing=True,
)

integrated_drivers = select_viral_references(
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

A derived selection is local: `select_viral_references(force=True)` rebuilds it from
the cached collection without downloading, and a damaged selection is rebuilt rather
than reported, since nothing in it is unrecoverable. Adding
`download_if_missing=True` to that call lets it repair the source collection too,
which does download.
