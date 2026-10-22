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
command. Focused command groups such as `orfdata viral coverage` remain useful
for domain-specific inspection, and `orfdata viral download` remains a
compatibility alias.

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
`missing`, `corrupt`, or `inaccessible`.
`downloadable_data()` and `Downloadable` remain aliases for early adopters of
this API.

## Viral collections and contexts

The published viral collection currently has these properties:

| Property | Value | Meaning |
| --- | --- | --- |
| Species | `human` | Target analysis species, NCBI taxon 9606 |
| Collection | `gdc-grch38-d1-vd1` | Stable identity of the 200-reference GDC set |
| Contexts | `cancer-transcriptome-screening` (`cancer`), `cancer-associated`, `cancer-causing` | Full screen plus IARC evidence tiers |

Run the default collection with:

```sh
orfdata download viral --species human
```

`--collection gdc-grch38-d1-vd1` is available for explicit configuration and
for scripts that will later choose among several collections. Context is
metadata rather than a sequence identity or a mutually exclusive selector. CMV and EBV,
for example, may be requested in both cancer and autoimmune workflows. A future
collection may carry several context tags, and the same exact accession/version
can belong to several collections without becoming several biological
references.

Pass `--context` to materialize a derived FASTA index from the shared source
bundle. Repeated contexts select their union; `--context-match all` selects
their intersection. It reuses source records already downloaded for the
collection, so overlapping contexts do not create conflicting copies of a
reference.

```sh
orfdata viral contexts
orfdata download viral --context cancer  # all 200 GDC screening references
orfdata download viral --context cancer-associated  # 37 IARC Group 1, 2A or 2B references
orfdata download viral --context cancer-causing  # `iarc-group-1` is an alias
# The other 180 GDC screening references are not described as non-carcinogenic:
orfdata download viral --exclude-context cancer-causing
# Once another listed context is published, select its union or intersection:
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context>
orfdata download viral --context cancer-transcriptome-screening --context <another-listed-context> --context-match all
```

The target species is an admission and organization choice for the reference
set. It does not replace a source record's viral taxon or assert that the source
declared a host. Unknown source-host metadata remains unknown.

The four facts stay separate in the downloaded manifest and FASTA headers:
the versioned viral accession and CDS feature identify the source product;
collection and context record where it was selected; the target species records
the analysis it supports; and a source-host taxon is populated only when the
source establishes one. This keeps the same reference idempotent when it is
obtained through several collections or exporter paths.

The word *panel* remains in `viral_panel_coverage()` and the GDC coverage page
because it describes the source inventory. The managed download API calls the
versioned artifact a collection, which is the value used for lookup and cache
identity.

Use `cancer` for the full 200-reference GDC screen. Use `cancer-associated`
for the 37 references in IARC Groups 1, 2A or 2B, and `cancer-causing` for
the 20 [Group 1 references](https://monographs.iarc.who.int/agents-classified-by-the-iarc/).
These selections represent evidence tiers: references outside one tier are not
evidence of no cancer association. Use `--exclude-context` only when that
operational selection is what a workflow needs. Respiratory and autoimmune
contexts can be added from their own versioned, scoped authorities and can
overlap this one.

Mechanism selectors are separate from contexts because they describe
virus--tumour biology, not a collection's intended analysis use. Named subsets
provide reusable selections composed from mechanisms. List them with `orfdata
viral subsets`. For a conservative set with evidence of both host-genome
integration and direct tumour-cell driving, use:

```sh
orfdata download viral --subset integrated-direct-driver
```

This derived index contains 15 references: the 12 IARC Group 1 HPV types,
HBV, HTLV-1 and MCV. It must not be read as a treatment-response prediction;
the relevant tumour, integration event, viral expression and clonal state need
sample-level evidence.

For the narrower evidence set where continued viral-oncogene expression has
been tested in named tumour settings, use:

```sh
orfdata download viral --subset tumor-cell-oncogene-dependence
```

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
    contexts=["cancer-causing"],
    species="human",
    download_if_missing=True,
)

integrated_drivers = index_viral_references(
    subsets=["integrated-direct-driver"],
    download_if_missing=True,
)
```

The default is `download_if_missing=False`, so opening or inspecting a resource
does not unexpectedly use the network. The name matches the lower-level
PyEnsembl/DataCache convention; high-level PyEnsembl and MHCflurry workflows
also expose explicit download or fetch commands.

Across ORFdata download APIs, repeated acquisition reuses a valid installation,
`force=True` refreshes or repairs it, status objects provide `to_dict()`, human
CLI output is the default, and `--json` provides machine-readable output.
