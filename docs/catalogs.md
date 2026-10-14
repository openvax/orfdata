# Optional ORF catalogues

The bundled human catalogue is available through `load_catalog()` without
downloads. `CatalogRegistry` manages additional versioned ORF catalogues, with
the same explicit discovery, download, build and offline-open workflow as
`DatasetRegistry`. Once a species' catalogue in the packaged registry is
downloaded, `load_catalog(species=...)` opens it, so `extra_orfs(genome)` and
`Proteome.with_extra_orfs()` use it without `catalog=`; pass
`load_catalog(species=..., registry=...)` as `catalog=` for your own registry.
A bundled catalogue is preferred, a catalogue's default release is opened
unless only another release covers the assembly, and two registered catalogues
for one species and assembly are never chosen between; open one by name
instead. It returns the existing `ORFCatalog`; query results still bind to the
original PyEnsembl transcript through `ReferenceORF`.

The default optional registry currently has no released entries. The examples
below use a caller-supplied registry; they do not imply a published
experimental catalogue for another species.

```python
import json
from orfdata import CatalogRegistry

with open("catalog-registry.json") as handle:
    registry = CatalogRegistry(json.load(handle))
releases = registry.list(species="yeast")
for manifest in releases:
    print(manifest.catalog_id, manifest.version, manifest.coverage)

# Select an entry from your populated registry.
name = releases[0].catalog_id
registry.download(name)             # Explicit, checksum-verified acquisition.
catalog = registry.open(name)       # Offline; no implicit download or repair.
records = catalog.query(evidence="ms", require_resolved=True)
orfs = catalog.for_transcript(transcript)  # An already loaded PyEnsembl transcript.
```

Use `info(name)` for immutable provenance and `download_metadata(name)` for
artifact URLs, checksums and access terms. `inspect(name)` checks the local
installation. Species, strain, assembly and annotation contexts remain attached
to accepted records; a matching accession alone cannot transfer evidence.

```sh
orfdata catalogs --registry catalog-registry.json list --species yeast
orfdata catalogs --registry catalog-registry.json download NAME
orfdata catalogs --registry catalog-registry.json query NAME --evidence ms --require-resolved
```

Queries share the ordinary ORF CLI's gene, transcript, genomic-position and
evidence filters, and JSON/TSV/FASTA exports. They cannot override the selected
release with a local `--catalog` argument. Use `orfdata --catalog FILE` for the
existing standalone local-file workflow.

## Source builds

`registry.build(name, source_paths={"source.json": path})` explicitly rebuilds
an allowlisted source adapter. `registry.open_build(name)` opens that local
result, and `inspect_build(name)` checks its input pins, transform identity and
output receipt without reopening original input files. The CLI equivalents are
`catalogs build NAME --source-path source.json=PATH` and
`catalogs query NAME --built`. A local receipt's checksum is distinct from an
independently pinned published artifact; `inspect_build().verified` remains false.

Builds are atomic. `force=True` creates a fresh bundle while existing readers'
files survive. Cache roots follow `ORFDATA_CACHE`, `OPENVAX_DATA_CACHE`, or the
explicit `cache_root`; catalogues, quantitative datasets and local source builds
use separate locations.

## Defining a release

A registry definition maps a name to `default_version` and `versions`. Each
version contains a `CatalogManifest`, `distribution` access terms and an
`assets` entry for raw `catalog.json.gz` bytes, pinned by URL, SHA-256 and size.
The manifest declares references, original `SourceAsset` inputs, citations,
processing and coverage. `coverage.records` is required; source importers also
report original assignment counts and unresolved categories. Distinct protein
sequences and evidence summaries are separate counts.

Managed files use the existing `orfdata.catalog.v1` JSON representation and embed
the full manifest in `metadata.catalog_manifest`. Every accepted record must
declare a complete genomic reference included in that manifest. Opening validates
the whole manifest, species/assembly, references and declared output counts.

The initial `catalog-copy-v1` adapter wraps an independently prepared,
unmanaged, pinned portable catalogue named by `processing.catalog_asset`; it
preserves native metadata and adds the managed manifest. It rejects an input
that already has a managed manifest instead of replacing its release provenance.
It does not infer missing references or validate new biological models on
behalf of a source importer.

Python developers can use `build_catalog(manifest, builder)`, where
`builder(inputs, output_path)` writes the managed file. DataCache verifies inputs and validates the completed
output before publishing the local build. Importer paths in metadata are never
executed as Python code.

[OpenProt source qualification](https://github.com/openvax/orfdata/issues/18)
and its [coordinate discrepancy](https://github.com/openvax/orfdata/issues/55)
remain separate from this data-management API.
