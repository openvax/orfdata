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
original PyEnsembl transcript through `ReferenceORF`. To use several catalogues
for one species together, see
[merging catalogues](proteome.md#merging-catalogues).

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

## Lifting a catalogue to another assembly

A catalogue built on one assembly can be lifted to another of the same species
with a chain file, for example GRCh37 to GRCh38 or mm10 to GRCm39:

```python
from orfdata import ChainFile, extra_orfs, lift_catalog, load_catalog

chain = ChainFile("GRCh38_to_GRCh37.chain.gz", source_build="GRCh38", target_build="GRCh37")
grch37 = lift_catalog(load_catalog(), chain)
print(len(grch37.records), len(grch37.metadata["liftover"]["unmapped"]), grch37.source)
```

```text
7247 17 orfdata-catalog:2026-10-11.1@GRCh37
```

- **Chain files.** `ChainFile` reads UCSC and Ensembl chain files, gzipped or
  not. Examples are Ensembl's
  [`GRCh37_to_GRCh38.chain.gz`](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/)
  and UCSC's [`mm10ToMm39.over.chain.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/).
  UCSC's chains to GCF assemblies, such as
  `canFam3ToGCF_014441545.1.over.chain.gz` for dog ROS_Cfam_1.0, name contigs
  by RefSeq accession, so pass `contig_names=` built from NCBI's assembly
  report to match the genome's names. Download chain files yourself and check
  their terms; orfdata never fetches one.
- **What lifts.** A record lifts only when every coding block lies inside one
  ungapped alignment block and all blocks land on one contig and strand, in
  order. Where a region also aligns to fix patches or alternate haplotypes,
  the same-named chromosome is used. The 17 bundled ORFs above that don't
  lift are listed in `metadata["liftover"]["unmapped"]` with the reason (and,
  when a block aligns to several places, their `contigs`): 15 have no
  alignment and 2 cross an alignment gap. UCSC's names for unplaced scaffolds,
  such as `chr1_GL456210v1_random`, are read as Ensembl's (`GL456210.1`).
- **Real data.** Lifting the 7,247 back to GRCh38 with
  `GRCh37_to_GRCh38.chain.gz` recovers the original coordinates of the 7,246
  that lift back. On Ensembl 75, 6,857 of the lifted ORFs place on a
  transcript with their exact coding path and protein.
- **Provenance.** Each lifted record keeps its original in `record.lifted_from`
  and the chain's name and SHA-256 in `record.liftover`. Its ID, protein and
  evidence must equal the original's. Its initiation flanks move with it when
  they map consistently and in order. A strain pin is kept, with the chain's
  SHA-256 in place of an annotation. `metadata["liftover"]["input"]` describes
  the input catalogue, including any earlier lift.
- **Placement.** `extra_orfs` on the target genome then checks each ORF's
  exact coding path and protein, so a lift that changes the sequence is
  excluded rather than trusted. A lifted record names no annotation release.
  ORFs that didn't lift are reported among the exclusions, as "Did not lift
  to GRCh37: …", when the evidence rule selects them; through every hop of
  an earlier lift and every merged input. On Ensembl 75, `"ms"` selects all
  2,186 MS-supported ORFs: 2,097 place, and 3 of the 89 excluded didn't lift.
- **Size.** Each lifted record carries its original, so a lifted catalogue is
  about twice the size of its input (33 MB for the bundled one as JSON).
- **Errors.** When only another assembly has a catalogue,
  `load_catalog(species=..., genome_build=...)` says which one to lift.

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
