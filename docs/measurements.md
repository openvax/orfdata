# RNA, ribosome and protein measurements

The quantitative layer of orfdata stores measurements of genes, transcripts,
reference ORFs, proteins, protein groups, peptides and source-defined regions.
A measurement stays at the resolution the experiment supports. Parent-gene RNA
abundance is useful context for an ORF; it is not a separate measurement of each
ORF in that gene.

Released datasets, each with its own guide:

- [Human, mouse and dog protein abundance](pfizer-protein.md) (Pfizer/Oziolor atlas)
- [Dog and mouse RNA abundance](pfizer-rna.md) (Pfizer/Oziolor atlas)
- [Matched human tissue RNA and protein](human-tissue.md) (Wang et al.)
- [Yeast protein and RNA abundance](copy-yeast.md) (CoPY)
- [Yeast RNA and ribosome occupancy](ingolia-yeast.md) (Ingolia et al.)
- [Dog MHC ligands and PeptideAtlas peptides](dog-peptides.md) (Kaabinejadian et al.; Kusebauch et al.)
- [Pfizer dog peptide and protein-group evidence](pfizer-peptides.md) (Oziolor et al.)

See [released coverage and gaps](coverage.md); the
[data roadmap](https://github.com/openvax/orfdata/issues/3) tracks further
imports. The bundled human ORF evidence catalogue remains available offline.
Its pooled PSM counts describe spectra, not sample-level protein abundance.

## Discover, download and query

Install `orfdata` to manage datasets. Creating a registry and listing
releases reads metadata only:

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()  # Metadata only; no download.
releases = registry.list(species="dog", assay="digest_ms")
```

For a registered release:

| Call | What it does |
| --- | --- |
| `registry.info(name, version)` | Returns its immutable manifest. |
| `registry.download(name, version)` | Fetches and verifies it. |
| `registry.open(name, version)` | Reads it offline. |
| `registry.inspect(name, version)` | Reports missing, available or invalid files without downloading or repairing them. |
| `registry.download_metadata(name, version)` | Returns the artifact's access terms and file hashes and sizes. |

An invalid cache requires an explicit `download(..., force=True)`; old reader
paths survive. Python download options such as timeout, retries and progress go
to DataCache.

Query an opened dataset by an exact, reference-scoped target. Here `transcript`
is a PyEnsembl Transcript in the dataset's exact reference context:

```python
dataset = registry.open(name, version)
target = ScopedTarget.from_pyensembl(
    transcript, kind="transcript", namespace="ensembl:110"
)
for measurement in dataset.measurements(target=target):
    print(measurement.value, measurement.unit, measurement.status)
```

Use `ScopedTarget.from_orf(bound_reference_orf, namespace=source_release)` for a
validated reference ORF. Namespaces and identifiers must match the source's
mapping table. The default PyEnsembl target uses its unversioned `id`; request
`versioned=True` when the source uses versioned identifiers. There is no implicit
mapping across annotation releases, gene isoforms, ORFs or species. A variant or
fusion product does not inherit a native reference protein measurement.

Use `ScopedTarget(reference, "region", source_id, source_namespace)` for a
quantified UTR, intron or other source-defined window. Keep its native feature
type and target length in measurement provenance. A region ID does not assert
known genomic coordinates or a translated ORF. RNA and Ribo measurements of the
same region can be compared through an explicit `same_target` relationship;
such a density ratio describes that region, not a particular uORF inside it.
Regions do not automatically bind to PyEnsembl or participate in gene/ORF orthology.

When a source also defines a candidate or partial coding model, connect the two
with `TargetRelation(..., relation="region_orf_context", ...)`: the left target
is the measured `region`; the right is an `orf` target created from its
`TranscriptORF` with `ScopedTarget.from_orf`. Annotated and source-defined ORFs
work without a catalogue record. Supply the source's mapping evidence, such as
the reported interval relationship, coordinates and source record IDs, and its
cardinality (`one_to_many` or `unknown` when appropriate). Keep each candidate
as its own relation. This preserves a region's relationship to coding models
without converting region measurements into ORF abundance or translation
evidence. `MatchedMeasurements.ratio` rejects this context relation even if the
source reports only one candidate. Partial boundary/frame details remain on the
referenced ORF; the relation does not fill them in.

For native sources with unreported genomic metadata, use explicit
[source references and mapping views](source-references.md) to retain the original
observation alongside its mapping to a PyEnsembl target.

Choose a native series explicitly when a source supplies several units:

```python
values = dataset.measurements(quantity="protein_abundance", unit="iBAQ", scale="linear")
dataset.export_parquet("ibaq.parquet", unit="iBAQ", scale="linear")
```

`unit`, `scale` and `log_base` also work with `mapped_measurements` and existing
storage v1-v4 releases. Units are exact and case-sensitive. Selecting
`scale="log", log_base=2` returns values already on that scale; nothing is transformed.
An unavailable series returns no observations. Contradictory linear/log requests
raise an error. Omitting selectors preserves all native series.

For example, Pfizer iBAQ and rIBAQ are two series from the same source rows.
Selecting or exporting both does not create additional donors, spectra or
independent detections. Source row/sample/assay lineage stays attached to each
measurement.

The CLI uses the same operations:

```sh
orfdata datasets list --species yeast
orfdata datasets info DATASET --version RELEASE
orfdata datasets download DATASET --version RELEASE
orfdata datasets inspect DATASET --version RELEASE
orfdata datasets query DATASET --version RELEASE --quantity protein_abundance
```

Queries show a human preview by default. Use `--format jsonl` to stream original
JSON lines, or `--json` for one complete JSON array. `--reference ref:SHA256` selects an exact context;
target queries also require `--target`, `--kind` and `--namespace`. Unqualified
sample/assay ID queries across multiple references are rejected. ORF catalogue queries use `orfdata query --gene TPST1`; add `--json` for the
complete structured export. See the [CLI migration guide](cli.md#migrating-scripts-from-042).

## What each value means

`Measurement` retains the dataset/version, reference-scoped target, sample and
assay IDs, value, native unit, linear/log scale and log base, uncertainty and its
definition, detection limit, status, origin and provenance. Measured zero,
below detection, not measured, not mapped and missing are distinct. Censored or
missing values use `None`; native sentinels can remain in provenance.

RNA counts/TPM/FPKM, ribosome occupancy, source-defined translation efficiency,
LFQ/iBAQ/TMT values, peptide intensity, PSM counts, copies per cell and turnover
remain distinct. No conversion is automatic. Derived/predicted values name a
method and retain source/input provenance. An estimate from calibration remains
an estimate; an intensity/TPM ratio does not acquire molecule-ratio units.

`Sample` records biological units, pools and source labels. `Assay` separately
records technical replicates, runs, fractions, TMT plex/channel/bridge,
normalization and calibration. A sample's strain/breed belongs in its attributes;
it need not be the strain used to assemble its reference genome. `TargetMapping`
retains all candidate assignments, an unmapped reason, or known candidates with
an explicit `unresolved` reason when source endpoints are missing. A protein-group
representative does not become a uniquely measured protein. See
[pairing and orthology comparisons](pairing.md) for explicit experimental
relationships and derived RNA/protein ratios.

The public query shapes are shared across sources and species:
`samples()`, `assays()`, `measurements()`, `mappings()`, `relations()` and
`reference_mappings()`. Source-specific fields stay in evidence and provenance;
native-table inspection belongs in the source modules. For example:

```python
assignments = dataset.mappings(source_identifier=peptide_sequence)
for assignment in assignments:
    print(assignment.status, assignment.targets, assignment.method)
```

`mappings()` also accepts an exact `source_namespace`, `reference` or destination
`target`. Selecting one destination preserves all candidates in each returned
mapping. Reported assignments and computed sequence compatibility remain
separate records with explicit methods and evidence. Mapping status describes
the assignment's cardinality or completeness, not independent protein detection.

## Importing a source

Define a `DatasetManifest` with exact reference contexts, raw source URLs,
SHA-256 hashes and byte sizes, citations, file-specific access/redistribution
terms, searched/measurable coverage, processing software/version/parameters and
row accounting. Unavailable metadata stay unknown. Article licensing does not
substitute for the data files' terms.

An importer calls `QuantitativeDataset.create(path, manifest, samples=...,
assays=..., measurements=..., mappings=...)`. These inputs may be generators;
measurements are written incrementally.

Every observation names a source asset or earlier input measurements. Derived
inputs precede outputs in the same dataset release. A derived measurement names
its inputs in one of two forms, which cannot be mixed:

- `input_measurements`: IDs within the output reference.
- `input_references`: `MeasurementReference` records with each input's dataset,
  version, reference and ID.

Cross-species derived abundance is rejected. Duplicate IDs, broken foreign keys,
wrong dataset releases and cross-reference links fail the import. Failed builds
do not replace data.

Use `build_dataset(manifest, builder, ...)` for a managed local import.
`builder(source_paths, output_path)` reads the verified inputs and creates that
SQLite dataset. DataCache publishes it atomically with input/output checksums
and a processing receipt. `source_paths={asset_id: local_file}` supports manually
acquired sources. Importing locally does not grant redistribution rights.

A `DatasetRegistry(definitions, cache_root=...)` registers releases as:

```python
definitions = {
    name: {
        "default_version": version,
        "versions": {version: {
            "manifest": manifest.to_dict(),
            "distribution": {"access_terms": terms, "redistributable": False},
            "assets": {"dataset.sqlite": {
                "url": url, "sha256": sha256, "size": size, "decompress": False,
            }},
        }},
    }
}
```

CLI users can supply this JSON with `orfdata datasets --registry FILE ...`.

## Rebuilding a release from source

Rebuild supported releases from their pinned original inputs with
`registry.build(name, version, source_paths={"asset_id": "/local/file"})` or
`orfdata datasets build NAME --source-path ASSET=PATH`. Local overrides still
require the published hashes. Rebuilds use a separate import cache and preserve
installed releases. Unsupported source importers fail explicitly.

Use `registry.open_build(name, version)` to reopen a source rebuild offline, or
`registry.inspect_build(name, version)` to inspect it. Normal `open` still selects
the published artifact. The CLI makes the same choice explicit:

```sh
orfdata datasets build DATASET --version RELEASE
orfdata datasets inspect DATASET --version RELEASE --built
orfdata datasets query DATASET --version RELEASE --built --unit counts
```

`pairs --built` also uses the local rebuild. Read operations never acquire or
repair data.

Build inspection checks source pins, transform identity and output bytes against
the saved receipt, even after local source overrides have been removed. Its
`verified` field remains false because a build receipt is not an independent
published output checksum. Both cache origins retain their original paths and
provenance; explicit force rebuilding preserves existing readers.

## Cache location and storage

The cache root is `ORFDATA_CACHE`, then `OPENVAX_DATA_CACHE`, then the
platform's OpenVax cache. orfdata owns its `orfdata/` child. Release cache
identities include manifest/reference and content hashes.

Datasets are SQLite files. Every table is partitioned by reference key, with
composite foreign keys, so human, mouse, dog, yeast and arbitrary other taxa use
the same schema. [Compact storage and compatibility](storage.md) describes new
SQLite v4 writes, legacy readers and actual-source size/query benchmarks.

The normal installation includes pandas and PyArrow. Use
`dataset.export_parquet(path)` to export query columns alongside the
full measurement JSON, preserving exact values and provenance even when a
floating-point analysis column cannot represent an integer exactly.
