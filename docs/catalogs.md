# Optional ORF catalogues

The bundled human catalogue is available through `load_catalog()` without
downloads. `CatalogRegistry` manages additional versioned ORF catalogues, with
the same explicit discovery, download, build and offline-open workflow as
`DatasetRegistry`. Once a species' catalogue in the packaged registry is
downloaded or built, `load_catalog(species=...)` opens it, so `extra_orfs(genome)` and
`Proteome.with_extra_orfs()` use it without `catalog=`; pass
`load_catalog(species=..., registry=...)` as `catalog=` for your own registry.
A bundled catalogue is preferred, a catalogue's default release is opened
unless only another release covers the assembly, and two registered catalogues
for one species and assembly are never chosen between; open one by name
instead. It returns the existing `ORFCatalog`; query results still bind to the
original PyEnsembl transcript through `ReferenceORF`. To use several catalogues
for one species together, see
[merging catalogues](proteome.md#merging-catalogues).

The packaged registry has one entry, the [OpenProt mouse catalogue](#mouse-orfs-from-openprot),
which you build from OpenProt's own files. `orfdata catalogs status` (in Python,
`catalog_status()`) lists every catalogue `load_catalog` can open, the bundled
one included, with its species, assembly, record count, whether it is
installed and which one is selected for each assembly:

```sh
orfdata catalogs status
```

The examples below use a caller-supplied registry.

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

## Mouse ORFs from OpenProt

`openprot-mouse` holds the alternative proteins (AltProts) and novel isoforms
that [OpenProt 2.2](https://openprot.org) reports for mouse with at least two
unique MS peptides, on Ensembl 106 GRCm39 transcripts. orfdata does not
redistribute OpenProt's files, and OpenProt states no licence for them. You
build the catalogue from OpenProt's own downloads (3.4 MB, about 5 seconds);
each file is checked against its pinned size and SHA-256:

```sh
pip install "orfdata[data]"
orfdata catalogs build openprot-mouse
```

After that, `load_catalog(species="mouse")`, `extra_orfs(genome)` and
`Proteome.build(genome)` use it for any GRCm39 Ensembl release:

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

proteome = Proteome.build(EnsemblRelease(115, species="mouse"))
print(proteome.summary()["kinds"], len(proteome.excluded))
```

```text
{'annotated': 66668, 'extra_orf': 2854} 122
```

- **Records.** OpenProt's file has 8,645 rows for 2,913 proteins; each row
  places one protein on one transcript. The 4,074 rows on Ensembl transcripts
  become 2,976 records, one per protein and spliced coding path, named by the
  OpenProt accession (`IP_` for AltProts, `II_` for isoforms). The 32 proteins
  found at several paths, mostly identical copies at different loci such as
  on chrY, have one record per path (96 in all), named
  `accession@transcript`; they share one evidence summary. A path on several
  transcripts is one record, placed on each of them. `orf_type` is OpenProt's
  protein type.
- **Not imported.** The 4,571 rows on RefSeq transcripts are listed in
  `metadata["unresolved"]`, because no RefSeq annotation is pinned. Every
  protein still has a record: the one accession with only RefSeq rows,
  IP_957586, has the same sequence and path as II_957586, which is imported.
- **Genes.** OpenProt gives Ensembl gene IDs, not names, so query the
  catalogue by gene ID (`--gene ENSMUSG00000020149`); proteome entries and
  placements take gene names from the PyEnsembl annotation. A split
  accession's bare ID, such as `IP_824771`, names all of its records.
- **Coordinates.** Coding blocks come from OpenProt's BED file. Its README
  calls the TSV's genomic coordinates 1-based, but they are one lower than the
  BED's and Ensembl's at both ends
  ([#55](https://github.com/openvax/orfdata/issues/55)). The build checks every
  row against that offset, its transcript span and its protein length, and
  fails rather than guessing if any disagrees.
- **Checked on Ensembl.** On Ensembl 106, OpenProt's own release, 2,864 records
  bind to their source transcript with their exact coding path and protein.
  The other 112 lie on transcripts that are in Ensembl's GTF but missing from
  its cDNA and ncRNA files, so there is no sequence to check. On Ensembl 115,
  2,854 place; of the rest, 109 lack a sequence, 10 transcripts are gone and
  3 splice paths changed. Each is listed with its reason.
- **Evidence.** OpenProt's MS score (a sum of unique peptides per study) and
  translation-event score (a count of studies) become aggregate summaries for
  the protein, not peptides or named studies. So `MS()` selects every record
  and `RiboSeq()` the 137 placed on Ensembl 115 with translation evidence,
  but `MS(min_peptides=2)` and `RiboSeq(min_studies=2)` select none. A sum of
  per-study peptide counts is not a count of distinct peptides, and a count
  of studies names no studies, so neither can show two independent
  observations.
- **Other assemblies.** GRCm38 and NCBIM37 genomes use the catalogue lifted
  automatically once the chains are downloaded (UCSC's GRCm39-to-GRCm38 chain
  needs `--accept-terms`); see [every assembly of a species](#every-assembly-of-a-species).
- **Dog.** OpenProt has no dog data, and orfdata has no dog catalogue. Mouse
  evidence is never transferred to dog through orthology.

Cite Leblanc S, et al. *OpenProt 2.0 builds a path to the functional
characterization of alternative proteins.* Nucleic Acids Res. 2024;52(D1):D522–D528,
[doi:10.1093/nar/gkad1050](https://doi.org/10.1093/nar/gkad1050).

## Every assembly of a species

Catalogues on another assembly of the species are lifted automatically. Download
the pinned chain files once, and `load_catalog`, `extra_orfs(genome)` and
`Proteome.build(genome)` work on every human and mouse assembly PyEnsembl
knows:

```sh
orfdata chains download                  # Ensembl's chains; add --accept-terms for UCSC's mouse ones
```

```python
from pyensembl import EnsemblRelease
from orfdata import Proteome

Proteome.build(EnsemblRelease(75))       # GRCh37: the GRCh38 catalogue, lifted
```

With evidence `None` (every catalogue ORF), each assembly gets:

| Genome | ORF source | Records | Placed | Not placed: main reasons |
| --- | --- | ---: | ---: | --- |
| Human Ensembl 115, GRCh38 | bundled | 7,264 | 7,239 | 24 splice paths changed |
| Human Ensembl 75, GRCh37 | lifted, 1 hop | 7,247 | 6,857 | 271 transcripts absent, 108 paths changed |
| Human Ensembl 54, NCBI36 | lifted, 1 hop | 7,253 | 4,195 | 2,710 transcripts absent |
| Mouse Ensembl 115, GRCm39 | OpenProt | 2,976 | 2,854 | 109 without sequence |
| Mouse Ensembl 102, GRCm38 | lifted, 1 hop | 2,975 | 2,856 | 112 without sequence |
| Mouse Ensembl 67, NCBIM37 | lifted, 2 hops | 2,922 | 1,620 | 1,187 transcripts absent |

Older annotations lack many of the transcripts the sources name, so fewer ORFs
place; the rest are listed with the reason, never silently dropped.

- **Which catalogues.** Each catalogue contributes one release: the one on the
  requested assembly if it has one, else its default, lifted. Each assembly
  contributes one catalogue (a bundled one if any), and catalogues from
  several assemblies are merged. This depends only on the registry and the
  pinned chains, never on what happens to be installed. A missing catalogue
  or chain raises an error naming the command that installs it; nothing is
  downloaded while loading.
- **Routes.** `chain_routes(species, build)` gives the fewest-hop route to
  each assembly the pinned chains reach, for example GRCm39 to NCBIM37
  through GRCm38. `lift_to(catalog, build)` lifts through it.
  `catalog_status()` lists, for each catalogue, the assemblies it reaches and
  whether their chains are downloaded.
- **Cost.** Lifting happens on each load: about 1.6 s for the human catalogue.

## Lifting a catalogue by hand

`orfdata lift` writes a lifted catalogue to a file, for any route the pinned
chains cover or with your own chain file:

```sh
orfdata lift --species human --to GRCh37 --output orfs-grch37.json.gz
```

```text
  "lifted": "orfdata-catalog:2026-10-11.1@GRCh37",
  ...
  "lifted_records": 7247,
  "unmapped": 17,
  "unmapped_reasons": {
    "A coding block crosses an alignment gap or edge": 2,
    "A coding block is not aligned": 15
  },
```

`orfdata lift` writes the lifted catalogue and prints a summary to standard
error. `--catalog FILE` lifts a local catalogue instead, and `--chain FILE` uses
your own chain file instead of a downloaded pin. Pass the result back with
`--catalog`, for example to `orfdata proteome`; JSON and TSV exports include
each record's original and the chain. Downloading chains needs the
`orfdata[data]` extra, which adds DataCache. In Python:

```python
from orfdata import download_chain, lift_catalog, load_catalog, load_chain

download_chain("GRCh38", "GRCh37")  # Once; checksum-verified. download_chains() gets them all.
grch37 = lift_catalog(load_catalog(), load_chain("GRCh38", "GRCh37"))
grch37.save("orfs-grch37.json.gz")
print(len(grch37.records), len(grch37.metadata["liftover"]["unmapped"]), grch37.source)
```

```text
7247 17 orfdata-catalog:2026-10-11.1@GRCh37
```

- **Pinned chain files.** Ensembl's human GRCh38 to and from
  [GRCh37 and NCBI36](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/),
  and mouse GRCm38 to and from
  [NCBIM37](https://ftp.ensembl.org/pub/assembly_mapping/mus_musculus/),
  have no usage restrictions. `orfdata chains list` shows all eight. UCSC's mouse
  [`mm10ToMm39`](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/README.txt)
  and [`mm39ToMm10`](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/README.txt)
  (GRCm38 and GRCm39) are free only for non-commercial use, so downloading
  them needs `--accept-terms` (`accept_terms=True`). `--source-path FILE`
  installs a copy you already have after checking it against the pin.
  Ensembl publishes no GRCm38/GRCm39 chain, and no chain is pinned for dog.
- **Other chain files.** `ChainFile(path, source_build=..., target_build=...)`
  reads any UCSC or Ensembl chain file, gzipped or not. A chain file doesn't
  name its assemblies, so you declare them. When the file is one of the
  pinned ones, its direction and species are checked, so a GRCh37-to-GRCh38
  file can't lift a GRCh38 catalogue; any other file is taken at your word.
  UCSC's chains to GCF
  assemblies, such as `canFam3ToGCF_014441545.1.over.chain.gz` for dog
  ROS_Cfam_1.0, name contigs by RefSeq accession, so pass `contig_names=` built
  from NCBI's assembly report to match the genome's names. Check each file's
  terms before using or sharing what it produces.
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
  and the chain's file name and SHA-256 in `record.liftover`, the same whether
  the file was downloaded as a pin or passed by hand. Its ID, protein and
  evidence must equal the original's. Its initiation flanks move with it when
  they map consistently and in order. A strain pin is kept, with the chain's
  SHA-256 in place of an annotation. `metadata["liftover"]["input"]` describes
  the input catalogue, including any earlier lift. `ORFCatalog.save` writes
  the same bytes for the same catalogue.
- **Placement.** `extra_orfs` on the target genome then checks each ORF's
  exact coding path and protein, so a lift that changes the sequence is
  excluded rather than trusted. A lifted record names no annotation release.
  ORFs that didn't lift are reported among the exclusions, as "Did not lift
  to GRCh37: …", when the evidence rule selects them; through every hop of
  an earlier lift and every merged input. On Ensembl 75, `"ms"` selects all
  2,186 MS-supported ORFs: 2,097 place, and 3 of the 89 excluded didn't lift.
- **Size.** Each lifted record carries its original, so a lifted catalogue is
  about twice the size of its input (33 MB for the bundled one as JSON).
- **Errors.** When no pinned chains reach an assembly, such as T2T-CHM13,
  `load_catalog(species=..., genome_build=...)` gives the `orfdata lift
  --chain` command and the Python call that lifts with your own chain file.

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
