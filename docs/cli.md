# Command-line guide

Run `orfdata` to see the command index and common starting points. It prints
help and exits without reading a catalogue or acquiring data.

```sh
orfdata
orfdata help data
orfdata data list --species dog
orfdata data status
orfdata data info pfizer-atlas-dog-protein
orfdata data fetch pfizer-atlas-dog-protein
orfdata data path pfizer-atlas-dog-protein
```

`list` shows registered resources, including resources that need a local build.
`status` shows resources with files in the selected local cache, including bundled
data and damaged installs that `data fetch` will repair. Neither command
downloads anything. Use `--cache-root DIR` if your data
lives outside the default cache, and `--verify` when you want to hash local files.

## Finding commands

`orfdata help COMMAND` and `orfdata COMMAND --help` are equivalent. Bare command
groups such as `orfdata datasets` show their subcommands. `download` and
`downloads` are aliases for `data`.

| Command | Purpose | Subcommands |
| --- | --- | --- |
| `data` | Discover and acquire resources across all sources | `list`, `status`, `info`, `path`, `fetch`, `catalog`, `dataset`, `reference`, `genome`, `chains`, `viral` |
| `query` | Preview or export ORF records | Use `--gene`, `--transcript`, `--orf-id`, `--species`, or `--catalog` |
| `catalogs` | Manage optional ORF catalogues | `list`, `status`, `info`, `inspect`, `download`, `build`, `query` |
| `datasets` | Manage quantitative measurements and their source references | `list`, `info`, `inspect`, `download`, `build`, `query`, `pairs`, `references`, `download-reference`, `inspect-reference`, `query-reference` |
| `proteome` | Export a reference proteome as FASTA | Use `--species` and `--release`, or `--genome` |
| `viral` | Inspect viral reference coverage and installed layers | `coverage`, `contexts`, `subsets`, `summary`, `completeness`, `inventory download` |
| `setup` | Install a species' default catalogues, chains and collections | Requires `--species` |
| `genomes` | Manage pinned annotations | `list`, `install` |
| `chains` | Manage assembly mapping files | `list`, `download` |
| `lift` | Write an ORF catalogue on another assembly | Requires `--to` and `--output` |

Commands that use a group-level cache or registry option show its placement in
help. For example:

```sh
orfdata datasets --cache-root /path/to/cache query pfizer-atlas-dog-protein
orfdata catalogs --registry registry.json list
```

## Human previews and complete exports

Tables and short field summaries are the default. Tables show up to 25 entries;
`--limit N` changes the preview and `--limit 0` shows every row. Tables wrap to
terminal width. Metadata summaries omit bulky pins, processing details and long
lists with a pointer to the full export.

Use explicit formats in scripts:

```sh
orfdata data list --species dog --json
orfdata catalogs info ncorf-mammals-mouse --json
orfdata query --gene TPST1 --json
orfdata query --species mouse --format fasta > mouse-ncorfs.fasta
orfdata datasets query pfizer-atlas-dog-protein --unit iBAQ --format jsonl
orfdata datasets pairs canine-mammary-stroma --basis same_specimen --json
```

`--json` emits one complete JSON document. Dataset `query` and `pairs` return
JSON arrays; `--format jsonl` streams one original record per line. ORF queries
also support `--format tsv`, `fasta`, and `gtf`. Viral coverage supports TSV.
The preview limit never truncates a machine export.

Measurement previews keep native units, scales and statuses, display values to
six significant digits, and abbreviate long identifiers. Exports retain full
identifiers and numeric precision. Reference
labels distinguish scoped identifiers that happen to have the same spelling.
Mapped measurements show the original observation and its mapping attribution;
pair previews show their pairing basis, temporal relationship and target
relationship. Full reference identities and provenance remain in JSON/JSONL.

Results and exports go to stdout. Progress and errors go to stderr. `data path`
and `viral inventory download` print a plain path by default. `proteome` writes
FASTA to stdout unless `--output` is supplied; its summary goes to stderr.
`lift` writes its catalogue to `--output` and also summarizes on stderr. For
these two commands, `--json` selects the summary format on stderr.

## Migrating scripts from 0.42 or 0.43

Bare `orfdata` now shows help. Use `orfdata query --json` for the previous
unfiltered JSON catalogue export. Existing root query flags (for example,
`orfdata --gene TPST1`) still select the same records, with a human preview by
default.

Add `--json` to management commands or ORF queries whose JSON output your script
parses. Add `--format jsonl` to dataset queries or paired-measurement exports to
retain their previous streaming representation. Explicit FASTA, TSV and GTF
exports retain their formats. Python APIs and stored data formats are unchanged.

The [CLI design and audit](cli-design.md) explains the conventions and the
fixture-backed command audit that checks them in CI.
