# CLI design and audit

ORFdata's command line is for discovering data, inspecting local state, acquiring
explicitly selected resources, and exporting records for other tools. Running
`orfdata` should explain these tasks without loading a catalogue or downloading
anything.

The design follows the [Command Line Interface Guidelines](https://clig.dev/#help):
show help for bare command groups, make subcommands discoverable, keep successful
output brief, and offer explicit machine formats. The output conventions also
draw on OpenVax's OncoRef inventory summaries and MHCtools' shared table formatter.

## Conventions

- `orfdata`, `orfdata help`, and `orfdata --help` show the command index.
- `orfdata help COMMAND` and `orfdata COMMAND --help` show command help.
- Human summaries and tables are the default. `--json` selects a complete JSON
  document. Export commands also offer explicit formats such as FASTA or JSONL.
- Human tables preview 25 rows; `--limit 0` shows all rows. Machine exports are
  complete regardless of the preview limit. Wrapping is only for human output.
- Results go to stdout; progress and errors go to stderr. Proteome FASTA stays on
  stdout and its build summary stays on stderr, as does the lift summary.
- Discovery and help never acquire resources. Downloads and builds require an
  explicit acquisition command and continue to enforce source terms and pins.
- Errors identify the problem and return a status a script can branch on: 1 when
  the command fails, 2 for a mistake in the command line itself, and 3 when a
  resource it needs isn't installed. Only command-line mistakes point to command
  help; a missing installation names the command that installs it. Unexpected
  internal errors keep their traceback rather than being reported as user errors.

## Audit

The command audit tests enumerate the parser tree, exercise help and bare groups,
and run each leaf with offline fixtures. Acquisition tests use tiny files or
replace the external acquisition boundary; they do not fetch production datasets
or rebuild genome indexes. Export tests check structured payloads separately from
human previews. This keeps the audit reproducible and small enough for CI.

Every Tests workflow run attaches `cli-audit-python-VERSION` artifacts containing
the actual help and command transcripts. Run the same audit locally with:

```sh
ORFDATA_CLI_AUDIT_DIR=/tmp/orfdata-cli-audit python -m pytest tests/test_cli_audit.py
```

The audit currently exercises 43 leaf commands in both human and JSON modes, plus
the top-level guide, version, command aliases, nested help, errors, bounded
previews, and complete machine exports. Its parser enumeration fails if a new
command has no execution case. A tiny FASTA inventory file replaces the inventory
download; other acquisitions use checksum-pinned fixtures and fixture viral
records through the real handlers. These checks are not production dataset
download or biological coverage validation.
