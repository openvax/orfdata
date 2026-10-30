"""Exercise every CLI command offline, including human output, exports and help.

The registries, genome pins and viral records are tiny test fixtures. Real
parsers, acquisition handlers, hashing, imports, queries and exports still run.
Only the inventory download and external source lookups use fixture boundaries.
"""

import argparse
import json
import os
import socket
from pathlib import Path

import pytest
from test_catalog_registry import definition as catalog_definition
from test_catalog_registry import setup_release
from test_genomes import slice_pins
from test_liftover import chain as chain
from test_liftover import pinned as pinned
from test_liftover import shifted
from test_measurements import definition as dataset_definition
from test_pairing import create
from test_pairing import paired as paired
from test_proteome import genome as genome
from test_proteome import make_catalog

from orfdata import CatalogRegistry, DatasetManifest, DatasetRegistry, ORFCatalog
from orfdata.cli import COMMANDS, run
from orfdata.cli_parser import ArgumentParser
from orfdata.sources.pfizer import import_protein

GROUPS = ("data", "catalogs", "datasets", "genomes", "chains", "viral")


def save_transcript(name, entries, tmp_path):
    directory = Path(os.environ.get("ORFDATA_CLI_AUDIT_DIR", tmp_path))
    directory.mkdir(parents=True, exist_ok=True)
    (directory / name).write_text("\n".join(entries))


def command_paths(monkeypatch):
    """Enumerate the actual parsers, so a new command must join the execution audit."""
    found = set()

    class Captured(Exception):
        pass

    def capture(parser, *args, **kwargs):
        def walk(parser, path):
            children = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
            if not children:
                found.add(tuple(path))
            for action in children:
                for name, child in action.choices.items():
                    walk(child, [*path, name])
        walk(parser, parser.prog.split()[1:])
        raise Captured

    with monkeypatch.context() as patch:
        patch.setattr(ArgumentParser, "parse_args", capture)
        for name in COMMANDS:
            if name != "help":
                with pytest.raises(Captured):
                    run([name, "--help"])
    return found


def invoke(argv):
    try:
        run(argv)
    except SystemExit as error:
        assert error.code == 0, argv


def test_every_command_is_discoverable_without_io(monkeypatch, capsys, tmp_path):
    def refuse(*args, **kwargs):
        pytest.fail("Help must not open catalogues, inspect caches, or acquire data")

    monkeypatch.setattr(CatalogRegistry, "__init__", refuse)
    monkeypatch.setattr(DatasetRegistry, "__init__", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    paths = command_paths(monkeypatch)
    transcript = []
    for argv in ([], ["--help"], ["-h"], ["help"]):
        run(argv, catalog_loader=refuse)
        output = capsys.readouterr()
        transcript.append("$ orfdata " + " ".join(argv) + "\n" + output.out)
        assert "data status" in output.out and "commands:" in output.out
        assert "--json" in output.out and len(output.out.splitlines()) <= 35
        assert not output.err
    for path in sorted(paths | {(name,) for name in GROUPS} | {("viral", "inventory")}):
        for argv in ([*path, "--help"], ["help", *path]):
            invoke(argv)
            output = capsys.readouterr()
            transcript.append("$ orfdata " + " ".join(argv) + "\n" + output.out)
            assert "usage: orfdata " + " ".join(path) in output.out
            assert "--help" in output.out and not output.err
    for path in [(name,) for name in GROUPS] + [("viral", "inventory")]:
        invoke(list(path))
        output = capsys.readouterr()
        assert "commands:" in output.out and not output.err
    for alias in ("download", "downloads"):
        invoke([alias])
        assert "fetch" in capsys.readouterr().out
    run(["--version"])
    from orfdata.version import __version__
    assert capsys.readouterr().out.strip() == f"orfdata {__version__}"
    save_transcript("help.txt", transcript, tmp_path)


@pytest.fixture
def audit_cases(tmp_path, monkeypatch, paired, chain, pinned, genome, viral_panel):
    # Refuse real network connections; local pinned files and fixture NCBI records suffice.
    def refuse(*args, **kwargs):
        pytest.fail("The CLI audit must not access the network")
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setenv("COLUMNS", "100")

    manifest, source, catalog = setup_release(tmp_path)
    artifact = tmp_path / "catalog.json.gz"
    build_registry = CatalogRegistry(catalog_definition(manifest), cache_root=tmp_path / "catalog-source")
    built = build_registry.build(manifest.catalog_id, source_paths={"source.json": source},
                                download_options={"show_progress": False})
    built.save(artifact)
    catalog_defs = catalog_definition(manifest, artifact)
    name = manifest.catalog_id
    original_catalog_init = CatalogRegistry.__init__
    monkeypatch.setattr(CatalogRegistry, "__init__", lambda self, definitions=None, **kwargs:
                        original_catalog_init(self, catalog_defs if definitions is None else definitions, **kwargs))

    paired_dataset = create(tmp_path, paired)
    root = Path(__file__).parent / "data/pfizer/dog"
    protein_manifest = DatasetManifest.from_dict(json.loads((root / "manifest.json").read_text()))
    paths = {asset.asset_id: root / asset.asset_id for asset in protein_manifest.sources}
    protein = import_protein(paths, tmp_path / "protein.sqlite", protein_manifest)
    dataset_defs = {**dataset_definition(paired_dataset), **dataset_definition(protein)}
    original_dataset_init = DatasetRegistry.__init__
    monkeypatch.setattr(DatasetRegistry, "__init__", lambda self, definitions=None, **kwargs:
                        original_dataset_init(self, dataset_defs if definitions is None else definitions, **kwargs))
    protein_name = protein.manifest.dataset_id
    source_args = [arg for key, path in paths.items() for arg in ("--source-path", f"{key}={path}")]
    reference_target = next(iter(DatasetRegistry().search_references(protein_name)))
    # Resolve a native accession directly from the pinned fixture FASTA.
    from orfdata import SearchReference
    search = SearchReference(paths["search.fasta"], reference_target)
    accession = next(search.entries()).target.identifier

    from orfdata import genomes, proteome_cli, viral_cli
    genome_paths, genome_pins = slice_pins(tmp_path)
    monkeypatch.setattr(genomes, "_pins", lambda: genome_pins)
    monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda *args: genome)
    inventory = tmp_path / "inventory.txt"
    inventory.write_text("fixture inventory\n")
    monkeypatch.setattr(viral_cli, "download_viral_inventory", lambda **kwargs: inventory)
    genome_args = [arg for key, path in genome_paths.items() for arg in ("--source-path", f"{key}={path}")]
    chain_args = ["GRCm38", "GRCm39", "--accept-terms", "--source-path", str(chain.path)]
    native = make_catalog()
    old = ORFCatalog([shifted(record) for record in native.records], dict(native.metadata, genome_build="GRCm38"))
    old_path = tmp_path / "old.json"
    old.save(old_path)

    # Order acquisitions before installed-state queries. Every command runs twice,
    # once as a human would use it and once with an explicit complete JSON export.
    return [
        (("query",), ["--catalog", str(source)]),
        (("catalogs", "list"), []),
        (("catalogs", "info"), [name]),
        (("catalogs", "download"), [name, "--source-path", f"catalog.json.gz={artifact}"]),
        (("catalogs", "build"), [name, "--source-path", f"source.json={source}"]),
        (("catalogs", "status"), []),
        (("catalogs", "inspect"), [name]),
        (("catalogs", "query"), [name]),
        (("datasets", "list"), []),
        (("datasets", "info"), [protein_name]),
        (("datasets", "download"), ["fixture", "--source-path", str(paired_dataset.path), "--quiet"]),
        (("datasets", "build"), [protein_name, *source_args, "--quiet"]),
        (("datasets", "inspect"), ["fixture"]),
        (("datasets", "query"), ["fixture"]),
        (("datasets", "pairs"), ["fixture", "--basis", "same_material"]),
        (("datasets", "references"), [protein_name]),
        (("datasets", "download-reference"), [protein_name, "--source-path", str(paths["search.fasta"]), "--quiet"]),
        (("datasets", "inspect-reference"), [protein_name]),
        (("datasets", "query-reference"), [protein_name, "--target", accession]),
        (("genomes", "list"), []),
        (("genomes", "install"), ["mouse-slice", *genome_args]),
        (("chains", "list"), []),
        (("chains", "download"), chain_args),
        (("lift",), ["--catalog", str(old_path), "--to", "GRCm39", "--chain", str(chain.path),
                      "--output", str(tmp_path / "lifted.json")]),
        (("proteome",), ["--species", "mouse", "--release", "106"]),
        (("setup",), ["--species", "mouse", "--catalogs", name, "--collections", "none", "--accept-terms"]),
        (("data", "list"), []),
        (("data", "status"), []),
        (("data", "info"), [name]),
        (("data", "path"), [name]),
        (("data", "fetch"), [name, "--accept-terms", "--source-path", f"catalog.json.gz={artifact}"]),
        (("data", "catalog"), [name, "--source-path", f"catalog.json.gz={artifact}"]),
        (("data", "dataset"), [protein_name, "--source-path", f"dataset.sqlite={protein.path}", "--quiet"]),
        (("data", "reference"), [protein_name, "--source-path", str(paths["search.fasta"]), "--quiet"]),
        (("data", "genome"), ["mouse-slice", *genome_args]),
        (("data", "chains"), chain_args),
        (("data", "viral"), []),
        (("viral", "coverage"), []),
        (("viral", "contexts"), []),
        (("viral", "subsets"), []),
        (("viral", "summary"), []),
        (("viral", "completeness"), []),
        (("viral", "inventory", "download"), []),
    ]


def test_every_command_executes_with_human_and_json_output(audit_cases, monkeypatch, capsys, tmp_path):
    assert {path for path, _ in audit_cases} == command_paths(monkeypatch)
    capsys.readouterr()  # Discard fixture construction progress.
    transcript = []
    for path, options in audit_cases:
        for machine in (False, True):
            argv = [*path, *options, *(["--json"] if machine else [])]
            invoke(argv)
            captured = capsys.readouterr()
            output = captured.err if path in (("lift",), ("proteome",)) else captured.out
            transcript.append("$ orfdata " + " ".join(argv) + "\n" + captured.out + captured.err)
            save_transcript("commands.txt", transcript, tmp_path)
            assert output.strip(), argv
            assert "Traceback" not in output, argv
            if machine:
                assert isinstance(json.loads(output), (list, dict)), argv
            elif path not in (("data", "path"), ("viral", "inventory", "download")):
                assert not output.lstrip().startswith(("{", "[")), argv
                assert len(output.splitlines()) <= 100, argv
            if path == ("proteome",):
                assert captured.out.startswith(">") and "{" not in captured.out
            elif path == ("lift",):
                assert not captured.out


@pytest.mark.parametrize("argv", [
    ["not-a-command"], ["data", "not-a-command"], ["data", "fetch"],
    ["catalogs", "info"], ["datasets", "query"], ["genomes", "install"],
    ["lift"], ["setup"], ["proteome"], ["query", "--limit", "-1"],
])
def test_errors_are_concise_and_point_to_help(argv, capsys):
    with pytest.raises(SystemExit) as error:
        run(argv)
    assert error.value.code == 2
    output = capsys.readouterr()
    assert not output.out and "error:" in output.err and "--help' for usage" in output.err
    assert "Traceback" not in output.err and len(output.err.splitlines()) <= 5


def test_machine_exports_are_complete_despite_human_preview_limit(tmp_path, capsys):
    catalog = make_catalog()
    path = tmp_path / "catalog.json"
    catalog.save(path)
    argv = ["query", "--catalog", str(path), "--limit", "1"]
    run(argv)
    human = capsys.readouterr().out
    assert f"{len(catalog.records)} matching ORFs" in human
    assert "Showing 1 of" in human and "--limit 0" in human
    run(argv + ["--json"])
    assert len(json.loads(capsys.readouterr().out)["records"]) == len(catalog.records)
    run(argv + ["--format", "fasta"])
    fasta = capsys.readouterr().out
    assert fasta.startswith(">") and fasta.count("\n>") + 1 == len(catalog.records)
    assert "Showing" not in fasta
    run(argv + ["--format", "tsv"])
    assert len(capsys.readouterr().out.splitlines()) == len(catalog.records) + 1
    run(argv + ["--format", "gtf"])
    gtf = capsys.readouterr().out
    assert "\tCDS\t" in gtf and "Showing" not in gtf


def test_streaming_previews_do_not_exhaust_large_query_results(capsys):
    from orfdata.cli_output import print_record_stream
    consumed = []

    def rows():
        for number in range(100):
            consumed.append(number)
            yield dict(value=number, unit="TPM", status="measured")

    print_record_stream(rows(), format="table", columns=(("VALUE", "value"),), limit=3)
    assert consumed == [0, 1, 2, 3]
    assert "Showing the first 3" in capsys.readouterr().out
    print_record_stream(rows(), format="json", columns=(), limit=3)
    exported = json.loads(capsys.readouterr().out)
    assert len(exported) == 100 and exported[0]["value"] == 0
    print_record_stream(rows(), format="jsonl", columns=(), limit=3)
    assert [json.loads(line) for line in capsys.readouterr().out.splitlines()] == exported


def test_position_preview_shows_annotations_and_no_translation_claim(capsys):
    run(["query", "--position", "7:66205480", "--genome-build", "GRCh38", "--ref", "A", "--alt", "C"])
    output = capsys.readouterr().out
    assert "CONTEXT POSITION" in output and "REFERENCE" in output and "ALTERNATE" in output
    assert "translation change is not measured" in output
    assert "GAGGCCAGG[ATG]CCGTCC" in output and "GAGGCCCGG[ATG]CCGTCC" in output


def test_tables_wrap_without_hiding_identifiers_or_values(monkeypatch, capsys):
    from orfdata.cli_output import print_table
    monkeypatch.setenv("COLUMNS", "80")
    print_table(("CONTEXT", "DESCRIPTION", "VALUE"), [
        ("cancer-transcriptome-screening", "Long explanatory text " * 8, 0),
        ("missing-value", "Unknown remains different from measured zero", None),
    ])
    output = capsys.readouterr().out
    assert max(map(len, output.splitlines())) <= 80
    assert "cancer-transcriptome-screening" in output and "  0" in output and "  -" in output


def test_aliases_produce_identical_results(capsys):
    results = []
    for alias in ("data", "download", "downloads"):
        run([alias, "list", "--kind", "viral", "--json"])
        results.append(json.loads(capsys.readouterr().out))
    assert results[0] == results[1] == results[2]


@pytest.mark.parametrize("argv,target", [
    (["catalogs", "download", "fixture"], "orfdata.catalogs.CatalogRegistry.download"),
    (["catalogs", "build", "fixture"], "orfdata.catalogs.CatalogRegistry.build"),
    (["datasets", "download", "fixture"], "orfdata.datasets.DatasetRegistry.download"),
    (["datasets", "build", "fixture"], "orfdata.datasets.DatasetRegistry.build"),
    (["genomes", "install", "fixture"], "orfdata.genomes_cli.install_genome"),
    (["chains", "download", "GRCm38", "GRCm39"], "orfdata.liftover_cli.download_chain"),
    (["setup", "--species", "mouse"], "orfdata.install.install_species"),
    (["viral", "inventory", "download"], "orfdata.viral_cli.download_viral_inventory"),
])
def test_acquisition_failures_show_actionable_errors(argv, target, monkeypatch, capsys):
    def fail(*args, **kwargs):
        raise RuntimeError("Source unavailable; retry after checking your connection")
    monkeypatch.setattr(target, fail)
    with pytest.raises(SystemExit) as error:
        run(argv)
    assert error.value.code == 2
    output = capsys.readouterr()
    assert not output.out and "Source unavailable" in output.err
    assert "--help' for usage" in output.err and "Traceback" not in output.err
