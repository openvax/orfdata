"""The full builder joins pinned local components and refuses incomplete inputs."""

import csv
import gzip
import json
import sqlite3
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_pfizer_evidence import build
from test_pfizer_evidence import source as source
from test_pfizer_lineage import quantification

from orfdata import QuantitativeDataset, SourceAsset
from orfdata.datasets import build_dataset
from orfdata.sources import importer, search_references
from orfdata.sources.pfizer import TABLE_COLUMNS, _audit, protein_manifest
from orfdata.sources.pfizer_evidence import EvidenceStore
from orfdata.sources.pfizer_native import IMPORTER, file_sha256, import_native, selection_sha256
from orfdata.sources.pfizer_vendor import NativeQuantification

DIRECTORY = Path(__file__).parent / "data/pfizer/vendor"


def pinned(name, path, **metadata):
    return SourceAsset(name, path.resolve().as_uri(), file_sha256(path), path.stat().st_size,
                       "Explicit regression slice; see tests/data/pfizer/vendor/NOTICE.md", metadata=metadata)


def repin(manifest, name, path, **metadata):
    return replace(manifest, sources=tuple(pinned(name, path, **metadata) if a.asset_id == name else a for a in manifest.sources))


@pytest.fixture
def inputs(source, tmp_path):
    native, result = build(source[0], tmp_path)
    evidence = EvidenceStore(native)
    lineage, tables, original = quantification()
    source_rows = {(r["WorkflowID"], r["FileID"]): r["_source_rowid"]
                   for r in evidence.records("WorkflowInputFiles")}
    for run in lineage["runs"]:
        run["vendor"]["_source_rowid"] = source_rows[run["workflow_id"], run["file_id"]]
    summary = tmp_path / "protein.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, sorted(TABLE_COLUMNS), delimiter="\t")
        writer.writeheader()
        for i, channel in enumerate(c for c in lineage["channels"] if c["summary_sample"]):
            writer.writerow(dict(channel["summary_sample"], ProteinGroup="A0A5F4CET6", Protein="A0A5F4CET6",
                                 IBAQ=("NA", "0", "23.5")[i % 3], rIBAQ=("NA", "0", "0.5")[i % 3],
                                 Abundance="-Inf" if i % 3 == 1 else "2", sumIBAQ="47"))
    lineage["native_export_configuration"] = evidence.metadata["source"]["selection"]
    lineage["source_sha256"]["protein_summary"] = file_sha256(summary)
    lineage_path = tmp_path / "lineage.json"
    lineage_path.write_text(json.dumps(lineage))
    paths = {"evidence.sqlite": native, "lineage.json": lineage_path,
             "search.fasta": DIRECTORY / "membership-search-slice.fasta", "protein.tsv": summary}
    assets = tuple(pinned(name, path) for name, path in paths.items())
    reference = replace(original.reference, sequence_sha256=file_sha256(paths["search.fasta"]))
    manifest = replace(protein_manifest("dog"), dataset_id="canine-evidence-regression", version="test",
                       title="Explicit fixture selection; not a complete canine source extraction",
                       sources=assets, references=(reference,),
                       coverage=dict(search_sequences=10, native_counts=result["counts"], scope="regression fixture"),
                       processing=dict(importer=IMPORTER,
                                       selection_sha256=selection_sha256(lineage["native_export_configuration"]["plan"]),
                                       limits=dict(output_bytes=32 * 1024 * 1024, scratch_bytes=128 * 1024 * 1024),
                                       protein_summary=dict(dataset_id="explicit-summary-regression", version="test",
                                                            source_sha256=file_sha256(summary), source_counts=_audit(summary)[2])))
    converter = NativeQuantification(original.decoder, lineage, reference, manifest.dataset_id, manifest.version,
                                     source_asset="evidence.sqlite")
    values = [m for table, columns in converter.COLUMNS.items() for row in tables[table]
              for column in columns for m in converter.peptide_abundances(table, row, column)]
    return paths, manifest, sorted(values, key=lambda m: m.measurement_id)


def test_complete_build_preserves_measurements_memberships_and_summary_context(inputs, tmp_path):
    paths, manifest, expected = inputs
    dataset = importer(IMPORTER)(paths, tmp_path / "dataset.sqlite", manifest)
    assert dataset.manifest == manifest
    assert list(dataset.measurements()) == expected
    assert dataset.counts["measurement"] == 429
    evidence = dataset.source_evidence()
    assert evidence.search_reference == search_references(manifest)[0]
    assert list(evidence.records("TargetPsms")) == list(EvidenceStore(paths["evidence.sqlite"]).records("TargetPsms"))
    group = evidence.summary_group("A0A5F4CET6")
    assert group["baseline"] == dict(manifest.processing["protein_summary"])
    assert group["status"] == "matched_native_records" and len(group["native_proteins"]) == 1
    protein = group["native_proteins"][0]
    assert protein["protein"]["record"]["Accession"] == "A0A5F4CET6"
    assert protein["group_links"][0]["left"]["record"]["ProteinGroupID"] == 1
    assert evidence.summary_group("not-in-summary") is None
    with sqlite3.connect(dataset.path) as db:
        receipt = json.loads(db.execute("SELECT value FROM metadata WHERE name='canine_build'").fetchone()[0])
        assert receipt["state"] == "complete" and receipt["summary"]["rows"] == 38
        assert receipt["summary"]["IBAQ:missing"] == 13 and receipt["summary"]["IBAQ:zero"] == 13
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
        assert db.execute("SELECT count(*) FROM sqlite_master WHERE name='summary_seen'").fetchone()[0] == 0
    assert not list(tmp_path.glob("canine-evidence-*"))
    with sqlite3.connect(dataset.path) as db:
        db.execute("DELETE FROM metadata WHERE name='canine_build'")
    with pytest.raises(ValueError, match="build is incomplete"):
        QuantitativeDataset(dataset.path)


def test_gzip_component_and_managed_cache_build_are_offline(inputs, tmp_path):
    paths, manifest, expected = inputs
    paths = dict(paths)
    original = paths["evidence.sqlite"]
    compressed = tmp_path / "component.sqlite.gz"
    compressed.write_bytes(gzip.compress(original.read_bytes(), mtime=0))
    manifest = repin(manifest, "evidence.sqlite", compressed, compression="gzip",
                     uncompressed_size=original.stat().st_size, uncompressed_sha256=file_sha256(original))
    paths["evidence.sqlite"] = compressed
    dataset = build_dataset(manifest, lambda p, o: import_native(p, o, manifest),
                            cache_root=tmp_path / "cache", source_paths=paths)
    assert list(dataset.measurements()) == expected
    assert dataset.source_evidence().summary_group("A0A5F4CET6")["status"] == "matched_native_records"
    assert not list(dataset.path.parent.glob("canine-evidence-*"))


@pytest.mark.parametrize("failure", ["expanded_limit", "expanded_hash", "truncated"])
def test_compressed_component_checks_expanded_bytes_and_cleans_private_spool(inputs, tmp_path, failure):
    paths, manifest, _ = inputs
    paths = dict(paths)
    original = paths["evidence.sqlite"]
    data = gzip.compress(original.read_bytes(), mtime=0)
    compressed = tmp_path / "component.sqlite.gz"
    compressed.write_bytes(data[:-8] if failure == "truncated" else data)
    manifest = repin(manifest, "evidence.sqlite", compressed, compression="gzip",
                     uncompressed_size=original.stat().st_size - (1 if failure == "expanded_limit" else 0),
                     uncompressed_sha256="0" * 64 if failure == "expanded_hash" else file_sha256(original))
    paths["evidence.sqlite"] = compressed
    output = tmp_path / "failed.sqlite"
    with pytest.raises((ValueError, EOFError)):
        import_native(paths, output, manifest)
    assert not output.exists()
    assert not list(tmp_path.glob("canine-evidence-*"))


def test_release_manifest_factory_refuses_a_partial_fixture_selection(inputs, tmp_path, monkeypatch):
    from tools import build_dog_dataset

    paths, manifest, _ = inputs
    options = dict(version="test", evidence_url="https://example.invalid/evidence.sqlite",
                   lineage_url="https://example.invalid/lineage.json", output_bytes=32 * 1024 * 1024,
                   scratch_bytes=128 * 1024 * 1024)
    with pytest.raises(ValueError, match="entire pinned"):
        build_dog_dataset.manifest_for(paths["evidence.sqlite"], paths["lineage.json"], **options)
    # Only this explicit test substitutes its own extraction plan and baseline.
    # The production factory always reads the entire 42-table pinned selection.
    evidence = EvidenceStore(paths["evidence.sqlite"])
    selection = tmp_path / "fixture-plan.json"
    selection.write_text(json.dumps(evidence.metadata["source"]["selection"]["plan"]))
    monkeypatch.setattr(build_dog_dataset, "PLAN_PATH", selection)
    baseline = replace(manifest, dataset_id="summary-regression", processing=dict(importer="pfizer-protein-v1"),
                       coverage=dict(search_sequences=10, source_counts=dict(manifest.processing["protein_summary"]["source_counts"])))
    monkeypatch.setattr(build_dog_dataset, "protein_manifest", lambda species: baseline)
    generated = build_dog_dataset.manifest_for(paths["evidence.sqlite"], paths["lineage.json"], **options)
    assert generated.coverage["accession"] == "MSV000099744"
    assert generated.coverage["proteomexchange_accession"] == "PXD070290"
    dataset = import_native(paths, tmp_path / "generated.sqlite", generated)
    assert dataset.counts["measurement"] == 429
    assert dataset.source_evidence().summary_group("A0A5F4CET6")["baseline"]["dataset_id"] == baseline.dataset_id

    from orfdata import DatasetRegistry
    from tools.package_sqlite_archive import archive_verified

    archived = archive_verified(paths["evidence.sqlite"], tmp_path / "parts", member="evidence.sqlite",
                                url_prefix="https://example.invalid/release/", part_bytes=65536, max_bytes=16 * 1024**2)
    distributed = build_dog_dataset.manifest_for(paths["evidence.sqlite"], paths["lineage.json"],
                                                evidence_archive=archived["archive"], **options)
    asset = next(a for a in distributed.sources if a.asset_id == "evidence.sqlite")
    assert asset.sha256 == file_sha256(paths["evidence.sqlite"])
    assert asset.to_dict()["metadata"]["archive"] == archived["archive"]
    definitions = {distributed.dataset_id: dict(default_version=distributed.version, versions={distributed.version: dict(
        manifest=distributed.to_dict(), assets={}, distribution=dict(access_terms="Fixture", redistributable=False))})}
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "distributed-cache")
    supplied = dict(paths)
    supplied["evidence.sqlite"] = [tmp_path / "parts" / p["url"].rsplit("/", 1)[-1] for p in archived["archive"]["parts"]]
    rebuilt = registry.build(distributed.dataset_id, source_paths=supplied, download_options=dict(show_progress=False))
    assert rebuilt.counts == dataset.counts
    assert list(rebuilt.measurements()) == list(dataset.measurements())


@pytest.mark.parametrize("failure", ["hash", "selection", "incomplete", "count", "schema", "seed_limit", "lineage", "lineage_rowid", "summary", "scratch", "output"])
def test_invalid_inputs_never_leave_a_consumer_dataset(inputs, tmp_path, failure):
    paths, manifest, _ = inputs
    processing = manifest.to_dict()["processing"]
    if failure == "hash":
        paths["protein.tsv"].write_text("changed")
    elif failure == "selection":
        processing["selection_sha256"] = "0" * 64
    elif failure in ("incomplete", "count", "schema", "seed_limit"):
        with sqlite3.connect(paths["evidence.sqlite"]) as db:
            if failure == "incomplete":
                db.execute("UPDATE evidence_metadata SET value='\"building\"' WHERE name='state'")
            elif failure == "count":
                db.execute("DELETE FROM evidence_record WHERE table_i=(SELECT i FROM evidence_table WHERE name='WorkflowInputFiles') "
                           "AND source_rowid=1")
            elif failure == "schema":
                db.execute("CREATE VIEW unexpected_source_sql AS SELECT * FROM evidence_record")
            else:
                db.execute("INSERT INTO evidence_metadata VALUES ('large_fixture',?)", (json.dumps('x' * 8 * 1024**2),))
                processing["limits"]["output_bytes"] = 4 * 1024**2
        manifest = repin(manifest, "evidence.sqlite", paths["evidence.sqlite"])
    elif failure in ("lineage", "lineage_rowid"):
        lineage = json.loads(paths["lineage.json"].read_text())
        lineage["runs"][0]["vendor"]["FileSize" if failure == "lineage" else "_source_rowid"] += 1
        paths["lineage.json"].write_text(json.dumps(lineage))
        manifest = repin(manifest, "lineage.json", paths["lineage.json"])
    elif failure == "summary":
        processing["protein_summary"]["source_counts"]["rows"] += 1
    elif failure == "scratch":
        processing["limits"]["scratch_bytes"] = 4096
    else:
        processing["limits"]["output_bytes"] = 4096
    manifest = replace(manifest, processing=processing)
    output = tmp_path / "failed.sqlite"
    with pytest.raises((ValueError, sqlite3.DatabaseError)):
        import_native(paths, output, manifest)
    assert not output.exists()
    assert not list(tmp_path.glob("canine-evidence-*"))


def test_interruption_and_existing_output_preservation(inputs, tmp_path):
    paths, manifest, _ = inputs
    output = tmp_path / "dataset.sqlite"

    def interrupt(state):
        assert not output.exists()
        raise KeyboardInterrupt("simulated interrupted copy")

    with pytest.raises(KeyboardInterrupt):
        import_native(paths, output, manifest, checkpoint=interrupt)
    assert not output.exists()
    output.write_bytes(b"preserve existing local work")
    with pytest.raises(FileExistsError):
        import_native(paths, output, manifest)
    assert output.read_bytes() == b"preserve existing local work"


def test_one_output_reservation_builds_privately_without_disk_journal(inputs, tmp_path, monkeypatch):
    from orfdata.sources import pfizer_native

    paths, manifest, expected = inputs
    output = tmp_path / "dataset.sqlite"
    cap = manifest.processing["limits"]["output_bytes"]
    # Enough for one capped output and small metadata journal, but not two.
    reserve = cap + pfizer_native.METADATA_BYTES + 1024**2
    assert reserve < 2 * cap
    monkeypatch.setattr(pfizer_native.shutil, "disk_usage", lambda _: SimpleNamespace(free=reserve))
    input_hashes = {name: file_sha256(path) for name, path in paths.items()}
    observed = []

    def inspect(state):
        assert not output.exists()
        assert not list(tmp_path.glob("canine-evidence-*/*-journal"))
        observed.append(state)

    dataset = import_native(paths, output, manifest, checkpoint=inspect)
    assert observed and list(dataset.measurements()) == expected
    assert {name: file_sha256(path) for name, path in paths.items()} == input_hashes
    with sqlite3.connect(output) as db:
        receipt = json.loads(db.execute("SELECT value FROM metadata WHERE name='canine_build'").fetchone()[0])
        assert receipt["storage"]["required_free_bytes"] == reserve
    assert not list(tmp_path.glob("canine-evidence-*"))
    monkeypatch.setattr(pfizer_native.shutil, "disk_usage", lambda _: SimpleNamespace(free=cap))
    with pytest.raises(ValueError, match="Insufficient free disk"):
        import_native(paths, tmp_path / "insufficient.sqlite", manifest)
    assert not (tmp_path / "insufficient.sqlite").exists()


def test_completed_build_never_replaces_a_racing_destination(inputs, tmp_path):
    paths, manifest, _ = inputs
    output = tmp_path / "dataset.sqlite"

    def race(state):
        if not output.exists():
            output.write_bytes(b"another writer owns this destination")

    with pytest.raises(FileExistsError):
        import_native(paths, output, manifest, checkpoint=race)
    assert output.read_bytes() == b"another writer owns this destination"
    assert not list(tmp_path.glob("canine-evidence-*"))


@pytest.mark.parametrize("cap", [True, 0, 4095, 4096.0, 4096])
def test_common_dataset_caps_apply_before_schema_creation(inputs, tmp_path, cap):
    _, manifest, _ = inputs
    output = tmp_path / "tiny.sqlite"
    with pytest.raises((ValueError, sqlite3.DatabaseError)):
        QuantitativeDataset.create(output, manifest, samples=(), assays=(), measurements=(), max_bytes=cap)
    assert not output.exists()
