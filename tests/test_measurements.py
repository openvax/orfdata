"""Units, censoring, provenance and reference boundaries must survive round trips."""

import hashlib
import json
import sqlite3
from dataclasses import replace
from types import SimpleNamespace

import pytest

from orfdata import (
    DatasetManifest,
    DatasetRegistry,
    Measurement,
    QuantitativeDataset,
    ScopedTarget,
    TargetMapping,
    build_dataset,
)
from orfdata.cli import run
from orfdata.measurements import MEASUREMENT_STATUSES


def write_dataset(path, experiment, **kwargs):
    manifest, samples, assays, observations = experiment
    return QuantitativeDataset.create(path, manifest, samples=samples, assays=assays,
                                      measurements=kwargs.pop("measurements", observations), **kwargs)


def definition(dataset):
    path = dataset.path
    return {dataset.manifest.dataset_id: dict(default_version=dataset.manifest.version, versions={
        dataset.manifest.version: dict(
            manifest=dataset.manifest.to_dict(),
            distribution=dict(access_terms="Synthetic Apache-2.0 test data", redistributable=True),
            assets={"dataset.sqlite": dict(url=path.as_uri(),
                                           sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                           size=path.stat().st_size, decompress=False)},
        )})}


def test_four_species_id_collisions_and_scoped_queries(tmp_path, experiment):
    dataset = write_dataset(tmp_path / "data.sqlite", experiment)
    observations = list(dataset.measurements())
    assert len(observations) == 4
    assert {o.target.reference.taxon_id for o in observations} == {9606, 10090, 9615, 4932}
    assert {o.dataset_id for o in observations} == {"fixture"}
    assert dataset.counts == dict(sample=4, assay=4, mapping=0, target=4, measurement=4,
                                  assay_pair=0, target_relation=0, reference_mapping=0)
    for expected in experiment[3]:
        assert list(dataset.measurements(target=expected.target)) == [expected]
        assert list(dataset.measurements(reference=expected.target.reference, sample_id="sample")) == [expected]
    with pytest.raises(ValueError, match="explicit reference"):
        list(dataset.measurements(assay_id="rna"))
    with pytest.raises(ValueError, match="not present"):
        list(dataset.measurements(target=replace(observations[0].target,
                                                reference=replace(observations[0].target.reference,
                                                                  annotation_version="111"))))
    with pytest.raises(ValueError, match="disagree"):
        list(dataset.measurements(target=observations[0].target, reference=observations[1].target.reference))


@pytest.mark.parametrize("status", MEASUREMENT_STATUSES)
def test_zero_and_missingness_roundtrip(tmp_path, experiment, status):
    original = experiment[3][0]
    observation = replace(original, value=0 if status == "measured" else None, status=status,
                          detection_limit=0.3, provenance=dict(source_asset="input.tsv", native="NA"))
    dataset = write_dataset(tmp_path / "missing.sqlite", experiment, measurements=[observation])
    assert list(dataset.measurements())[0] == observation
    assert list(dataset.measurements(status=status)) == [observation]
    if status != "measured":
        with pytest.raises(ValueError, match="must be None"):
            replace(observation, value=0)


@pytest.mark.parametrize("quantity,unit,kind", (
    ("rna_abundance", "counts", "gene"), ("rna_abundance", "FPKM", "transcript"),
    ("ribosome_occupancy", "RPF_reads_per_nt", "orf"),
    ("translation_efficiency", "source_defined_RPF_RNA_ratio", "orf"),
    ("protein_abundance", "iBAQ", "protein_group"), ("protein_abundance", "copies/cell", "protein"),
    ("protein_abundance", "TMT_ratio", "protein"), ("protein_abundance", "LFQ", "protein"),
    ("peptide_abundance", "intensity", "peptide"), ("spectrum_count", "PSM", "orf"),
    ("turnover", "1/hour", "protein"),
))
def test_native_quantities_and_targets(experiment, quantity, unit, kind):
    observation = replace(experiment[3][0], quantity=quantity, unit=unit,
                          target=replace(experiment[3][0].target, kind=kind),
                          value=1.25, uncertainty=0.1, uncertainty_type="standard_error")
    assert Measurement.from_dict(json.loads(json.dumps(observation.to_dict()))) == observation
    log = replace(observation, scale="log", log_base=2, value=-2)
    assert Measurement.from_dict(log.to_dict()) == log
    assert log.value == -2  # No implicit exponentiation or unit conversion.


@pytest.mark.parametrize("changes", (
    {"value": float("nan")}, {"value": float("inf")}, {"value": True},
    {"value": None}, {"scale": "log", "log_base": 1}, {"log_base": 2},
    {"uncertainty": 0.1}, {"uncertainty": -1, "uncertainty_type": "sd"},
    {"origin": "derived"}, {"provenance": {"bad": float("nan")}},
    {"input_measurements": ("value",)}, {"derivation": {"method": "wrong"}},
))
def test_invalid_values_are_rejected(experiment, changes):
    with pytest.raises((ValueError, TypeError)):
        replace(experiment[3][0], **changes)


def test_metadata_is_deeply_immutable_and_json_independent(experiment):
    data = {"source_asset": "input.tsv", "rows": [1, 2], "nested": {"a": [3]}}
    measurement = replace(experiment[3][0], provenance=data)
    data["rows"].append(99)
    data["nested"]["a"].append(4)
    assert measurement.provenance["rows"] == (1, 2)
    with pytest.raises(TypeError):
        measurement.provenance["nested"]["a"] = (9,)
    exported = measurement.to_dict()
    exported["provenance"]["nested"]["a"].append(9)
    assert measurement.provenance["nested"]["a"] == (3,)
    assert DatasetManifest.from_dict(experiment[0].to_dict()) == experiment[0]


def test_protein_group_mapping_remains_ambiguous(tmp_path, experiment):
    context = experiment[0].references[0]
    candidates = tuple(ScopedTarget(context, "protein", identifier, "source-fasta:1") for identifier in ("p1", "p2"))
    mapping = TargetMapping("m", "vendor", "p1;p2", context, candidates,
                            "ambiguous", "shared_peptide", evidence={"peptide": "MK"})
    assert TargetMapping.from_dict(mapping.to_dict()) == mapping
    group = ScopedTarget(context, "protein_group", "p1;p2", "fixture:1")
    observation = replace(experiment[3][0], target=group, quantity="protein_abundance", unit="iBAQ")
    dataset = write_dataset(tmp_path / "groups.sqlite", experiment,
                            measurements=[observation], mappings=[mapping])
    assert list(dataset.mappings()) == [mapping]
    assert list(dataset.measurements(target=candidates[0])) == []
    assert list(dataset.measurements(target=group)) == [observation]
    with pytest.raises(ValueError, match="disagrees"):
        replace(mapping, status="unique")
    with pytest.raises(ValueError, match="cross reference"):
        replace(mapping, targets=(candidates[0], replace(candidates[1], reference=experiment[0].references[1])))


def test_mapping_selectors_use_the_same_models_across_species(tmp_path, experiment):
    mappings = []
    for context in experiment[0].references:
        targets = tuple(ScopedTarget(context, "protein", name, "source-fasta") for name in ("p1", "p2"))
        mappings.append(TargetMapping("same", "amino-acid-sequence", "PEPTIDE", context, targets,
                                      "ambiguous", "source-reported assignment", evidence={"basis": "reported"}))
    dataset = write_dataset(tmp_path / "mappings.sqlite", experiment, mappings=mappings)
    assert set(m.reference.taxon_id for m in dataset.mappings(source_identifier="PEPTIDE")) == {
        context.taxon_id for context in experiment[0].references}
    for mapping in mappings:
        selected = list(dataset.mappings(reference=mapping.reference, source_identifier="PEPTIDE",
                                         source_namespace="amino-acid-sequence", target=mapping.targets[0]))
        assert selected == [mapping]
        assert len(selected[0].targets) == 2
        assert list(dataset.mappings(target=replace(mapping.targets[0], identifier="missing"))) == []
    with pytest.raises(ValueError, match="disagree"):
        list(dataset.mappings(reference=mappings[0].reference, target=mappings[1].targets[0]))
    assert list(dataset.mappings(source_namespace="another-source")) == []


def test_incomplete_source_assignment_retains_known_candidates_and_reason(experiment):
    context = experiment[0].references[0]
    target = ScopedTarget(context, "protein", "p1", "source-fasta")
    mapping = TargetMapping("partial", "vendor", "record-1", context, (target,), "unresolved",
                            "source-reported assignment", "Another reported endpoint is missing")
    assert TargetMapping.from_dict(mapping.to_dict()) == mapping
    with pytest.raises(ValueError, match="requires a reason"):
        replace(mapping, reason=None)


@pytest.mark.parametrize("changes", (
    {"assay_id": "unknown"}, {"sample_id": "unknown"},
    {"dataset_version": "2"}, {"provenance": {"source_asset": "unknown"}}, {"provenance": {}},
))
def test_invalid_import_is_not_published(tmp_path, experiment, changes):
    path = tmp_path / "bad.sqlite"
    with pytest.raises((ValueError, sqlite3.IntegrityError)):
        write_dataset(path, experiment, measurements=[replace(experiment[3][0], **changes)])
    assert not path.exists()


def test_derived_lineage_cannot_dangle_or_cross_species(tmp_path, experiment):
    original = experiment[3][0]
    derived = replace(original, measurement_id="derived", origin="derived", value=original.value,
                      input_measurements=("value",), derivation={"method": "identity for test"}, provenance={})
    dataset = write_dataset(tmp_path / "derived.sqlite", experiment, measurements=[original, derived])
    assert list(dataset.measurements(target=original.target)) == [derived, original]
    for observations in ([derived, original], [experiment[3][1], derived]):
        with pytest.raises(ValueError, match="precede outputs"):
            write_dataset(tmp_path / "invalid.sqlite", experiment, measurements=observations)


def test_existing_and_missing_paths_are_never_overwritten(tmp_path, experiment):
    path = tmp_path / "existing"
    path.write_text("owned by user")
    with pytest.raises(FileExistsError):
        write_dataset(path, experiment)
    assert path.read_text() == "owned by user"
    missing = tmp_path / "missing"
    with pytest.raises(sqlite3.OperationalError):
        QuantitativeDataset(missing)
    assert not missing.exists()


def test_registry_download_corruption_repair_and_offline_open(tmp_path, experiment, monkeypatch):
    from datacache import FileValidationError

    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    registry = DatasetRegistry(definition(dataset), cache_root=tmp_path / "cache")
    assert registry.list(species="dog") == (experiment[0],)
    assert registry.list(taxon_id=9685) == ()
    assert registry.list(assay="mhc_ms") == ()
    assert not (tmp_path / "cache").exists()
    assert registry.inspect("fixture").status == "missing"
    with pytest.raises(FileNotFoundError):
        registry.open("fixture")
    paths = registry.download("fixture", show_progress=False)
    installed = registry.open("fixture")
    assert list(installed.measurements()) == list(dataset.measurements())
    old_path = installed.path
    with old_path.open("r+b") as handle:
        handle.seek(-1, 2)
        last = handle.read(1)
        handle.seek(-1, 2)
        handle.write(bytes([last[0] ^ 1]))  # Same size, different content.
    assert registry.inspect("fixture").status == "invalid"
    with pytest.raises(FileValidationError):
        registry.open("fixture")
    with pytest.raises(FileValidationError):
        registry.download("fixture", show_progress=False)
    replacement = registry.download("fixture", force=True, show_progress=False)
    assert replacement != paths
    assert old_path.exists()
    monkeypatch.setattr("datacache.download.fetch_file", lambda *a, **k: pytest.fail("network/cache write"))
    assert registry.inspect("fixture").status == "available"
    assert len(list(registry.open("fixture").measurements())) == 4
    assert registry.download("fixture", show_progress=False) == replacement


def test_registry_rejects_wrong_manifest_and_changed_release(tmp_path, experiment):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    definitions = definition(dataset)
    definitions["fixture"]["versions"]["1"]["manifest"]["title"] = "Changed source scope"
    registry = DatasetRegistry(definitions, cache_root=tmp_path / "cache")
    registry.download("fixture", show_progress=False)
    with pytest.raises(ValueError, match="manifest disagrees"):
        registry.open("fixture")
    correct = DatasetRegistry(definition(dataset), cache_root=tmp_path / "cache")
    assert correct.inspect("fixture").status == "missing"
    assert correct.download_metadata("fixture")["assets"]["dataset.sqlite"]["size"] == dataset.path.stat().st_size


def test_cross_reference_import_and_duplicate_ids_fail_atomically(tmp_path, experiment):
    manifest, samples, assays, observations = experiment
    foreign = replace(manifest.references[0], annotation_version="111")
    for changed_samples, changed_assays, changed_observations in (
        ([replace(samples[0], reference=foreign)], assays, observations),
        (samples, [replace(assays[0], reference=foreign)], observations),
        (samples, assays, [replace(observations[0], target=replace(observations[0].target, reference=foreign))]),
        (samples, assays, [observations[0], observations[0]]),
    ):
        with pytest.raises(sqlite3.IntegrityError):
            QuantitativeDataset.create(tmp_path / "bad.sqlite", manifest,
                                        samples=changed_samples, assays=changed_assays,
                                        measurements=changed_observations)
        assert not (tmp_path / "bad.sqlite").exists()


def test_materialize_import_receipt_reuse_and_failure(tmp_path, experiment):
    calls = []

    def importer(inputs, output):
        assert inputs["input.tsv"]
        calls.append(output)
        write_dataset(output, experiment)

    dataset = build_dataset(experiment[0], importer, cache_root=tmp_path / "cache",
                            download_options={"show_progress": False})
    assert dataset.counts["measurement"] == 4
    receipt = json.loads((dataset.path.parent / ".datacache-manifest.json").read_text())
    assert receipt["sources"]["input.tsv"]["sha256"] == experiment[0].sources[0].sha256
    assert receipt["outputs"]["dataset.sqlite"]["sha256"] == hashlib.sha256(dataset.path.read_bytes()).hexdigest()
    assert build_dataset(experiment[0], importer, cache_root=tmp_path / "cache").path == dataset.path
    assert len(calls) == 1

    def broken(inputs, output):
        raise RuntimeError("import failed")

    with pytest.raises(RuntimeError, match="import failed"):
        build_dataset(experiment[0], broken, cache_root=tmp_path / "cache", force=True,
                      download_options={"show_progress": False})
    assert dataset.path.exists()
    assert len(list(QuantitativeDataset(dataset.path).measurements())) == 4


def test_cli_matches_python_queries_and_original_catalogue(tmp_path, experiment, capsys):
    dataset = write_dataset(tmp_path / "source.sqlite", experiment)
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(definition(dataset)))
    args = ["datasets", "--registry", str(registry_path), "--cache-root", str(tmp_path / "cache")]
    run(args + ["list", "--species", "yeast"])
    assert json.loads(capsys.readouterr().out)[0]["dataset_id"] == "fixture"
    run(args + ["download", "fixture", "--quiet"])
    capsys.readouterr()
    run(args + ["query", "fixture", "--reference", experiment[0].references[0].key])
    assert json.loads(capsys.readouterr().out)["value"] == 0
    run(["--orf-id", "c7riboseqorf56"])
    assert json.loads(capsys.readouterr().out)["records"][0]["orf_id"] == "c7riboseqorf56"


def test_pyensembl_targets_keep_reference_and_requested_id_version(experiment):
    from pyensembl import EnsemblRelease

    for reference in experiment[0].references:
        genome = EnsemblRelease(110, species=reference.species)
        feature = SimpleNamespace(id="gene", versioned_id="gene.2", genome=genome)
        plain = ScopedTarget.from_pyensembl(feature, kind="gene", namespace="ensembl:110")
        # PyEnsembl can update the assembly metadata for a release. Scope the
        # target to the Genome's declared reference, not this synthetic fixture's build.
        assert plain.reference.taxon_id == reference.taxon_id
        assert plain.reference.species == genome.species.latin_name
        assert plain.reference.genome_build == genome.reference_name
        assert plain.reference.annotation_name == genome.annotation_name
        assert plain.reference.annotation_version == reference.annotation_version
        assert plain.identifier == "gene"
        assert ScopedTarget.from_pyensembl(feature, kind="gene", namespace="ensembl:110",
                                          versioned=True).identifier == "gene.2"
        with pytest.raises(ValueError):
            ScopedTarget.from_pyensembl(feature, kind="gene", namespace="ensembl:110",
                                       context=replace(plain.reference, annotation_version="111"))


def test_orf_target_uses_bound_context_and_rejects_modified_objects(experiment):
    from orfdata import ORFModel, ORFRecord

    context = experiment[0].references[0]
    transcript = SimpleNamespace(
        id="rna", versioned_id="rna.1", contig="1", strand="+", sequence="ATGAAATAA",
        exons=[SimpleNamespace(start=1, end=9)],
        genome=SimpleNamespace(reference_name=context.genome_build, species=context.species,
                               annotation_name=context.annotation_name, annotation_version="110"),
    )
    model = ORFModel("orf", "rna.1", "gene", "GENE", "1", "+", ((1, 9),), "MK", "ncRNA",
                     context.genome_build, context.taxon_id, context.species,
                     reference_context=context, source="synthetic:1")
    bound = ORFRecord.from_model(model).bind(transcript)
    assert ScopedTarget.from_orf(bound, namespace="test:1") == ScopedTarget(context, "orf", "orf", "test:1")
    with pytest.raises(TypeError, match="bound ReferenceORF"):
        ScopedTarget.from_orf(SimpleNamespace(reference=bound), namespace="test:1")


def test_parquet_export_retains_columns_exact_numbers_and_provenance(tmp_path, experiment):
    pq = pytest.importorskip("pyarrow.parquet")
    exact = replace(experiment[3][0], value=2 ** 53 + 1)
    dataset = write_dataset(tmp_path / "data.sqlite", experiment, measurements=[exact])
    path = tmp_path / "values.parquet"
    dataset.export_parquet(path)
    row = pq.read_table(path).to_pylist()[0]
    assert row["target_kind"] == "gene"
    assert row["unit"] == "TPM"
    assert row["dataset_id"] == "fixture"
    assert json.loads(row["measurement_json"])["value"] == exact.value
