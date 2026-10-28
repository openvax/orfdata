"""Assemble the complete canine artifact from bounded, normalized local inputs."""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

from orfdata import DatasetManifest, SourceAsset
from orfdata.artifact_archives import archive_options
from orfdata.sources.pfizer import protein_manifest
from orfdata.sources.pfizer_attribution import CanineEvidence
from orfdata.sources.pfizer_evidence import EvidenceStore
from orfdata.sources.pfizer_native import IMPORTER, file_sha256, import_native, selection_sha256

PLAN_PATH = Path(__file__).with_name("dog_result_plan.json")


def manifest_for(evidence_path, lineage_path, *, version, evidence_url, lineage_url,
                 output_bytes, scratch_bytes, compressed_evidence=None, evidence_archive=None):
    """Pin a full selection; regression slices cannot become release manifests."""
    evidence = EvidenceStore(evidence_path)
    source = evidence.metadata["source"]
    plan = json.loads(PLAN_PATH.read_text())
    if source["selection"]["plan"] != plan:
        raise ValueError("Release requires the entire pinned canine extraction selection")
    counts = evidence.metadata["reconciliation"]["counts"]
    progress = source["extraction_progress"]
    if (set(counts) != {t["table"] for t in plan["tables"]} or set(progress) != set(counts)
            or any(p["complete"] != 1 or p["rows"] != counts[name] for name, p in progress.items())):
        raise ValueError("Release requires all selected native tables to be complete")
    baseline = protein_manifest("dog")
    baseline_assets = {s.asset_id: s for s in baseline.sources}
    native_path = Path(compressed_evidence or evidence_path)
    metadata = dict(scope="Complete selected native source facts; no spectra/peak arrays",
                    integrity_scope=source["integrity_scope"])
    if evidence_archive is not None:
        if compressed_evidence is not None:
            raise ValueError("Choose either a gzip component or an archive of the normalized component")
        archive_options(evidence_archive, "evidence.sqlite", dict(
            sha256=file_sha256(evidence_path), size=Path(evidence_path).stat().st_size))
        metadata["archive"] = evidence_archive
    if compressed_evidence is not None:
        metadata.update(compression="gzip", uncompressed_size=Path(evidence_path).stat().st_size,
                        uncompressed_sha256=file_sha256(evidence_path))
    native_asset = SourceAsset("evidence.sqlite", evidence_url, file_sha256(native_path), native_path.stat().st_size,
                               "MassIVE MSV000099744 deposited analysis facts, CC0-1.0; Oziolor et al., Pfizer",
                               "https://creativecommons.org/publicdomain/zero/1.0/", True, metadata)
    lineage_path = Path(lineage_path)
    lineage_asset = SourceAsset("lineage.json", lineage_url, file_sha256(lineage_path), lineage_path.stat().st_size,
                                "Compiled factual run/channel/tissue associations; deposited metadata CC0-1.0. "
                                "Publication supplement CC BY-NC-ND 4.0; original workbook not redistributed.",
                                redistributable=True, metadata=dict(
                                    derivation="Literal factual associations with source hashes/cell locations; no original workbook",
                                    source_article="https://doi.org/10.1038/s41597-026-07512-9"))
    return DatasetManifest(
        "pfizer-atlas-dog-peptide-evidence", version,
        "Pfizer canine native peptide, PSM and protein-group evidence",
        baseline.references, (native_asset, lineage_asset, baseline_assets["search.fasta"], baseline_assets["protein.tsv"]),
        baseline.citations, ("digest_ms",),
        coverage=dict(accession="MSV000099744", proteomexchange_accession="PXD070290",
                      native_counts=counts, search_sequences=baseline.coverage["search_sequences"],
                      channels=44, tissue_pools=38, bridge_channels=4, empty_channels=2, runs=96,
                      individual_positive_dogs="unresolved tissue pools", normal_tissue_dla_presentation=False,
                      attribution="All source memberships and exact/I-L sequence-compatible candidates retained separately",
                      searched_reference="Deposited no-isoform search FASTA; not a complete canine isoform background"),
        processing=dict(importer=IMPORTER, selection_sha256=selection_sha256(plan),
                        limits=dict(output_bytes=output_bytes, scratch_bytes=scratch_bytes),
                        protein_summary=dict(dataset_id=baseline.dataset_id, version=baseline.version,
                                             source_sha256=baseline_assets["protein.tsv"].sha256,
                                             source_counts=dict(baseline.coverage["source_counts"])),
                        quantities="Native peptide-group and PSM columns, unfiltered and untransformed",
                        protein_summary_use="Same-study context only; original abundance series remain in the pinned baseline",
                        instrument_status="Publication Eclipse/deposit Lumos conflict retained with native run reports"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path, help="Complete normalized evidence.sqlite")
    parser.add_argument("lineage", type=Path, help="Checked compiled lineage JSON")
    parser.add_argument("search_fasta", type=Path)
    parser.add_argument("protein_summary", type=Path)
    parser.add_argument("output_directory", type=Path, help="New directory owned by this build")
    parser.add_argument("--version", required=True)
    parser.add_argument("--evidence-url", required=True, help="Stable intended source-component distribution URL")
    parser.add_argument("--lineage-url", required=True, help="Stable intended compiled-lineage distribution URL")
    distribution = parser.add_mutually_exclusive_group()
    distribution.add_argument("--compressed-evidence", type=Path, help="Optional gzip of the exact normalized input")
    distribution.add_argument("--evidence-archive", type=Path, help="Pinned split tar.gz archive JSON for the normalized input")
    parser.add_argument("--output-mib", required=True, type=int)
    parser.add_argument("--scratch-mib", required=True, type=int)
    args = parser.parse_args()
    started = time.monotonic()
    manifest = manifest_for(args.evidence, args.lineage, version=args.version, evidence_url=args.evidence_url,
                            lineage_url=args.lineage_url, output_bytes=args.output_mib * 1024 * 1024,
                            scratch_bytes=args.scratch_mib * 1024 * 1024, compressed_evidence=args.compressed_evidence,
                            evidence_archive=None if args.evidence_archive is None else json.loads(args.evidence_archive.read_text()))
    args.output_directory.mkdir(parents=True, exist_ok=False)
    with (args.output_directory / "manifest.json").open("x") as handle:
        json.dump(manifest.to_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    inputs = {"evidence.sqlite": args.compressed_evidence or args.evidence, "lineage.json": args.lineage,
              "search.fasta": args.search_fasta, "protein.tsv": args.protein_summary}
    dataset = import_native(inputs, args.output_directory / "dataset.sqlite", manifest,
                            checkpoint=lambda state: print(json.dumps(state, sort_keys=True), flush=True))
    evidence = CanineEvidence(dataset.path)
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    receipt = dict(dataset_id=manifest.dataset_id, version=manifest.version, dataset_sha256=file_sha256(dataset.path),
                   dataset_bytes=dataset.path.stat().st_size, counts=dataset.counts,
                   attribution=evidence.metadata["attribution_counts"], candidate_digest=evidence.metadata["candidate_digest"],
                   peak_process_rss_bytes=peak_rss, elapsed_seconds=time.monotonic() - started,
                   tool_sha256=file_sha256(__file__),
                   integrity_scope="Complete normalized selected facts; no authoritative whole-vendor-file checksum")
    with (args.output_directory / "build-receipt.json").open("x") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
