"""Package a verified canine artifact and prepare its registry metadata locally."""

import argparse
import gzip
import hashlib
import json
import os
import zlib
from pathlib import Path

from orfdata import DatasetRegistry, QuantitativeDataset
from orfdata.sources.pfizer_attribution import CanineEvidence
from orfdata.sources.pfizer_native import IMPORTER, file_sha256, selection_sha256
from tools import verify_dog_dataset
from tools.package_sqlite_archive import archive_verified

PLAN_PATH = Path(__file__).with_name("dog_result_plan.json")


class _BoundedOutput:
    def __init__(self, handle, limit):
        self.handle, self.limit, self.count = handle, limit, 0

    def write(self, block):
        if self.count + len(block) > self.limit:
            raise ValueError("Compressed artifact exceeds its explicit byte limit")
        written = self.handle.write(block)
        self.count += written
        return written


def compress_verified(source, output, *, max_bytes):
    """Deterministic-header gzip, bounded writes and streamed round-trip checks.

    This primitive alone does not establish dataset acceptance. `package`
    checks the full selection and successful source-verification receipt first.
    """
    if type(max_bytes) is not int or max_bytes < 1:
        raise ValueError("Compressed artifact byte limit must be a positive integer")
    source, output = Path(source), Path(output)
    source_hash, source_bytes = file_sha256(source), source.stat().st_size
    created = False
    try:
        with output.open("xb") as destination:
            created = True
            bounded = _BoundedOutput(destination, max_bytes)
            digest, count = hashlib.sha256(), 0
            with source.open("rb") as original, gzip.GzipFile(filename="", fileobj=bounded, mode="wb", mtime=0, compresslevel=6) as compressed:
                while block := original.read(1024 * 1024):
                    count += len(block)
                    if count > source_bytes:
                        raise ValueError("Artifact grew during compression")
                    digest.update(block)
                    compressed.write(block)
            destination.flush()
            os.fsync(destination.fileno())
        if count != source_bytes or digest.hexdigest() != source_hash:
            raise ValueError("Artifact changed during compression")
        restored_hash, restored_bytes = hashlib.sha256(), 0
        with gzip.open(output, "rb") as restored:
            while block := restored.read(min(1024 * 1024, source_bytes - restored_bytes + 1)):
                restored_bytes += len(block)
                if restored_bytes > source_bytes:
                    raise ValueError("Compressed artifact expands beyond its source byte limit")
                restored_hash.update(block)
        if restored_bytes != source_bytes or restored_hash.hexdigest() != source_hash or file_sha256(source) != source_hash:
            raise ValueError("Compressed artifact does not reproduce the unchanged original bytes")
        return dict(installed_size=source_bytes, installed_sha256=source_hash,
                    download_size=output.stat().st_size, download_sha256=file_sha256(output),
                    compression=dict(format="gzip", level=6, mtime=0, filename="", zlib=zlib.ZLIB_RUNTIME_VERSION))
    except BaseException:
        if created:
            output.unlink()
        raise


def package(dataset_path, acceptance_path, output, *, url, data_release, max_bytes, part_bytes=None):
    """Create a local archive/registry entry only after full-source verification.

    Does not upload, publish, mutate the packaged registry or replace an existing
    archive. The stored-byte and compressed-byte fingerprints remain separate.
    """
    dataset_hash = file_sha256(dataset_path)
    dataset = QuantitativeDataset(dataset_path)
    manifest = dataset.manifest
    evidence = CanineEvidence(dataset.path)
    plan = json.loads(PLAN_PATH.read_text())
    if (manifest.processing.get("importer") != IMPORTER
            or manifest.processing["selection_sha256"] != selection_sha256(plan)
            or evidence.metadata["source"]["selection"]["plan"] != plan):
        raise ValueError("Packaging requires the complete pinned canine selection")
    acceptance_path = Path(acceptance_path)
    if acceptance_path.stat().st_size > 1024 * 1024:
        raise ValueError("Canine acceptance receipt exceeds its metadata bound")
    acceptance_hash = file_sha256(acceptance_path)
    acceptance = json.loads(acceptance_path.read_text())
    if (acceptance.get("state") != "verified" or acceptance.get("dataset_sha256") != dataset_hash
            or acceptance.get("tool_sha256") != file_sha256(verify_dog_dataset.__file__)
            or acceptance.get("native_counts") != dict(manifest.coverage["native_counts"])
            or acceptance.get("candidate_digest") != evidence.metadata["candidate_digest"]
            or acceptance.get("native_group_count_mismatches") != []
            or acceptance.get("logical_cells") != dataset.counts["measurement"]):
        raise ValueError("Artifact lacks its matching successful current full-source acceptance receipt")
    if part_bytes is None and (not isinstance(url, str) or not url.startswith("https://") or not url.endswith(".gz")):
        raise ValueError("Supply the intended stable HTTPS gzip distribution URL")
    if not isinstance(data_release, str) or not data_release.strip():
        raise ValueError("Data release identity is required")
    if part_bytes is None:
        result = compress_verified(dataset.path, output, max_bytes=max_bytes)
    else:
        result = archive_verified(dataset.path, output, member="dataset.sqlite", url_prefix=url,
                                  part_bytes=part_bytes, max_bytes=max_bytes)
        result.update(download_size=result["archive"]["size"], download_sha256=result["archive"]["sha256"])
    if result["installed_sha256"] != dataset_hash or file_sha256(acceptance_path) != acceptance_hash:
        if part_bytes is None:
            Path(output).unlink()  # Only the new archive successfully created above.
        else:
            for part in result["archive"]["parts"]:
                (Path(output) / part["url"].rsplit("/", 1)[-1]).unlink()
        raise ValueError("Verified inputs changed while preparing publication")
    distribution = dict(
        access_terms="Native deposited analysis facts: MassIVE CC0-1.0, Oziolor et al./Pfizer. "
                     "Full searched FASTA and UniProt annotations: CC-BY-4.0, attribution to UniProt Consortium. "
                     "Compiled factual lineage cites the CC BY-NC-ND 4.0 supplement; original workbook excluded. "
                     "ORFdata transformation code: Apache-2.0.",
        redistributable=True, download_size=result["download_size"], download_sha256=result["download_sha256"],
        installed_size=result["installed_size"], data_release=data_release,
        acceptance_sha256=acceptance_hash)
    asset = dict(size=result["installed_size"], sha256=result["installed_sha256"])
    if part_bytes is None:
        asset.update(url=url, decompress=True)
    release = dict(manifest=manifest.to_dict(), distribution=distribution, assets={"dataset.sqlite": asset})
    if part_bytes is not None:
        release["archive"] = result["archive"]
    definitions = {manifest.dataset_id: dict(default_version=manifest.version, versions={manifest.version: release})}
    DatasetRegistry(definitions)  # Metadata validation only; no acquisition.
    return dict(registry=definitions, compression=result, tool_sha256=file_sha256(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("acceptance", type=Path)
    parser.add_argument("directory", type=Path, help="New directory for reviewed publication artifacts")
    parser.add_argument("--url", required=True, help="HTTPS gzip URL, or release URL prefix ending in / with --part-mib")
    parser.add_argument("--data-release", required=True)
    parser.add_argument("--max-compressed-mib", required=True, type=int)
    parser.add_argument("--part-mib", type=int, help="Create tar.gz parts smaller than 2048 MiB, retaining one logical dataset.sqlite")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=False)
    output = args.directory / ("archive" if args.part_mib is not None else "dataset.sqlite.gz")
    result = package(args.dataset, args.acceptance, output, url=args.url,
                     data_release=args.data_release, max_bytes=args.max_compressed_mib * 1024 * 1024,
                     part_bytes=None if args.part_mib is None else args.part_mib * 1024 * 1024)
    for name, value in (("registry-entry.json", result["registry"]),
                        ("compression-receipt.json", {k: v for k, v in result.items() if k != "registry"})):
        with (args.directory / name).open("x") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(result["compression"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
