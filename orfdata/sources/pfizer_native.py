"""Build the canine evidence artifact from pinned local, normalized inputs.

The explicit build joins source facts and complete searched-FASTA candidates.
It never reads the vendor database or copies protein-summary abundances into
peptide observations. Ordinary dataset queries require none of the build tools.
"""

import gzip
import hashlib
import json
import math
import shutil
import sqlite3
import tempfile
from collections import Counter
from contextlib import closing
from pathlib import Path

from ..dataset import QuantitativeDataset, _json
from ..search_reference import SearchReference, SearchReferenceSpec
from .pfizer import SAMPLE_COLUMNS, VALUE_COLUMNS, _rows, _value, protein_sample_id
from .pfizer_attribution import CanineEvidence, build_attribution
from .pfizer_dataset import copy_evidence, install_projection
from .pfizer_evidence import EvidenceStore
from .pfizer_lineage import lineage_records

IMPORTER = "pfizer-canine-evidence-v1"
ASSETS = {"evidence.sqlite", "lineage.json", "search.fasta", "protein.tsv"}
BUILD_MODULES = ("pfizer_native.py", "pfizer_dataset.py", "pfizer_attribution.py", "pfizer_evidence.py",
                 "pfizer_vendor.py", "pfizer_lineage.py", "pfizer.py", "../dataset.py", "../numeric_vector.py",
                 "../search_reference.py")


def file_sha256(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def selection_sha256(plan):
    """Canonical selection hash, independent of JSON whitespace."""
    return hashlib.sha256(_json(plan).encode()).hexdigest()


def search_references(manifest):
    if manifest.processing.get("importer") != IMPORTER or len(manifest.references) != 1:
        raise ValueError("Not a single-reference canine evidence manifest")
    reference, = manifest.references
    if reference.taxon_id != 9615:
        raise ValueError("Canine evidence requires a canine reference")
    asset = next(a for a in manifest.sources if a.asset_id == "search.fasta")
    return (SearchReferenceSpec(manifest.dataset_id, manifest.version, reference, asset,
                                "uniprot", "uniprot-fasta", manifest.coverage["search_sequences"]),)


def _verify_inputs(inputs, assets):
    for name, asset in assets.items():
        path = Path(inputs[name])
        if path.stat().st_size != asset.size or file_sha256(path) != asset.sha256:
            raise ValueError("Canine input size/SHA-256 disagrees with manifest: " + name)


def _evidence_input(inputs, asset, directory):
    if not asset.metadata.get("compression"):
        return Path(inputs[asset.asset_id])
    if asset.metadata["compression"] != "gzip":
        raise ValueError("Unsupported canine evidence compression")
    expected = asset.metadata["uncompressed_size"]
    path = Path(directory) / "evidence.sqlite"
    count, digest = 0, hashlib.sha256()
    with gzip.open(inputs[asset.asset_id], "rb") as source, path.open("xb") as output:
        while block := source.read(min(1024 * 1024, expected - count + 1)):
            count += len(block)
            if count > expected:
                raise ValueError("Expanded canine evidence exceeds its pinned byte limit")
            digest.update(block)
            output.write(block)
    if count != expected or digest.hexdigest() != asset.metadata["uncompressed_sha256"]:
        raise ValueError("Expanded canine evidence size/SHA-256 disagrees with manifest")
    return path


def _validate_component(evidence, manifest, lineage):
    source = evidence.metadata["source"]
    selection = source["selection"]
    if selection_sha256(selection["plan"]) != manifest.processing["selection_sha256"]:
        raise ValueError("Native selection differs from the manifest's complete source selection")
    expected = dict(manifest.coverage["native_counts"])
    progress = source["extraction_progress"]
    if (set(expected) != set(evidence.tables) or set(progress) != set(expected)
            or evidence.metadata["reconciliation"]["counts"] != expected
            or any(p["complete"] != 1 or p["rows"] != expected[name] for name, p in progress.items())):
        raise ValueError("Native extraction completeness/counts disagree with the manifest")
    # A lineage receipt can predate an additive table selection. Only its native
    # source identity must match; the selection is independently pinned above.
    if selection["identity"] != lineage["native_export_configuration"]["identity"]:
        raise ValueError("Lineage and evidence identify different native source databases")
    if lineage["counts"] != dict(channels=44, runs=96, tissue_pool=38, bridge=4, empty=2):
        raise ValueError("Canine lineage must retain all channels, controls and 96 runs")
    files = {(r["WorkflowID"], r["FileID"]): {k: v for k, v in r.items() if k != "_source_rowid"}
             for r in evidence.records("WorkflowInputFiles") if r["FileType"] == "XcaliburRawfile"}
    if len(files) != len(lineage["runs"]) or len(files) != 96:
        raise ValueError("Native raw run coverage differs from lineage")
    for run in lineage["runs"]:
        if files.get((run["workflow_id"], run["file_id"])) != run["vendor"]:
            raise ValueError("Native run record differs from checked lineage")
    with closing(evidence.connect()) as db:
        for name, (index, _) in evidence.tables.items():
            actual = db.execute("SELECT count(*) FROM evidence_record WHERE table_i=?", (index,)).fetchone()[0]
            if actual != expected[name]:
                raise ValueError("Native evidence row count differs from manifest: " + name)


def _summary_context(db, path, manifest, lineage):
    """Reconcile all summary rows and retain every native representative match.

    Summary measurements live in the independently pinned baseline dataset.
    This table records their representative's native identities, never unique
    protein attribution. A disk-bounded key table detects duplicate source rows.
    """
    baseline = manifest.to_dict()["processing"]["protein_summary"]
    asset = next(a for a in manifest.sources if a.asset_id == "protein.tsv")
    if baseline["source_sha256"] != asset.sha256 or lineage["source_sha256"]["protein_summary"] != asset.sha256:
        raise ValueError("Protein summary fingerprint differs from baseline or lineage")
    samples = {c["summary_sample_id"]: c["summary_sample"] for c in lineage["channels"] if c["summary_sample_id"]}
    db.executescript("""
        CREATE TABLE summary_representative(group_id TEXT PRIMARY KEY, accession TEXT NOT NULL) WITHOUT ROWID;
        CREATE TABLE summary_native(group_id TEXT NOT NULL REFERENCES summary_representative(group_id),
            record_i INTEGER NOT NULL REFERENCES evidence_record(i), PRIMARY KEY(group_id,record_i)) WITHOUT ROWID;
        CREATE TABLE summary_seen(group_id TEXT NOT NULL, sample_id TEXT NOT NULL,
            PRIMARY KEY(group_id,sample_id)) WITHOUT ROWID;
    """)
    counts, observed_samples = Counter(), set()
    for _, row in _rows(path):
        sample = protein_sample_id(row)
        if samples.get(sample) != {key: row[key] for key in SAMPLE_COLUMNS}:
            raise ValueError("Protein summary sample differs from checked lineage")
        observed_samples.add(sample)
        group, accession = row["ProteinGroup"], row["Protein"]
        previous = db.execute("SELECT accession FROM summary_representative WHERE group_id=?", (group,)).fetchone()
        if previous and previous[0] != accession:
            raise ValueError("Conflicting summary protein-group representatives")
        db.execute("INSERT OR IGNORE INTO summary_representative VALUES (?,?)", (group, accession))
        try:
            db.execute("INSERT INTO summary_seen VALUES (?,?)", (group, sample))
        except sqlite3.IntegrityError:
            raise ValueError("Duplicate summary protein-group/sample row") from None
        counts["rows"] += 1
        for column in VALUE_COLUMNS:
            value = _value(row[column])
            counts[column + ":" + ("missing" if value is None else "zero" if value == 0 else "positive")] += 1
        if not math.isfinite(float(row["sumIBAQ"])) or float(row["sumIBAQ"]) <= 0:
            raise ValueError("Invalid source normalization sum")
        if counts["rows"] % 5000 == 0:
            db.commit()
    groups = db.execute("SELECT count(*) FROM summary_representative").fetchone()[0]
    counts.update(groups=groups, assays=len(observed_samples), absent_matrix_cells=groups * len(samples) - counts["rows"])
    if observed_samples != set(samples) or dict(counts) != dict(baseline["source_counts"]):
        raise ValueError("Protein summary source counts differ from pinned baseline")
    db.execute("INSERT INTO summary_native SELECT s.group_id,p.record_i FROM summary_representative s "
               "JOIN native_protein_search p ON p.accession=s.accession "
               "JOIN evidence_record r ON r.i=p.record_i JOIN evidence_table t ON t.i=r.table_i "
               "WHERE t.name='TargetProteins'")
    db.execute("DROP TABLE summary_seen")
    db.execute("INSERT INTO evidence_metadata VALUES ('protein_summary',?)", (_json(dict(baseline)),))
    db.commit()
    return dict(counts)


def import_native(inputs, output, manifest, *, checkpoint=None):
    """Build a new complete artifact, with byte limits and owned-file cleanup.

    Four inputs are explicitly pinned: normalized evidence, compiled lineage,
    the entire searched FASTA and the original protein summary. Optional gzip
    evidence also pins its expanded size/hash. No network calls are made here.
    """
    spec, = search_references(manifest)
    assets = {asset.asset_id: asset for asset in manifest.sources}
    if set(assets) != ASSETS or set(inputs) != ASSETS:
        raise ValueError("Canine evidence requires exactly its four pinned input assets")
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    limits = manifest.processing["limits"]
    cap, scratch = limits["output_bytes"], limits["scratch_bytes"]
    if any(type(n) is not int or n < 4096 for n in (cap, scratch)):
        raise ValueError("Canine output and scratch byte limits must be positive page-sized integers")
    asset = assets["evidence.sqlite"]
    expanded = asset.metadata.get("uncompressed_size", 0) if asset.metadata.get("compression") else 0
    if type(expanded) is not int or expanded < 0:
        raise ValueError("Expanded evidence must have a nonnegative pinned size")
    # Reserve a full output, rollback journal and index/temp allowance, plus
    # all pinned inputs and the optional expanded component. No hidden download.
    working = expanded + 3 * cap + 1024 * 1024
    if sum(a.size for a in assets.values()) + working > scratch:
        raise ValueError("Canine scratch limit cannot hold pinned inputs and bounded output")
    if shutil.disk_usage(output.parent).free < working:
        raise ValueError("Insufficient free disk for the bounded canine build")
    if assets["lineage.json"].size > 1024 * 1024:
        raise ValueError("Compiled canine lineage exceeds its 1 MiB metadata limit")
    _verify_inputs(inputs, assets)
    lineage = json.loads(Path(inputs["lineage.json"]).read_text())
    created, db = False, None
    try:
        with tempfile.TemporaryDirectory(prefix="canine-evidence-", dir=output.parent) as temporary:
            native_path = _evidence_input(inputs, asset, temporary)
            evidence = EvidenceStore(native_path)
            _validate_component(evidence, manifest, lineage)
            searched = SearchReference(inputs["search.fasta"], spec)
            samples, assays = lineage_records(lineage, spec.reference)
            QuantitativeDataset.create(output, manifest, samples=samples, assays=assays, measurements=(), max_bytes=cap)
            created = True
            db = sqlite3.connect(output)
            db.execute(f"PRAGMA max_page_count={cap // 4096}")
            db.execute("PRAGMA foreign_keys=ON")
            db.execute("PRAGMA cache_size=-8192")
            db.execute("PRAGMA temp_store=MEMORY")
            db.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
            copy_evidence(db, native_path, checkpoint=checkpoint)
            attribution = build_attribution(db, searched, checkpoint=checkpoint)
            summary = _summary_context(db, inputs["protein.tsv"], manifest, lineage)
            projection = install_projection(db, manifest, lineage, source_asset=asset.asset_id, checkpoint=checkpoint)
            if db.execute("PRAGMA foreign_key_check").fetchone() is not None or db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("Canine artifact failed SQLite integrity checks")
            _verify_inputs(inputs, assets)
            receipt = dict(state="complete", importer=IMPORTER, attribution=attribution, summary=summary,
                           projection=projection, selection_sha256=manifest.processing["selection_sha256"],
                           input_sha256={name: a.sha256 for name, a in assets.items()},
                           software_sha256={name: file_sha256(Path(__file__).parent / name) for name in BUILD_MODULES},
                           limits=dict(limits), sqlite_version=sqlite3.sqlite_version)
            db.execute("INSERT INTO metadata VALUES ('canine_build',?)", (_json(receipt),))
            db.commit()
            db.close()
            db = None
        dataset = QuantitativeDataset(output)
        CanineEvidence(output)
        return dataset
    except BaseException:
        if db is not None:
            db.close()
            db = None
        if created:
            output.unlink()
        raise
    finally:
        if db is not None:
            db.close()
