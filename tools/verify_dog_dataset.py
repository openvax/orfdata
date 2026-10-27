"""Verify a canine artifact against every selected native row and pinned inputs.

This is an explicit local release audit. It never opens the remote vendor file.
All native rows, membership indexes and quantitative layouts are checked; public
measurement queries and independent exhaustive FASTA searches use a declared,
deterministically spread sample. A receipt distinguishes these scopes.
"""

import argparse
import hashlib
import json
import math
import resource
import sqlite3
import struct
import sys
import time
from collections import Counter
from contextlib import closing
from functools import lru_cache
from itertools import zip_longest
from pathlib import Path

from orfdata import QuantitativeDataset
from orfdata.dataset import _json
from orfdata.search_reference import SearchReference
from orfdata.sources.pfizer import _audit
from orfdata.sources.pfizer_attribution import ENDPOINTS
from orfdata.sources.pfizer_dataset import METADATA_TABLES
from orfdata.sources.pfizer_evidence import EvidenceStore, record_keys, source_key
from orfdata.sources.pfizer_lineage import lineage_records
from orfdata.sources.pfizer_native import IMPORTER, file_sha256, search_references, selection_sha256
from orfdata.sources.pfizer_vendor import DistributionDecoder, NativeQuantification
from tools.export_sqlite_rows import canonical_row, quoted

PLAN_PATH = Path(__file__).with_name("dog_result_plan.json")


def _spread(total, limit):
    count = min(total, limit)
    return {i * (total - 1) // max(1, count - 1) for i in range(count)}


def _native_rows(db, name, columns):
    db.row_factory = sqlite3.Row
    return db.execute('SELECT _source_rowid,' + ','.join(quoted(c) for c in columns) +
                      ' FROM ' + quoted(name) + ' ORDER BY _source_rowid')


def _candidate_sample(evidence, searched, sequence):
    def expected():
        found = False
        for protein in sorted(searched.entries(), key=lambda p: p.target.identifier):
            haystack, needle = protein.sequence.replace("I", "L"), sequence.replace("I", "L")
            start = 0
            while (start := haystack.find(needle, start)) >= 0:
                end = start + len(sequence)
                found = True
                yield dict(peptide=sequence, accession=protein.target.identifier, start=start, end=end,
                           status="candidate", basis="exact_residues" if protein.sequence[start:end] == sequence
                           else "isoleucine_leucine_equivalent")
                start += 1
        if not found:
            yield dict(peptide=sequence, accession=None, start=None, end=None,
                       status="no_search_match", basis="exact_or_il_equivalent")

    count = 0
    for original, retained in zip_longest(expected(), evidence.candidates(sequence)):
        if original != retained:
            raise ValueError("Independent full-FASTA search differs for sampled peptide: " + sequence)
        count += original["status"] == "candidate"
    return count


def _check_numeric_bytes(converter, table, row, column, observed):
    """Independent struct reads; do not reuse NumericVector or its cell decoder."""
    _, distribution, boxes = converter.decoder.distribution(table, column)
    format_string = {1: "<i", 2: "<q", 3: "<d"}[distribution["CustomDataType"]]
    width = struct.calcsize(format_string)
    blob = row[column]
    if blob is not None and len(blob) != len(boxes) * (width + 1):
        raise ValueError("Original native numeric byte count differs from its definition")
    values = {m.provenance["native_distribution"]["box_id"]: m for m in observed.values()}
    if len(values) != len(boxes):
        raise ValueError("Public numeric cells do not cover all original boxes")
    for position, box in enumerate(boxes):
        number, present = None, False
        if blob is not None:
            offset = position * (width + 1)
            number = struct.unpack_from(format_string, blob, offset)[0]
            flag = blob[offset + width]
            if flag not in (0, 1):
                raise ValueError("Original native numeric presence byte is invalid")
            present = bool(flag)
        nonfinite = isinstance(number, float) and not math.isfinite(number)
        expected_value = number if present and not nonfinite else None
        expected_status = "measured" if expected_value is not None else "missing"
        measurement = values[box["BoxID"]]
        original_number = ("NaN" if math.isnan(number) else "Infinity" if number > 0 else "-Infinity") if nonfinite else number
        provenance = measurement.provenance["native_distribution"]
        native_fields = dict(native_present=present, native_value=original_number,
                             status="missing" if not present else "nonfinite" if nonfinite else "measured",
                             missing_reason="null_vector" if blob is None else "null_element" if not present else None)
        bits_match = (not isinstance(expected_value, float) or measurement.value is not None
                      and struct.pack("<d", measurement.value) == struct.pack("<d", expected_value))
        if (measurement.value != expected_value or measurement.status != expected_status or not bits_match
                or any(provenance.get(k) != v for k, v in native_fields.items())):
            raise ValueError("Independent native numeric bytes differ from public value/missingness")
    return len(boxes)


def verify(dataset_path, native_path, fasta_path, summary_path, *, lineage_path, sample_records=100, checkpoint=None):
    if type(sample_records) is not int or not 1 <= sample_records <= 1000:
        raise ValueError("Audit sample must be between 1 and 1000 records per quantitative table")
    dataset = QuantitativeDataset(dataset_path)
    manifest = dataset.manifest
    if manifest.processing.get("importer") != IMPORTER:
        raise ValueError("Expected a completed canine evidence generation")
    evidence = dataset.source_evidence()
    plan = json.loads(PLAN_PATH.read_text())
    native_path = Path(native_path).resolve()
    native_hash, dataset_hash = file_sha256(native_path), file_sha256(dataset.path)
    if native_hash != evidence.metadata["source"]["source_file_sha256"]:
        raise ValueError("Original selected native file differs from the normalized source pin")
    assets = {a.asset_id: a for a in manifest.sources}
    if Path(summary_path).stat().st_size != assets["protein.tsv"].size or file_sha256(summary_path) != assets["protein.tsv"].sha256:
        raise ValueError("Original protein summary differs from the manifest")
    if Path(lineage_path).stat().st_size != assets["lineage.json"].size or file_sha256(lineage_path) != assets["lineage.json"].sha256:
        raise ValueError("Checked lineage input differs from the manifest")
    lineage = json.loads(Path(lineage_path).read_text())
    spec, = search_references(manifest)
    searched = SearchReference(fasta_path, spec)
    tables, selected, digests, sequences = {}, [], {}, set()
    counters = Counter(native_rows=0, links=0, unresolved_links=0, projection_vectors=0, logical_cells=0)
    with closing(sqlite3.connect(native_path.as_uri() + "?mode=ro", uri=True)) as native, closing(evidence.connect()) as db:
        for connection in (native, db):
            connection.execute("PRAGMA query_only=ON")
            connection.execute("PRAGMA cache_size=-8192")
            connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 16 * 1024 * 1024)
            connection.execute("BEGIN")
        configuration = json.loads(native.execute("SELECT value FROM _export_metadata WHERE key='configuration'").fetchone()[0])
        if (configuration["plan"] != plan or configuration != evidence.metadata["source"]["selection"]
                or selection_sha256(plan) != manifest.processing["selection_sha256"]):
            raise ValueError("Release audit requires the complete pinned native selection")
        progress = {r[0]: dict(rows=r[1], digest=r[2], complete=r[3], missing_keys=r[4]) for r in native.execute(
            "SELECT name,rows,digest,complete,missing_keys FROM _export_progress")}
        if progress != evidence.metadata["source"]["extraction_progress"] or any(p["complete"] != 1 for p in progress.values()):
            raise ValueError("Native extraction progress differs or is incomplete")
        if set(evidence.tables) != set(progress) or set(progress) != {t["table"] for t in plan["tables"]}:
            raise ValueError("Selected native table coverage differs")
        storage = json.loads(db.execute("SELECT value FROM metadata WHERE name='measurement_storage'").fetchone()[0])
        if storage["lineage"] != lineage:
            raise ValueError("Retained lineage differs from the checked lineage input")
        samples, assays = lineage_records(lineage, spec.reference)
        if ({s.sample_id: s for s in dataset.samples()} != {s.sample_id: s for s in samples}
                or {a.assay_id: a for a in dataset.assays()} != {a.assay_id: a for a in assays}):
            raise ValueError("Common samples/assays differ from checked lineage")
        metadata = [[dict(row) for row in _native_rows(native, name, evidence.tables[name][1])] for name in METADATA_TABLES]
        converter = NativeQuantification(DistributionDecoder(*metadata), storage["lineage"], spec.reference,
                                         manifest.dataset_id, manifest.version, source_asset="evidence.sqlite")
        layouts = {index: tuple(db.execute(
            "SELECT x.position,x.box_id,a.id,a.sample_id FROM peptide_axis x JOIN assay a ON a.i=x.assay_i "
            "WHERE x.layout_i=? ORDER BY x.position", (index,))) for (index,) in db.execute("SELECT i FROM peptide_layout")}

        @lru_cache(maxsize=8192)
        def endpoint(name, values):
            found = db.execute("SELECT i FROM evidence_record WHERE table_i=? AND source_key=?",
                               (evidence.tables[name][0], source_key(values))).fetchmany(2)
            if len(found) > 1:
                raise ValueError("Native endpoint identity is not unique")
            return found[0][0] if found else None

        for item in plan["tables"]:
            name, columns = item["table"], item["columns"]
            index, stored_columns = evidence.tables[name]
            if tuple(columns) != stored_columns or db.execute("SELECT role FROM evidence_table WHERE i=?", (index,)).fetchone()[0] != item["role"]:
                raise ValueError("Native columns or partition differ: " + name)
            source_rows = _native_rows(native, name, columns)
            artifact_rows = db.execute("SELECT i,source_rowid,attributes,source_key,sequence FROM evidence_record "
                                       "WHERE table_i=? ORDER BY source_rowid", (index,))
            count, digest = 0, bytes(32)
            sample_positions = _spread(progress[name]["rows"], sample_records) if name in converter.COLUMNS else set()
            for original, record in zip_longest(source_rows, artifact_rows):
                if original is None or record is None:
                    raise ValueError("Artifact/native row count differs: " + name)
                row, restored = dict(original), EvidenceStore.restore(db, record[:3], columns=columns)
                encoded = canonical_row(tuple(original))
                if set(row) != set(restored) or encoded != canonical_row(tuple(restored[k] for k in original.keys())):
                    raise ValueError(f"Artifact/native typed fields differ: {name} row {row['_source_rowid']}")
                if record[3] != source_key(tuple(row[k] for k in record_keys(name))) or record[4] != row.get("Sequence"):
                    raise ValueError("Native identity/sequence index differs: " + name)
                digest = hashlib.sha256(digest + encoded).digest()
                if name in ENDPOINTS:
                    expected = tuple(endpoint(t, tuple(row[t + k] for k in record_keys(t))) for t in ENDPOINTS[name])
                    found = db.execute("SELECT left_i,right_i FROM evidence_link WHERE record_i=?", (record[0],)).fetchone()
                    if found != expected:
                        raise ValueError("Indexed membership differs from original native endpoint IDs")
                    counters["links"] += 1
                    counters["unresolved_links"] += int(None in expected)
                if name in ("TargetProteins", "DecoyProteins"):
                    protein = searched.get(row["Accession"]) if name == "TargetProteins" else None
                    expected = (protein.target.identifier if protein else None,
                                "decoy_partition" if name == "DecoyProteins" else "not_in_search" if protein is None
                                else "exact_sequence" if protein.sequence == row["Sequence"] else "sequence_mismatch")
                    found = db.execute("SELECT accession,status FROM native_protein_search WHERE record_i=?", (record[0],)).fetchone()
                    if found != expected:
                        raise ValueError("Native protein/search identity index differs")
                    counters["native_protein_search_rows"] += 1
                if name in converter.COLUMNS:
                    projections = {p[0]: p[1:] for p in db.execute(
                        "SELECT p.column_name,p.layout_i,p.prefix,p.unit,t.identifier,t.reference_key,t.kind,t.namespace "
                        "FROM peptide_projection p JOIN target t ON t.i=p.target_i WHERE p.record_i=?", (record[0],))}
                    if set(projections) != set(converter.COLUMNS[name]):
                        raise ValueError("Native quantity projection coverage differs")
                    source_id = ":".join(str(row[k]) for k in record_keys(name))
                    for column, projected in projections.items():
                        axes = converter.peptide_axes(name, row, column)
                        expected_axes = tuple((i, axis["box_id"], axis["assay_id"], axis["sample_id"]) for i, axis in enumerate(axes))
                        expected_fields = (f"{name}:{source_id}:{column}:", f"PD2.4 {name}.{column}", row["Sequence"],
                                           spec.reference.key, "peptide", "amino-acid-sequence")
                        if projected[1:] != expected_fields or layouts[projected[0]] != expected_axes:
                            raise ValueError("Native quantitative target, unit or sample axes differ")
                        counters["projection_vectors"] += 1
                        counters["logical_cells"] += len(axes)
                    if count in sample_positions:
                        selected.append((name, row["_source_rowid"]))
                        sequences.add(row["Sequence"])
                count += 1
            tables[name], digests[name] = count, digest.hex()
            if count != progress[name]["rows"] or digest.hex() != evidence.metadata["reconciliation"]["native_rowid_order_sha256"][name]:
                raise ValueError("Native row count/digest differs from normalization: " + name)
            counters["native_rows"] += count
            if checkpoint is not None:
                checkpoint(dict(table=name, rows=count, native_sha256=digest.hex()))
        if tables != dict(manifest.coverage["native_counts"]):
            raise ValueError("Artifact table counts differ from the release manifest")
        for table, key in (("evidence_link", "links"), ("peptide_projection", "projection_vectors"),
                           ("native_protein_search", "native_protein_search_rows")):
            if db.execute("SELECT count(*) FROM " + table).fetchone()[0] != counters[key]:
                raise ValueError("Unexpected extra indexed records in " + table)
        if counters["logical_cells"] != dataset.counts["measurement"]:
            raise ValueError("Logical measurement count differs from original native axes")
        missing = native.execute("SELECT name,source_key,reason FROM _export_missing ORDER BY name,source_key")
        retained = db.execute("SELECT table_name,source_key,reason FROM evidence_missing ORDER BY table_name,source_key")
        for a, b in zip_longest(missing, retained):
            if a is None or b is None or tuple(a) != tuple(b):
                raise ValueError("Native missing-acquisition reports differ")
        # Every original search entry is retained, including unreported proteins.
        found_proteins = 0
        for accession, serialized in db.execute("SELECT accession,record FROM search_protein"):
            protein = searched.get(accession)
            expected = protein.to_dict() if protein else {}
            expected.pop("target", None)
            if json.loads(serialized) != expected:
                raise ValueError("Retained search protein differs from the complete FASTA")
            found_proteins += 1
        if found_proteins != len(searched):
            raise ValueError("Incomplete searched-reference entries")
        # Full candidate-index completion/digest, independently of sampled searches.
        peptide_tables = tuple(evidence.tables[name][0] for name in converter.COLUMNS)
        for table_i in peptide_tables:
            if db.execute("SELECT 1 FROM evidence_record r LEFT JOIN search_peptide p ON p.sequence=r.sequence "
                          "WHERE r.table_i=? AND p.sequence IS NULL LIMIT 1", (table_i,)).fetchone():
                raise ValueError("A native target peptide is absent from the candidate index")
        if db.execute("SELECT 1 FROM search_peptide p WHERE NOT EXISTS (SELECT 1 FROM evidence_record r "
                      "WHERE r.table_i IN (?,?) AND r.sequence=p.sequence) LIMIT 1", peptide_tables).fetchone():
            raise ValueError("Candidate index contains a sequence absent from native target records")
        candidate_digest, candidate_counts = hashlib.sha256(), Counter()
        for sequence, status in db.execute("SELECT sequence,status FROM search_peptide ORDER BY sequence"):
            if status not in ("matched", "no_search_match"):
                raise ValueError("Candidate search is incomplete")
            candidate_counts["peptide_sequences"] += 1
            hits = db.execute("SELECT accession,start,end,basis FROM search_peptide_hit WHERE sequence=? ORDER BY accession,start,end", (sequence,))
            present = False
            for hit in hits:
                present = True
                candidate_counts["exact_placements" if hit[3] == "exact_residues" else "il_equivalent_placements"] += 1
                encoded = _json((sequence, status, *hit)).encode()
                candidate_digest.update(len(encoded).to_bytes(8, "little") + encoded)
            if present != (status == "matched"):
                raise ValueError("Candidate status and retained placements disagree")
            if not present:
                candidate_counts["no_search_match"] += 1
                encoded = _json((sequence, status, None, None, None, None)).encode()
                candidate_digest.update(len(encoded).to_bytes(8, "little") + encoded)
        if candidate_digest.hexdigest() != evidence.metadata["candidate_digest"] or any(
                candidate_counts[k] != evidence.metadata["attribution_counts"][k] for k in
                ("peptide_sequences", "no_search_match", "exact_placements", "il_equivalent_placements")):
            raise ValueError("Candidate index digest/counts differ from the completed build")
        sampled_cells, independent_numeric_cells = 0, 0
        for name, rowid in selected:
            row = dict(native.execute("SELECT * FROM " + quoted(name) + " WHERE _source_rowid=?", (rowid,)).fetchone())
            for column in converter.COLUMNS[name]:
                expected = {m.measurement_id: m for m in converter.peptide_abundances(name, row, column)}
                first = next(iter(expected.values()))
                observed = {m.measurement_id: m for m in dataset.measurements(target=first.target, unit=first.unit)
                            if m.measurement_id in expected}
                if observed != expected:
                    raise ValueError("Public common measurements differ from original native values")
                sampled_cells += len(expected)
                independent_numeric_cells += _check_numeric_bytes(converter, name, row, column, observed)
        _, groups, summary_counts = _audit(summary_path)
        if summary_counts != dict(manifest.processing["protein_summary"]["source_counts"]):
            raise ValueError("Original summary counts differ from pinned baseline")
        if dict(db.execute("SELECT group_id,accession FROM summary_representative")) != groups:
            raise ValueError("Summary group/representative identities differ")
        native_accessions = {}
        for accession, identifier in native.execute("SELECT Accession,UniqueSequenceID FROM TargetProteins"):
            native_accessions.setdefault(accession, []).append(identifier)
        summary_matches = 0
        for group, accession in groups.items():
            expected = {endpoint("TargetProteins", (identifier,)) for identifier in native_accessions.get(accession, ())}
            observed = {row[0] for row in db.execute("SELECT record_i FROM summary_native WHERE group_id=?", (group,))}
            if observed != expected:
                raise ValueError("Summary representative/native identity index differs")
            summary_matches += len(expected)
        if db.execute("SELECT count(*) FROM summary_native").fetchone()[0] != summary_matches:
            raise ValueError("Unexpected extra summary/native identity records")
        group_sizes, group_mismatches = Counter(), []
        link_table = evidence.tables["TargetProteinGroupsTargetProteins"][0]
        for group in native.execute("SELECT ProteinGroupID,ProteinCount FROM TargetProteinGroups"):
            group_record = endpoint("TargetProteinGroups", (group[0],))
            count = db.execute("SELECT count(*) FROM evidence_link l JOIN evidence_record r ON r.i=l.record_i "
                               "WHERE l.left_i=? AND r.table_i=?", (group_record, link_table)).fetchone()[0]
            group_sizes[count] += 1
            if count != group[1]:
                group_mismatches.append(dict(group_id=group[0], declared_members=group[1], membership_rows=count))
        if db.execute("PRAGMA foreign_key_check").fetchone() is not None or db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Artifact SQLite integrity check failed")
    candidate_placements = sum(_candidate_sample(evidence, searched, sequence) for sequence in sorted(sequences))
    if (file_sha256(native_path) != native_hash or file_sha256(dataset.path) != dataset_hash
            or file_sha256(fasta_path) != spec.source.sha256 or file_sha256(summary_path) != assets["protein.tsv"].sha256
            or file_sha256(lineage_path) != assets["lineage.json"].sha256):
        raise ValueError("Artifact or pinned local input changed during verification")
    return dict(state="verified", dataset_sha256=dataset_hash, native_file_sha256=native_hash,
                lineage_sha256=assets["lineage.json"].sha256, search_fasta_sha256=spec.source.sha256,
                protein_summary_sha256=assets["protein.tsv"].sha256,
                native_counts=tables, native_rowid_order_sha256=digests, **dict(counters),
                search_proteins_compared=found_proteins, candidate_digest=candidate_digest.hexdigest(),
                summary_source_counts=summary_counts, sampled_quantitative_records=len(selected),
                summary_native_matches=summary_matches, native_group_size_histogram=dict(group_sizes),
                native_group_count_mismatches=group_mismatches,
                sampled_common_measurements=sampled_cells, independently_searched_peptides=len(sequences),
                independently_decoded_native_cells=independent_numeric_cells,
                independently_checked_candidate_placements=candidate_placements, sample_records_per_table=sample_records,
                verification_scope="All selected typed native rows, membership indexes, quantity layouts and search entries; "
                                   "public measurement expansion and exhaustive full-FASTA search checked on the declared spread sample",
                source_integrity_scope="Selected local source facts; no authoritative checksum of the whole vendor database")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("native_export", type=Path)
    parser.add_argument("lineage", type=Path)
    parser.add_argument("search_fasta", type=Path)
    parser.add_argument("protein_summary", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--sample-records", type=int, default=100)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    started = time.monotonic()
    result = verify(args.dataset, args.native_export, args.search_fasta, args.protein_summary,
                    lineage_path=args.lineage, sample_records=args.sample_records,
                    checkpoint=lambda state: print(json.dumps(state, sort_keys=True), flush=True))
    result.update(tool_sha256=file_sha256(__file__), elapsed_seconds=time.monotonic() - started,
                  peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024))
    with args.receipt.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
