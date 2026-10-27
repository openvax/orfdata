"""Project native canine indexes into the source-independent public query models.

No extra artifact copy, inference, sequence search or source acquisition occurs.
Reported assignments retain their native record identity. Complete searched-FASTA
candidates are a separate mapping with an explicit sequence-compatibility basis.
"""

import hashlib
import json
from contextlib import closing
from functools import lru_cache

from ..measurements import TargetMapping
from ..pairing import TargetRelation
from ..reference import ScopedTarget
from .pfizer_attribution import CanineEvidence
from .pfizer_evidence import source_key

_REPORTS = {
    "TargetPeptideGroups": "TargetPeptideGroupsTargetProteins",
    "TargetPsms": "TargetPsmsTargetProteins",
    "TargetProteinGroups": "TargetProteinGroupsTargetProteins",
}


def _native_key(identifier):
    try:
        value = int(identifier)
    except ValueError:
        return None
    return source_key((value,)) if str(value) == identifier else None


def _has_attribution(connection):
    return connection.execute("SELECT 1 FROM evidence_metadata WHERE name='attribution_state'").fetchone() is not None


def counts(connection):
    """Logical common-query rows; the source indexes remain stored only once."""
    if not _has_attribution(connection):
        return {}
    tables = dict(connection.execute("SELECT name,i FROM evidence_table"))
    mappings = connection.execute("SELECT count(*) FROM search_peptide").fetchone()[0]
    mappings += sum(connection.execute("SELECT count(*) FROM evidence_record WHERE table_i=?", (tables[name],)).fetchone()[0]
                    for name in _REPORTS if name in tables)
    if _REPORTS['TargetProteinGroups'] not in tables:
        return dict(mapping=mappings, target_relation=0)
    relations = connection.execute(
        "SELECT count(*) FROM (SELECT l.left_i,CASE WHEN p.status='exact_sequence' THEN p.accession ELSE r.source_key END target,"
        "CASE WHEN p.status='exact_sequence' THEN 0 ELSE 1 END namespace "
        # Keep the sparse membership table first even before ANALYZE has run.
        "FROM evidence_record x CROSS JOIN evidence_link l ON l.record_i=x.i "
        "CROSS JOIN evidence_record r ON r.i=l.right_i LEFT JOIN native_protein_search p ON p.record_i=r.i "
        "WHERE x.table_i=? AND r.table_i=? AND l.left_i IS NOT NULL GROUP BY l.left_i,target,namespace)",
        (tables[_REPORTS['TargetProteinGroups']], tables['TargetProteins'])).fetchone()[0]
    return dict(mapping=mappings, target_relation=relations)


class _Reader:
    def __init__(self, dataset, connection):
        self.dataset, self.db = dataset, connection
        self.native = CanineEvidence(dataset.path)
        self.reference = self.native.search_reference.reference
        self.protein_namespace = self.native.search_reference.namespace
        self.group_namespace = dataset.manifest.dataset_id + ":TargetProteinGroups"
        self.native_protein_namespace = dataset.manifest.dataset_id + ":TargetProteins"
        self.tables = self.native.tables

    @lru_cache(maxsize=4096)
    def protein(self, identifier):
        row = self.db.execute(
            "SELECT r.source_key,p.accession,p.status FROM evidence_record r "
            "LEFT JOIN native_protein_search p ON p.record_i=r.i WHERE r.i=? AND r.table_i=?",
            (identifier, self.tables["TargetProteins"][0])).fetchone()
        if row is None:
            return None
        key, accession, status = row
        if status == "exact_sequence":
            return ScopedTarget(self.reference, "protein", accession, self.protein_namespace)
        return ScopedTarget(self.reference, "protein", json.loads(key)[0], self.native_protein_namespace)

    def protein_ids(self, target):
        """SQL identity selection for a common destination target, without remapping it."""
        if target.kind != "protein":
            return "SELECT NULL WHERE 0", ()
        if target.namespace == self.protein_namespace:
            return ("SELECT record_i FROM native_protein_search WHERE accession=? AND status='exact_sequence'",
                    (target.identifier,))
        if target.namespace == self.native_protein_namespace:
            key = _native_key(target.identifier)
            if key is None:
                return "SELECT NULL WHERE 0", ()
            return ("SELECT i FROM evidence_record WHERE table_i=? AND source_key=?",
                    (self.tables["TargetProteins"][0], key))
        return "SELECT NULL WHERE 0", ()

    def assignments(self, table, *, source_identifier, target):
        namespace = self.group_namespace if table == "TargetProteinGroups" else "amino-acid-sequence"
        link_table = self.tables.get(_REPORTS[table], (None,))[0]
        conditions, params = ["r.table_i=?"], [self.tables[table][0]]
        if source_identifier is not None:
            if table == "TargetProteinGroups":
                key = _native_key(source_identifier)
                if key is None:
                    return
                conditions.append("r.source_key=?")
                params.append(key)
            else:
                conditions.append("r.sequence=?")
                params.append(source_identifier)
        if target is not None:
            query, bindings = self.protein_ids(target)
            conditions.append("r.i IN (SELECT l.left_i FROM evidence_link l WHERE l.right_i IN (" + query + "))")
            params.extend(bindings)
        query = ("SELECT r.i,r.source_rowid,r.source_key,r.sequence FROM evidence_record r WHERE "
                 + " AND ".join(conditions) + " ORDER BY r.source_rowid")
        for identifier, rowid, key, sequence in self.db.execute(query, params):
            links = tuple(self.db.execute(
                "SELECT l.record_i,l.right_i,r.source_rowid FROM evidence_link l "
                "JOIN evidence_record r ON r.i=l.record_i WHERE l.left_i=? AND r.table_i=? ORDER BY r.source_rowid",
                (identifier, link_table)))
            targets, missing = {}, 0
            for _, native_id, _ in links:
                protein = self.protein(native_id)
                if protein is None:
                    missing += 1
                else:
                    targets[protein.key] = protein
            found = tuple(sorted(targets.values(), key=lambda value: value.key))
            incomplete = missing or link_table is None
            status = "unresolved" if incomplete else "unique" if len(found) == 1 else "ambiguous" if found else "unmapped"
            reason = ("Native membership table is absent from the selected source" if link_table is None else
                      "One or more reported protein endpoints are absent from the selected source" if missing else
                      "No reported protein membership" if not found else None)
            source_id = json.loads(key)[0] if table == "TargetProteinGroups" else sequence
            yield TargetMapping(
                "reported:" + table + ":" + ":".join(json.loads(key)), namespace, source_id,
                self.reference, found, status, "source-reported protein membership", reason,
                evidence=dict(basis="reported", partition="target", unresolved_links=missing if link_table is not None else None,
                              membership_records=len(links) if link_table is not None else None,
                              unique_protein_detection=False),
                provenance=dict(source_asset="evidence.sqlite", native_table=table, native_key=json.loads(key),
                                native_rowid=rowid, link_table=_REPORTS[table],
                                link_rowids=[row[2] for row in links]))

    def candidates(self, *, source_identifier, target):
        conditions, params = [], []
        if source_identifier is not None:
            conditions.append("p.sequence=?")
            params.append(source_identifier)
        if target is not None:
            if target.kind != "protein" or target.namespace != self.protein_namespace:
                return
            conditions.append("p.sequence IN (SELECT sequence FROM search_peptide_hit WHERE accession=?)")
            params.append(target.identifier)
        query = "SELECT p.sequence,p.status FROM search_peptide p"
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        for sequence, state in self.db.execute(query + " ORDER BY p.sequence", params):
            placements = [dict(accession=a, start=start, end=end, basis=basis) for a, start, end, basis in self.db.execute(
                "SELECT accession,start,end,basis FROM search_peptide_hit WHERE sequence=? ORDER BY accession,start,end",
                (sequence,))]
            targets = tuple(ScopedTarget(self.reference, "protein", name, self.protein_namespace)
                            for name in sorted({p["accession"] for p in placements}))
            if state not in ("matched", "no_search_match") or bool(targets) != (state == "matched"):
                raise ValueError("Incomplete or inconsistent searched-reference candidates")
            yield TargetMapping(
                "sequence:" + sequence, "amino-acid-sequence", sequence, self.reference, targets,
                "unique" if len(targets) == 1 else "ambiguous" if targets else "unmapped",
                "exact and isoleucine/leucine-compatible searched-reference placements",
                None if targets else "No compatible placement in the complete deposited searched FASTA",
                evidence=dict(basis="sequence_compatibility", placements=placements,
                              coordinates="zero-based half-open", searched_reference_complete=True,
                              unique_protein_detection=False),
                provenance=dict(source_asset="search.fasta", search_sha256=self.native.search_reference.source.sha256))


def mappings(dataset, *, scope, source_identifier, source_namespace, target):
    with closing(dataset._connect()) as connection:
        if not _has_attribution(connection):
            return
        reader = _Reader(dataset, connection)
        if scope is not None and scope != reader.reference.key:
            return
        try:
            for table in _REPORTS:
                if table not in reader.tables:
                    continue
                namespace = reader.group_namespace if table == "TargetProteinGroups" else "amino-acid-sequence"
                if source_namespace is None or source_namespace == namespace:
                    yield from reader.assignments(table, source_identifier=source_identifier, target=target)
            if source_namespace in (None, "amino-acid-sequence"):
                yield from reader.candidates(source_identifier=source_identifier, target=target)
        finally:
            reader.protein.cache_clear()


def relations(dataset, *, scope, ids, left_target, right_target):
    """Reported native group memberships as ordinary, non-quantifying relations."""
    with closing(dataset._connect()) as connection:
        if not _has_attribution(connection) or ids == ():
            return
        reader = _Reader(dataset, connection)
        if scope is not None and scope != reader.reference.key:
            return
        if left_target is not None and (left_target.kind != "protein_group"
                                        or left_target.namespace != reader.group_namespace):
            return
        selected = None if ids is None else set(ids)
        try:
            for mapping in reader.assignments("TargetProteinGroups", source_identifier=(
                    None if left_target is None else left_target.identifier), target=right_target):
                group = ScopedTarget(reader.reference, "protein_group", mapping.source_identifier, reader.group_namespace)
                for protein in mapping.targets:
                    if right_target is not None and protein.key != right_target.key:
                        continue
                    identity = hashlib.sha256(json.dumps(protein.key, separators=(",", ":")).encode()).hexdigest()
                    identifier = "member:" + mapping.source_identifier + ":" + identity
                    if selected is not None and identifier not in selected:
                        continue
                    yield TargetRelation(
                        identifier, dataset.manifest.dataset_id, dataset.manifest.version,
                        group, protein, "shared_group_member", "many_to_many",
                        dict(mapping.evidence, source_provenance=dict(mapping.provenance),
                             membership_status=mapping.status,
                             protein_summary=reader.native.metadata.get("protein_summary")))
        finally:
            reader.protein.cache_clear()
