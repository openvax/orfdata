"""Native membership links and complete searched-reference peptide candidates.

Reported source links and sequence-compatible candidates are independent facts.
Neither converts a shared peptide or protein group into a uniquely measured
protein. Decoy records retain their partition and are not biological assignments.
Only the explicit builder uses the optional sequence-matching dependency.
"""

import hashlib
import json
from contextlib import closing
from functools import lru_cache

from ..dataset import _json
from ..measurements import DatasetManifest
from ..reference import ScopedTarget
from ..search_reference import SearchReference, SearchReferenceSpec, SearchSequence
from .pfizer_evidence import EvidenceStore, record_keys, source_key
from .pfizer_vendor import peptide_candidates

ENCODING = "orfdata.pfizer-attribution.v1"
# These endpoint columns were verified in the deposited SQLite schema and are
# retained in dog_result_plan.json. This is not table-name-based schema discovery.
ENDPOINTS = {
    "TargetPeptideGroupsTargetProteins": ("TargetPeptideGroups", "TargetProteins"),
    "TargetProteinGroupsTargetProteins": ("TargetProteinGroups", "TargetProteins"),
    "TargetProteinGroupsTargetPeptideGroups": ("TargetProteinGroups", "TargetPeptideGroups"),
    "TargetProteinGroupsTargetPsms": ("TargetProteinGroups", "TargetPsms"),
    "TargetPsmsTargetPeptideGroups": ("TargetPsms", "TargetPeptideGroups"),
    "TargetPsmsTargetProteins": ("TargetPsms", "TargetProteins"),
    "FoundModificationsTargetPsms": ("FoundModifications", "TargetPsms"),
    "TargetPsmsMSnSpectrumInfo": ("TargetPsms", "MSnSpectrumInfo"),
    "TargetPsmsQuanSpectrumInfo": ("TargetPsms", "QuanSpectrumInfo"),
    "QuanSpectrumInfoTargetPeptideGroups": ("QuanSpectrumInfo", "TargetPeptideGroups"),
    "DecoyPeptideGroupsDecoyProteins": ("DecoyPeptideGroups", "DecoyProteins"),
    "DecoyProteinGroupsDecoyProteins": ("DecoyProteinGroups", "DecoyProteins"),
    "DecoyProteinGroupsDecoyPeptideGroups": ("DecoyProteinGroups", "DecoyPeptideGroups"),
    "DecoyProteinGroupsDecoyPsms": ("DecoyProteinGroups", "DecoyPsms"),
    "DecoyPsmsDecoyPeptideGroups": ("DecoyPsms", "DecoyPeptideGroups"),
    "DecoyPsmsDecoyProteins": ("DecoyPsms", "DecoyProteins"),
    "FoundModificationsDecoyPsms": ("FoundModifications", "DecoyPsms"),
    "DecoyPsmsMSnSpectrumInfo": ("DecoyPsms", "MSnSpectrumInfo"),
    "QuanSpectrumInfoReporterQuanResults": ("QuanSpectrumInfo", "ReporterQuanResults"),
}
SCHEMA = """
CREATE TABLE evidence_link(record_i INTEGER PRIMARY KEY REFERENCES evidence_record(i),
    left_i INTEGER REFERENCES evidence_record(i), right_i INTEGER REFERENCES evidence_record(i));
CREATE INDEX evidence_link_left ON evidence_link(left_i,record_i);
CREATE INDEX evidence_link_right ON evidence_link(right_i,record_i);
CREATE TABLE search_protein(accession TEXT PRIMARY KEY, record TEXT NOT NULL) WITHOUT ROWID;
CREATE TABLE native_protein_search(record_i INTEGER PRIMARY KEY REFERENCES evidence_record(i),
    accession TEXT REFERENCES search_protein(accession), status TEXT NOT NULL);
CREATE INDEX native_protein_search_accession ON native_protein_search(accession,record_i);
CREATE TABLE search_peptide(sequence TEXT PRIMARY KEY,
    status TEXT NOT NULL CHECK(status IN ('pending','matched','no_search_match'))) WITHOUT ROWID;
CREATE TABLE search_peptide_hit(sequence TEXT NOT NULL REFERENCES search_peptide(sequence),
    accession TEXT NOT NULL REFERENCES search_protein(accession), start INTEGER NOT NULL,
    end INTEGER NOT NULL, basis TEXT NOT NULL,
    PRIMARY KEY(sequence,accession,start,end), CHECK(start>=0 AND end>start),
    CHECK(basis IN ('exact_residues','isoleucine_leucine_equivalent'))) WITHOUT ROWID;
CREATE INDEX search_peptide_hit_protein ON search_peptide_hit(accession,sequence,start);
"""


def build_attribution(connection, searched, *, max_residues=250000, checkpoint=None):
    """Index complete local evidence and every placement in its verified FASTA.

    The caller owns the new output, its SQLite/disk limits and failure cleanup.
    Source records remain unchanged. Dangling native links and unmatched target
    proteins are reported explicitly; they must be reviewed in release acceptance.
    """
    if not isinstance(searched, SearchReference):
        raise TypeError("Attribution requires an already checksum-verified SearchReference")
    if searched.reference.taxon_id != 9615:
        raise ValueError("Native canine evidence requires a canine searched reference")
    if type(max_residues) is not int or max_residues < 1:
        raise ValueError("Peptide batch residue limit must be positive")
    state = connection.execute("SELECT value FROM evidence_metadata WHERE name='state'").fetchone()
    if state is None or json.loads(state[0]) != "complete":
        raise ValueError("Finish evidence normalization before building attributions")
    if connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='metadata'").fetchone():
        stored = connection.execute("SELECT value FROM metadata WHERE name='manifest'").fetchone()
        manifest = DatasetManifest.from_dict(json.loads(stored[0])) if stored else None
        if manifest is None or searched.reference not in manifest.references:
            raise ValueError("Searched FASTA reference differs from the quantitative dataset")
    tables = {name: index for index, name in connection.execute("SELECT i,name FROM evidence_table")}
    connection.executescript(SCHEMA)
    connection.execute("INSERT INTO evidence_metadata VALUES ('attribution_state','\"building\"')")
    connection.execute("INSERT INTO evidence_metadata VALUES ('attribution_encoding',?)", (_json(ENCODING),))
    connection.execute("INSERT INTO evidence_metadata VALUES ('search_reference',?)", (_json(searched.spec.to_dict()),))
    counts = dict(native_links=0, unresolved_links=0, searched_proteins=0, target_proteins=0,
                  exact_native_proteins=0, unmatched_native_proteins=0, decoy_proteins=0,
                  peptide_sequences=0, no_search_match=0, exact_placements=0, il_equivalent_placements=0)
    for protein in searched.entries():
        record = protein.to_dict()
        del record["target"]
        connection.execute("INSERT INTO search_protein VALUES (?,?)", (protein.target.identifier, _json(record)))
        counts["searched_proteins"] += 1
    for table in ("TargetProteins", "DecoyProteins"):
        if table not in tables:
            continue
        for index, attributes in connection.execute("SELECT i,attributes FROM evidence_record WHERE table_i=? ORDER BY source_rowid",
                                                    (tables[table],)):
            native = json.loads(attributes)
            if table == "DecoyProteins":
                accession, status = None, "decoy_partition"
                counts["decoy_proteins"] += 1
            else:
                protein = searched.get(native["Accession"])
                status = ("not_in_search" if protein is None else "exact_sequence" if protein.sequence == native["Sequence"]
                          else "sequence_mismatch")
                accession = protein.target.identifier if protein else None
                counts["target_proteins"] += 1
                counts["exact_native_proteins" if status == "exact_sequence" else "unmatched_native_proteins"] += 1
            connection.execute("INSERT INTO native_protein_search VALUES (?,?,?)", (index, accession, status))
    connection.commit()

    @lru_cache(maxsize=8192)
    def endpoint(table, key):
        if table not in tables:
            raise ValueError("Native link endpoint table is absent from the selected evidence: " + table)
        rows = connection.execute("SELECT i FROM evidence_record WHERE table_i=? AND source_key=?",
                                  (tables[table], source_key(key))).fetchmany(2)
        if len(rows) > 1:
            raise ValueError("Native endpoint identity is not unique: " + table)
        return rows[0][0] if rows else None

    for table, endpoints in ENDPOINTS.items():
        if table not in tables:
            continue
        for identifier, attributes in connection.execute(
                "SELECT i,attributes FROM evidence_record WHERE table_i=? ORDER BY source_rowid", (tables[table],)):
            native = json.loads(attributes)
            linked = [endpoint(name, tuple(native[name + key] for key in record_keys(name))) for name in endpoints]
            connection.execute("INSERT INTO evidence_link VALUES (?,?,?)", (identifier, *linked))
            counts["native_links"] += 1
            counts["unresolved_links"] += int(None in linked)
            if counts["native_links"] % 5000 == 0:
                connection.commit()
                if checkpoint is not None:
                    checkpoint(dict(counts))
        connection.commit()
    endpoint.cache_clear()
    # A disk-backed primary key deduplicates across modified peptide groups and
    # PSMs without accumulating the full peptide universe in Python memory.
    for table in ("TargetPeptideGroups", "TargetPsms"):
        if table in tables:
            connection.execute("INSERT OR IGNORE INTO search_peptide SELECT sequence,'pending' FROM evidence_record "
                               "WHERE table_i=? AND sequence IS NOT NULL", (tables[table],))
    connection.commit()

    def batch(peptides):
        matched = set()
        for hit in peptide_candidates(peptides, searched, max_residues=max_residues):
            if hit["status"] == "candidate":
                connection.execute("INSERT INTO search_peptide_hit VALUES (?,?,?,?,?)",
                                   tuple(hit[key] for key in ("peptide", "accession", "start", "end", "basis")))
                matched.add(hit["peptide"])
                counts["exact_placements" if hit["basis"] == "exact_residues" else "il_equivalent_placements"] += 1
        connection.executemany("UPDATE search_peptide SET status=? WHERE sequence=?",
                               (("matched" if peptide in matched else "no_search_match", peptide) for peptide in peptides))
        counts["peptide_sequences"] += len(peptides)
        counts["no_search_match"] += len(peptides) - len(matched)
        connection.commit()
        if checkpoint is not None:
            checkpoint(dict(counts))

    peptides, residues = [], 0
    for (sequence,) in connection.execute("SELECT sequence FROM search_peptide ORDER BY sequence"):
        if len(sequence) > max_residues:
            raise ValueError("Single native peptide exceeds the explicit batch residue limit")
        if peptides and residues + len(sequence) > max_residues:
            batch(peptides)
            peptides, residues = [], 0
        peptides.append(sequence)
        residues += len(sequence)
    if peptides:
        batch(peptides)
    if connection.execute("SELECT 1 FROM search_peptide WHERE status='pending' LIMIT 1").fetchone():
        raise ValueError("Not all native target peptide sequences were searched")
    if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
        raise ValueError("Native attribution indexes have invalid references")
    candidate_digest = hashlib.sha256()
    for row in connection.execute(
            "SELECT p.sequence,p.status,h.accession,h.start,h.end,h.basis FROM search_peptide p "
            "LEFT JOIN search_peptide_hit h ON h.sequence=p.sequence ORDER BY p.sequence,h.accession,h.start,h.end"):
        encoded = _json(row).encode()
        candidate_digest.update(len(encoded).to_bytes(8, "little") + encoded)
    connection.execute("INSERT INTO evidence_metadata VALUES ('attribution_counts',?)", (_json(counts),))
    connection.execute("INSERT INTO evidence_metadata VALUES ('candidate_digest',?)", (_json(candidate_digest.hexdigest()),))
    connection.execute("UPDATE evidence_metadata SET value='\"complete\"' WHERE name='attribution_state'")
    connection.execute("ANALYZE")
    connection.commit()
    return counts


class CanineEvidence(EvidenceStore):
    """Local native records, membership links and complete FASTA candidates.

    Exclusion and confidence values remain source facts; no implicit detection
    threshold is applied. Callers must use those facts for their own evidence
    policy. Bulk measurements provide no DLA-presentation observation.
    """

    def __init__(self, path):
        super().__init__(path)
        if (self.metadata.get("attribution_encoding") != ENCODING
                or self.metadata.get("attribution_state") != "complete"):
            raise ValueError("Canine membership/search evidence is incomplete or unsupported")
        self.search_reference = SearchReferenceSpec.from_dict(self.metadata["search_reference"])

    @staticmethod
    def _record(db, identifier):
        if identifier is None:
            return None
        row = db.execute("SELECT r.i,r.source_rowid,r.attributes,t.name,t.role FROM evidence_record r "
                         "JOIN evidence_table t ON t.i=r.table_i WHERE r.i=?", (identifier,)).fetchone()
        return dict(table=row[3], role=row[4], record=EvidenceStore.restore(db, row[:3]))

    def linked(self, table, key, *, via=None, direction="outgoing"):
        """Return each original link and both endpoints, including absent endpoints.

        Incoming links are useful for finding protein groups or peptide records
        mentioning a protein. `via` names a native link table, never an inferred
        relationship. The source link fields retain all native positions/usage.
        """
        if direction not in ("outgoing", "incoming"):
            raise ValueError("Native link direction must be outgoing or incoming")
        if via is not None and via not in ENDPOINTS:
            raise ValueError("Unknown native relationship table")
        column = "left_i" if direction == "outgoing" else "right_i"
        with closing(self.connect()) as db:
            rows = db.execute("SELECT i FROM evidence_record WHERE table_i=? AND source_key=?",
                              (self.tables[table][0], source_key(tuple(key)))).fetchmany(2)
            if not rows:
                return
            if len(rows) != 1:
                raise ValueError("Native record identity is not unique")
            query = ("SELECT l.record_i,l.left_i,l.right_i FROM evidence_link l "
                     "JOIN evidence_record r ON r.i=l.record_i WHERE l." + column + "=?")
            params = [rows[0][0]]
            if via is not None:
                query += " AND r.table_i=?"
                if via not in self.tables:
                    return
                params.append(self.tables[via][0])
            query += " ORDER BY r.table_i,r.source_rowid"
            for link, left, right in db.execute(query, params):
                yield dict(source_link=self._record(db, link), left=self._record(db, left), right=self._record(db, right),
                           status="resolved" if left is not None and right is not None else "unresolved")

    def candidates(self, sequence):
        """All exact/I-L-compatible placements, or explicit no-match/unsearched state."""
        with closing(self.connect()) as db:
            row = db.execute("SELECT status FROM search_peptide WHERE sequence=?", (sequence,)).fetchone()
            if row is None or row[0] == "no_search_match":
                yield dict(peptide=sequence, accession=None, start=None, end=None,
                           status="not_searched" if row is None else "no_search_match", basis="exact_or_il_equivalent")
                return
            for accession, start, end, basis in db.execute(
                    "SELECT accession,start,end,basis FROM search_peptide_hit WHERE sequence=? ORDER BY accession,start,end", (sequence,)):
                yield dict(peptide=sequence, accession=accession, start=start, end=end, status="candidate", basis=basis)

    def search_protein(self, accession):
        """An original searched sequence; presence in the FASTA is not detection."""
        with closing(self.connect()) as db:
            row = db.execute("SELECT record FROM search_protein WHERE accession=?", (accession,)).fetchone()
            if row is None:
                return None
            values = json.loads(row[0])
            values["target"] = ScopedTarget(self.search_reference.reference, "protein", accession, self.search_reference.namespace)
            return SearchSequence(**values)

    def protein_search_match(self, table, key):
        """Validate native target protein identity, retaining decoys and mismatches."""
        if table not in ("TargetProteins", "DecoyProteins"):
            raise ValueError("Search identity validation requires a native protein record")
        with closing(self.connect()) as db:
            rows = db.execute("SELECT p.accession,p.status FROM native_protein_search p "
                              "JOIN evidence_record r ON r.i=p.record_i WHERE r.table_i=? AND r.source_key=?",
                              (self.tables[table][0], source_key(tuple(key)))).fetchmany(2)
            if len(rows) != 1:
                raise ValueError("Expected one native protein search reconciliation")
            return dict(accession=rows[0][0], status=rows[0][1],
                        search_fasta_sha256=self.search_reference.source.sha256)
