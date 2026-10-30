"""Exact source slices and labelled mutations test the pinned reference join."""

import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import Proteome, ReferenceContext

FIXTURE = Path(__file__).parent / "data" / "ensembl"
RNA = "ENST00000303766.12"
PROTEIN = "ENSP00000307155.7"
GENE = "ENSG00000169800"


def inputs(tmp_path, *, plain=False):
    paths = {}
    for role in ("protein", "gtf", "rna", "ncrna"):
        raw = (FIXTURE / f"{role}.gz").read_bytes()
        paths[role] = tmp_path / role
        paths[role].write_bytes(gzip.decompress(raw) if plain else raw)
    return paths


def context(**kwargs):
    return replace(ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115",
        assembly_accession="GCA_000001405.29"), **kwargs)


def build(paths, *, reference=None, hashes=None):
    hashes = hashes or {role: hashlib.sha256(path.read_bytes()).hexdigest() for role, path in paths.items()}
    return Proteome.from_ensembl(paths["protein"], paths["gtf"], paths["rna"], paths["ncrna"],
        reference=reference or context(), expected_sha256=hashes)


def mutate(paths, role, change):
    raw = paths[role].read_bytes()
    text = (gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw).decode()
    paths[role].write_text(change(text))


def fasta_records(text):
    return {part.split()[0]: ">" + part for part in text.split(">")[1:]}


@pytest.mark.parametrize("plain", [False, True])
@pytest.mark.parametrize("annotation_name", ["Ensembl", "ensembl"])
def test_pinned_real_ensembl_pairs_and_all_source_identities_replay_offline(tmp_path, plain, annotation_name):
    paths = inputs(tmp_path, plain=plain)
    if not plain:
        pins = json.loads((FIXTURE / "manifest.json").read_text())["fixture_sha256"]
        assert pins == {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in paths.items()}
    result = build(paths, reference=context(annotation_name=annotation_name))
    assert len(result) == 13
    assert result.reference.annotation_name == annotation_name
    first = result.by_accession()[PROTEIN]
    assert (first.gene_id, first.transcript_id, len(first.sequence)) == (GENE, RNA, 496)
    duplicate_sequence = result.by_accession()["ENSP00000250831.5"]
    assert first.sequence == duplicate_sequence.sequence
    assert first.gene_id != duplicate_sequence.gene_id
    assert {o.transcript_id for e in result for o in e.occurrences} >= {RNA, "ENST00000454978.2"}
    assert all((o.start, o.end, o.coding_blocks) == (None, None, ()) for e in result for o in e.occurrences)
    layer = result.layers[0]
    assert layer["annotated_transcripts"] == layer["rna_records"] == 19
    assert set(layer["content_provenance"]) == {"protein", "gtf", "rna", "ncrna"}
    assert len(result.reference.source_sha256) == len(result.reference.sequence_sha256) == 64
    saved = tmp_path / "proteome.json.gz"
    result.save(saved)
    for path in paths.values():
        path.unlink()
    replay = Proteome.load(saved)
    assert (replay.identity, replay.proteins, replay.layers) == (result.identity, result.proteins, result.layers)


@pytest.mark.parametrize("role", ["protein", "gtf", "rna", "ncrna"])
def test_bad_source_hash_fails_before_an_inventory_is_returned(tmp_path, role):
    paths = inputs(tmp_path)
    pins = {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in paths.items()}
    pins[role] = "0" * 64
    with pytest.raises(ValueError, match="checksum mismatch"):
        build(paths, hashes=pins)


@pytest.mark.parametrize("change", [dict(assembly_accession="GCA_000001405.28"),
    dict(annotation_name="NCBI RefSeq"), dict(source_sha256="0" * 64), dict(sequence_sha256="0" * 64)])
def test_reference_metadata_must_match_pinned_content(tmp_path, change):
    with pytest.raises(ValueError, match="assembly|Ensembl annotation|checksums"):
        build(inputs(tmp_path), reference=context(**change))


@pytest.mark.parametrize("role,old,new", [
    ("protein", PROTEIN, "ENSP00000307155.8"), ("rna", RNA, "ENST00000303766.13"),
    ("rna", "gene:ENSG00000169800.14", "gene:ENSG00000226941.9"),
    ("protein", "transcript:" + RNA, "transcript:ENST00000250831.6"),
    ("gtf", 'transcript_version "12"', 'transcript_version "13"'),
    ("gtf", 'protein_version "7"', 'protein_version "8"')])
def test_exact_versions_and_gene_assignments_are_not_heuristic_joins(tmp_path, role, old, new):
    paths = inputs(tmp_path)
    def change(text):
        assert old in text
        return text.replace(old, new)
    mutate(paths, role, change)
    with pytest.raises(ValueError, match="inventories disagree|identity conflict|exact RNA identity"):
        build(paths)


@pytest.mark.parametrize("role,identifier", [("protein", PROTEIN), ("rna", RNA)])
def test_missing_coding_sequences_fail_and_conflicting_duplicates_cannot_overwrite(tmp_path, role, identifier):
    paths = inputs(tmp_path)
    def missing(text):
        records = fasta_records(text)
        records.pop(identifier)
        return "".join(records.values())
    mutate(paths, role, missing)
    with pytest.raises(ValueError, match="inventories disagree|exact RNA identity"):
        build(paths)
    paths = inputs(tmp_path)
    def duplicate(text):
        record = fasta_records(text)[identifier]
        return text + record.rstrip() + "A\n"
    mutate(paths, role, duplicate)
    with pytest.raises(ValueError, match="Conflicting Ensembl FASTA identity"):
        build(paths)


def test_identical_fasta_duplicates_are_recorded_without_merging_distinct_ids(tmp_path):
    paths = inputs(tmp_path)
    mutate(paths, "protein", lambda text: text + fasta_records(text)[PROTEIN])
    result = build(paths)
    assert len(result) == 13
    assert result.layers[0]["identical_protein_fasta_duplicates"] == 1


def test_absent_untranslated_rna_is_explicit_and_cannot_become_a_protein(tmp_path):
    paths = inputs(tmp_path)
    original = build(paths)
    coding = {o.transcript_id for e in original for o in e.occurrences}
    removed = []
    def change(text):
        records = fasta_records(text)
        removed.append(next(key for key in records if key not in coding))
        records.pop(removed[0])
        return "".join(records.values())
    mutate(paths, "rna", change)
    result = build(paths)
    assert result.proteins == original.proteins
    assert result.layers[0]["untranslated_transcripts_without_rna"] == removed


def test_repeated_gtf_loci_survive_as_distinct_occurrences_and_spans(tmp_path):
    paths = inputs(tmp_path)
    def change(text):
        extra = []
        for line in text.splitlines():
            if 'transcript_id "ENST00000303766";' not in line:
                continue
            fields = line.split("\t")
            fields[3], fields[4] = str(int(fields[3]) + 1000000), str(int(fields[4]) + 1000000)
            extra.append("\t".join(fields))
        return text + "\n".join(extra) + "\n"
    mutate(paths, "gtf", change)
    result = build(paths)
    entry = result.by_accession()[PROTEIN]
    assert len(entry.occurrences) == 2
    spans = result.layers[0]["transcript_spans"][PROTEIN]
    assert spans[1] == [spans[0][0] + 1000000, spans[0][1] + 1000000]
    saved = tmp_path / "repeated.json.gz"
    result.save(saved)
    assert Proteome.load(saved).by_accession()[PROTEIN].occurrences == entry.occurrences


def test_cds_outside_its_transcript_and_conflicting_metadata_are_rejected(tmp_path):
    paths = inputs(tmp_path)
    def outside(text):
        lines = text.splitlines()
        for index, line in enumerate(lines):
            if '\tCDS\t' in line and 'transcript_id "ENST00000303766";' in line:
                fields = line.split("\t")
                fields[3], fields[4] = "1", "3"
                lines[index] = "\t".join(fields)
                break
        return "\n".join(lines) + "\n"
    mutate(paths, "gtf", outside)
    with pytest.raises(ValueError, match="containing transcript locus"):
        build(paths)
    paths = inputs(tmp_path)
    def conflict(text):
        row = next(line for line in text.splitlines() if '\ttranscript\t' in line and 'transcript_id "ENST00000303766";' in line)
        return text + row.replace('gene_name "RBMY1F"', 'gene_name "conflicting fixture"') + "\n"
    mutate(paths, "gtf", conflict)
    with pytest.raises(ValueError, match="Conflicting Ensembl transcript metadata"):
        build(paths)


def test_parsing_reads_the_open_pinned_file_even_if_its_path_is_replaced(tmp_path):
    from orfdata.ensembl import _source
    path = tmp_path / "source"
    content = b"first\nsecond\n"
    path.write_bytes(gzip.compress(content, mtime=0))
    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    with _source(path, expected) as stream:
        assert stream.readline() == "first\n"
        replacement = tmp_path / "replacement"
        replacement.write_text("different bytes\n")
        replacement.replace(path)
        assert stream.read() == "second\n"


def test_real_ensembl_internal_stops_remain_reference_membership():
    root = FIXTURE / "internal-stop"
    pins = json.loads((root / "manifest.json").read_text())["fixture_sha256"]
    result = build({role: root / f"{role}.gz" for role in pins}, hashes=pins)
    entry, = result.proteins
    assert entry.accession == "ENSP00000488168.1"
    raw = gzip.decompress((root / "protein.gz").read_bytes()).decode()
    assert entry.sequence == "".join(raw.splitlines()[1:])
    assert len(entry.sequence) == 6 and [i for i, value in enumerate(entry.sequence) if value == "*"] == [0, 4]
    assert entry.evidence == ()


def test_source_headers_allow_utf8_while_sequence_identity_stays_exact(tmp_path):
    paths = inputs(tmp_path)
    mutate(paths, "protein", lambda text: text.replace("description:", "description:β "))
    assert len(build(paths)) == 13
