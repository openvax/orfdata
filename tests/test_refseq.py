"""RefSeq source accession joins preserve identities without a symbol heuristic."""

import csv
import gzip
import hashlib
from dataclasses import replace

import pytest

from orfdata import Proteome, ReferenceContext


def source(tmp_path, *, rows=None, fasta=None, compressed=False):
    fields = ['# feature', 'class', 'assembly', 'product_accession', 'related_accession',
              'GeneID', 'symbol', 'genomic_accession', 'strand', 'product_length', 'start', 'end']
    if rows is None:
        rows = [dict(zip(fields, ['CDS', 'with_protein', 'GCF_011100685.1', 'XP_123.1', 'XM_123.1',
                    '123', 'GENE', 'NC_123.1', '+', '6', '1', '200'])),
                dict(zip(fields, ['CDS', 'with_protein', 'GCF_011100685.1', 'XP_123.1', 'XM_124.2',
                    '124', 'OTHER', 'NC_124.1', '-', '6', '2', '201']))]
    protein, features = tmp_path / 'protein.faa', tmp_path / 'feature_table.txt'
    protein.write_text(fasta if fasta is not None else '>XP_123.1 source description\nACDEFG\n')
    with features.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    if compressed:
        for path in (protein, features):
            path.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
    hashes = {role: hashlib.sha256(path.read_bytes()).hexdigest()
              for role, path in (('protein', protein), ('feature_table', features))}
    context = ReferenceContext(9615, 'dog', 'UU_Cfam_GSD_1.0', 'NCBI RefSeq', '106',
                               assembly_accession='GCF_011100685.1')
    return protein, features, context, hashes, rows


def build(args):
    protein, features, reference, hashes, _ = args
    return Proteome.from_refseq(protein, features, reference=reference, expected_sha256=hashes)


@pytest.mark.parametrize('compressed', [False, True])
def test_refseq_keeps_all_gene_occurrences_and_offline_provenance(tmp_path, compressed):
    args = source(tmp_path, compressed=compressed)
    result = build(args)
    entry, = result.entries
    assert entry.sequence == 'ACDEFG'
    assert [o.gene_id for o in entry.occurrences] == ['123', '124']
    assert [o.transcript_id for o in entry.occurrences] == ['XM_123.1', 'XM_124.2']
    assert all((o.start, o.end, o.coding_blocks, o.transcript_biotype) == (None, None, (), None)
               for o in entry.occurrences)
    assert dict(result.layers[0]['content_provenance']) == args[3]
    assert result.layers[0]['genomic_spans']['XP_123.1'] == [[1, 200], [2, 201]]
    saved = tmp_path / 'proteome.json.gz'
    result.save(saved)
    args[0].unlink()
    args[1].unlink()
    replay = Proteome.load(saved)
    assert replay.identity == result.identity
    assert replay.entries == result.entries
    assert replay.layers == result.layers


def test_refseq_preserves_repeated_placements_and_untranslated_cds_count(tmp_path):
    args = source(tmp_path)
    rows = args[-1]
    rows.append(dict(rows[0], start='1001', end='1200'))
    rows.append(dict(rows[0], **{'class': 'without_protein', 'product_accession': '', 'related_accession': ''}))
    result = build(source(tmp_path, rows=rows))
    assert len(result.entries[0].occurrences) == 3
    assert result.layers[0]['cds_without_protein'] == 1
    assert result.layers[0]['genomic_spans']['XP_123.1'] == [[1, 200], [1001, 1200], [2, 201]]


@pytest.mark.parametrize('change', [dict(assembly='GCF_011100685.2'), dict(GeneID='GENE'),
    dict(related_accession='XM_123'), dict(product_accession='XP_123'), dict(product_length='7'),
    dict(strand='?'), dict(start='200', end='1'), dict(genomic_accession='')])
def test_refseq_rejects_invalid_and_conflicting_source_fields(tmp_path, change):
    args = source(tmp_path)
    rows = args[-1]
    rows[0].update(change)
    with pytest.raises(ValueError):
        build(source(tmp_path, rows=rows))


@pytest.mark.parametrize('fasta', ['>XP_123.1\nACDEFG\n>XP_123.1\nACDEFG\n',
    '>XP_123\nACDEFG\n', 'ACDEFG\n', '>XP_123.1\nAC*EFG\n',
    '>XP_123.1\nACDEFG\n>XP_999.1\nACDEFG\n', ''])
def test_refseq_rejects_invalid_or_unmapped_fasta_proteins(tmp_path, fasta):
    with pytest.raises(ValueError):
        build(source(tmp_path, fasta=fasta))


def test_refseq_rejects_checksum_changes_and_missing_sources(tmp_path):
    args = source(tmp_path)
    args[0].write_text('changed')
    with pytest.raises(ValueError, match='checksum mismatch'):
        build(args)
    args = source(tmp_path)
    with pytest.raises(ValueError, match='requires protein and feature_table'):
        Proteome.from_refseq(args[0], args[1], reference=args[2], expected_sha256={'protein': args[3]['protein']})
    with pytest.raises(ValueError, match='checksums disagree'):
        Proteome.from_refseq(args[0], args[1], reference=replace(args[2], source_sha256='0' * 64),
                            expected_sha256=args[3])


def test_refseq_rejects_one_transcript_assigned_to_different_gene(tmp_path):
    args = source(tmp_path)
    rows = args[-1]
    rows[1]['related_accession'] = rows[0]['related_accession']
    with pytest.raises(ValueError, match='Conflicting RefSeq transcript'):
        build(source(tmp_path, rows=rows))
