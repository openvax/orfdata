"""Producer entries compose without inventing reference protein parents."""

import gzip
import json
from dataclasses import replace

import pytest

from orfdata import DerivedProtein, ProducerInfo, Proteome, ProteomeEntry, ReferenceContext


@pytest.fixture
def reference():
    return ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")


def produced(sequence="MPEPTIDE", *, taxon_id=9606, **options):
    fields = dict(origin="nc_orf", method="rna_assembly", tool="isovar",
                  native_id="ASM42", sample_id="tumor-1")
    fields.update(options)
    info = ProducerInfo(**fields)
    return ProteomeEntry.from_producer(sequence, taxon_id=taxon_id, producer=info)


def test_reference_and_two_producers_compose_without_parents(reference, tmp_path):
    base = Proteome(reference, [ProteomeEntry("ENSP1", "MPEPTIDE", "annotated", 9606,
                                             "Ensembl:115"),
                                ProteomeEntry("ORF1", "MORF", "extra_orf", 9606,
                                             "catalogue:1")],
                    [dict(kind="annotated", role="host", source="Ensembl:115",
                          reference=reference.key, added=1),
                     dict(kind="extra_orf", role="host", source="catalogue:1",
                          reference=reference.key, added=1)])
    rna = produced(reference=reference, gene_names=("GENE1",),
                   transcript_ids=("ENST1",), start_complete=False, stop_complete=None)
    dna = ProteomeEntry.from_producer("MPEPTIDE", taxon_id=9606, producer=ProducerInfo(
        origin="reference_cds", method="dna_prediction", tool="varcode", native_id="V1",
        sample_id="tumor-1", reference=reference, gene_names=("GENE1",)))
    joint = base.merge(Proteome.from_entries((rna, dna), reference=reference))
    assert len(joint) == 4 and len(joint.by_sequence()["MPEPTIDE"]) == 3
    assert rna.occurrences == dna.occurrences == ()
    assert rna.producer.start_complete is False and rna.producer.stop_complete is None
    assert rna.gene_name == "GENE1" and dna.gene_name == "GENE1"
    assert joint.with_entries((rna, dna)).identity == joint.identity
    path = tmp_path / "joint.proteome.json.gz"
    joint.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == joint.identity
    assert loaded.by_accession()[rna.accession].producer == rna.producer


def test_source_qualified_ids_distinguish_samples_tools_and_runs(reference):
    first = produced(reference=reference)
    second = produced(reference=reference, sample_id="tumor-2")
    other_tool = ProteomeEntry.from_producer("MPEPTIDE", taxon_id=9606, producer=ProducerInfo(
        origin="nc_orf", method="rna_assembly", tool="other", native_id="ASM42",
        sample_id="tumor-1", reference=reference))
    other_run = produced(reference=reference, run_id="run-2")
    assert len({item.accession for item in (first, second, other_tool, other_run)}) == 4
    combined = Proteome.from_entries((first, second, other_tool, other_run), reference=reference)
    assert len(combined) == 4
    assert combined.with_entries(first).identity == combined.identity
    with pytest.raises(ValueError, match="Conflicting proteome entry identity"):
        combined.with_entries(produced("MCHANGED", reference=reference))


def test_unknown_reference_requires_an_explicit_policy(reference, tmp_path):
    unknown = produced()
    with pytest.raises(ValueError, match="unscoped proteome"):
        Proteome.from_entries([unknown])
    mixed = Proteome.from_entries([unknown], compatibility="mixed")
    assert mixed.reference is None
    with pytest.raises(ValueError, match="explicit compatibility"):
        Proteome(reference, ()).merge(mixed)
    joint = Proteome(reference, ()).merge(mixed, compatibility="allow_unknown")
    assert joint.reference_policy == "allow_unknown"
    path = tmp_path / "unknown.proteome.json.gz"
    joint.save(path)
    assert Proteome.load(path).identity == joint.identity


def test_gene_labels_and_incomplete_boundaries_are_not_genomic_occurrences(reference):
    fusion = produced(reference=reference, gene_names=("GENE1", "GENE2"),
                      gene_ids=("G1", "G2"), transcript_ids=("T1", "T2"),
                      orf_ids=("O1", "O2"), start_complete=None, stop_complete=False)
    missing = produced(reference=reference, native_id="NO_GENE", gene_names=())
    assert fusion.gene_name is None and fusion.occurrences == ()
    assert fusion.producer.gene_names == ("GENE1", "GENE2")
    assert fusion.producer.stop_complete is False
    assert missing.gene_name is None
    snapshot = Proteome.from_entries((fusion, missing), reference=reference)
    assert len(snapshot) == 2


def test_wrong_taxon_or_reference_and_conflicting_merge_fail(reference):
    another = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "114")
    with pytest.raises(ValueError, match="another reference"):
        Proteome.from_entries([produced(reference=another)], reference=reference)
    with pytest.raises(ValueError, match="another species"):
        Proteome.from_entries([produced(reference=None, taxon_id=10090)],
                              reference=reference, compatibility="allow_unknown")
    first = Proteome.from_entries([produced(reference=reference)], reference=reference)
    changed = replace(first.entries[0], sequence="MCHANGED")
    with pytest.raises(ValueError, match="Conflicting proteome entry identity"):
        first.merge(Proteome.from_entries([changed], reference=reference))
    with pytest.raises(ValueError, match="different reference"):
        first.merge(Proteome.from_entries([produced(reference=another)], reference=another))


def test_legacy_v2_snapshot_still_loads(reference, tmp_path):
    ordinary = Proteome(reference, [ProteomeEntry("ENSP1", "MAA", "annotated", 9606,
                                                 "Ensembl:115")])
    path = tmp_path / "legacy.json.gz"
    ordinary.save(path)
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    snapshot["schema"] = "orfdata.proteome.v2"
    path.write_bytes(gzip.compress(json.dumps(snapshot).encode(), mtime=0))
    assert Proteome.load(path).identity == ordinary.identity


def test_external_viral_origin_keeps_the_host_boundary(reference):
    source_reference = ReferenceContext(999999, "synthetic_virus", "toy-v1", "source", "1")
    entry = ProteomeEntry.from_producer("MAA", taxon_id=999999, producer=ProducerInfo(
        origin="viral", method="rna_assembly", tool="fixture", native_id="A",
        reference=source_reference, host_taxon_ids=(9606,), host_source="synthetic fixture"))
    host = Proteome.from_entries([entry], reference=reference)
    assert len(host.host()) == 0
    with pytest.raises(ValueError, match="viral parents"):
        host.with_derived_proteins(DerivedProtein(
            "changed", "MAB", (entry.accession,), "sample-1", ("change-1",), "fixture"))
    mouse = ReferenceContext(10090, "mouse", "GRCm39", "Ensembl", "115")
    with pytest.raises(ValueError, match="recorded host"):
        Proteome.from_entries([entry], reference=mouse)
