"""Producer entries compose without inventing reference protein parents."""

import gzip
import json
from dataclasses import replace

import pytest

from orfdata import (
    DerivedProtein,
    ProducerInfo,
    ProducerLink,
    ProducerRecord,
    ProteinSequence,
    Proteome,
    ReferenceContext,
)
from orfdata.proteome import _sha256


@pytest.fixture
def reference():
    return ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")


def produced(sequence="MPEPTIDE", *, taxon_id=9606, **options):
    fields = dict(origin="nc_orf", method="rna_assembly", tool="isovar",
                  native_id="ASM42", sample_id="tumor-1")
    fields.update(options)
    info = ProducerInfo(**fields)
    return ProteinSequence.from_producer(sequence, taxon_id=taxon_id, producer=info)


def test_reference_and_two_producers_compose_without_parents(reference, tmp_path):
    base = Proteome(reference, [ProteinSequence("ENSP1", "MPEPTIDE", "annotated", 9606,
                                             "Ensembl:115"),
                                ProteinSequence("ORF1", "MORF", "extra_orf", 9606,
                                             "catalogue:1")],
                    [dict(kind="annotated", role="host", source="Ensembl:115",
                          reference=reference.key, added=1),
                     dict(kind="extra_orf", role="host", source="catalogue:1",
                          reference=reference.key, added=1)])
    rna = produced(reference=reference, gene_names=("GENE1",),
                   transcript_ids=("ENST1",), start_complete=False, stop_complete=None)
    dna = ProteinSequence.from_producer("MPEPTIDE", taxon_id=9606, producer=ProducerInfo(
        origin="reference_cds", method="dna_prediction", tool="varcode", native_id="V1",
        sample_id="tumor-1", reference=reference, gene_names=("GENE1",)))
    joint = base.merge(Proteome.from_proteins((rna, dna), reference=reference))
    assert len(joint) == 4 and len(joint.by_sequence()["MPEPTIDE"]) == 3
    assert rna.occurrences == dna.occurrences == ()
    assert rna.producer.start_complete is False and rna.producer.stop_complete is None
    assert rna.gene_name == "GENE1" and dna.gene_name == "GENE1"
    assert joint.with_proteins((rna, dna)).identity == joint.identity
    path = tmp_path / "joint.proteome.json.gz"
    joint.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == joint.identity
    assert loaded.by_accession()[rna.accession].producer == rna.producer


def test_source_qualified_ids_distinguish_samples_tools_and_runs(reference):
    first = produced(reference=reference)
    second = produced(reference=reference, sample_id="tumor-2")
    other_tool = ProteinSequence.from_producer("MPEPTIDE", taxon_id=9606, producer=ProducerInfo(
        origin="nc_orf", method="rna_assembly", tool="other", native_id="ASM42",
        sample_id="tumor-1", reference=reference))
    other_run = produced(reference=reference, run_id="run-2")
    assert len({item.accession for item in (first, second, other_tool, other_run)}) == 4
    combined = Proteome.from_proteins((first, second, other_tool, other_run), reference=reference)
    assert len(combined) == 4
    assert combined.with_proteins(first).identity == combined.identity
    with pytest.raises(ValueError, match="Conflicting protein identity"):
        combined.with_proteins(produced("MCHANGED", reference=reference))


def test_unknown_reference_requires_an_explicit_policy(reference, tmp_path):
    unknown = produced()
    with pytest.raises(ValueError, match="unscoped proteome"):
        Proteome.from_proteins([unknown])
    mixed = Proteome.from_proteins([unknown], compatibility="mixed")
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
    snapshot = Proteome.from_proteins((fusion, missing), reference=reference)
    assert len(snapshot) == 2


def test_wrong_taxon_or_reference_and_conflicting_merge_fail(reference):
    another = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "114")
    with pytest.raises(ValueError, match="another reference"):
        Proteome.from_proteins([produced(reference=another)], reference=reference)
    with pytest.raises(ValueError, match="another species"):
        Proteome.from_proteins([produced(reference=None, taxon_id=10090)],
                              reference=reference, compatibility="allow_unknown")
    first = Proteome.from_proteins([produced(reference=reference)], reference=reference)
    changed = replace(first.proteins[0], sequence="MCHANGED")
    with pytest.raises(ValueError, match="Conflicting protein identity"):
        first.merge(Proteome.from_proteins([changed], reference=reference))
    with pytest.raises(ValueError, match="different reference"):
        first.merge(Proteome.from_proteins([produced(reference=another)], reference=another))


def test_legacy_v2_snapshot_still_loads(reference, tmp_path):
    ordinary = Proteome(reference, [ProteinSequence("ENSP1", "MAA", "annotated", 9606,
                                                 "Ensembl:115")])
    path = tmp_path / "legacy.json.gz"
    ordinary.save(path)
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    snapshot["schema"] = "orfdata.proteome.v2"
    path.write_bytes(gzip.compress(json.dumps(snapshot).encode(), mtime=0))
    assert Proteome.load(path).identity == ordinary.identity


def test_external_viral_origin_keeps_the_host_boundary(reference):
    source_reference = ReferenceContext(999999, "synthetic_virus", "toy-v1", "source", "1")
    entry = ProteinSequence.from_producer("MAA", taxon_id=999999, producer=ProducerInfo(
        origin="viral", method="rna_assembly", tool="fixture", native_id="A",
        reference=source_reference, host_taxon_ids=(9606,), host_source="synthetic fixture"))
    host = Proteome.from_proteins([entry], reference=reference)
    assert len(host.host()) == 0
    with pytest.raises(ValueError, match="viral parents"):
        host.with_derived_proteins(DerivedProtein(
            "changed", "MAB", (entry.accession,), "sample-1", ("change-1",), "fixture"))
    mouse = ReferenceContext(10090, "mouse", "GRCm39", "Ensembl", "115")
    with pytest.raises(ValueError, match="recorded host"):
        Proteome.from_proteins([entry], reference=mouse)


def test_export_paths_merge_without_splitting_source_identity(reference, tmp_path):
    def exported(exporter, *, sequence="MPEPTIDE", sample="tumor-1", read_set="bam-1"):
        protein = ProteinSequence.from_producer(
            sequence,
            taxon_id=9606,
            producer=ProducerInfo(
                origin="nc_orf",
                method="rna_assembly",
                tool="isovar",
                native_id="hypothesis-1",
                sample_id=sample,
                reference=reference,
                exporter=exporter,
            ),
        )
        evidence = ProducerRecord(
            namespace="isovar",
            schema_version="isovar.protein_hypotheses.v2",
            kind="evidence",
            native_id="read-set-1",
            sample_id=sample,
            exporter=exporter,
            support="rna_sequence",
            scope={"read_set": read_set},
        )
        return Proteome.from_proteins(
            protein, reference=reference, records=(evidence,)
        )

    direct = exported("isovar")
    wrapped = exported("vaxrank")
    assert direct.proteins[0].accession == wrapped.proteins[0].accession
    assert direct.records[0].key == wrapped.records[0].key
    joint = direct.merge(wrapped)
    assert len(joint.proteins) == len(joint.records) == 1
    assert joint.proteins[0].producer.exporters == ("isovar", "vaxrank")
    assert joint.proteins[0].producer.exporter is None
    assert joint.records[0].exporters == ("isovar", "vaxrank")
    assert joint.records[0].exporter is None
    assert joint.merge(direct, wrapped).identity == joint.identity
    path = tmp_path / "wrapped.json.gz"
    joint.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == joint.identity
    assert loaded.proteins[0].producer.exporters == ("isovar", "vaxrank")
    assert loaded.records[0].exporters == ("isovar", "vaxrank")
    assert loaded.fasta_index() == joint.fasta_index()

    assert len(direct.merge(exported("vaxrank", sample="tumor-2")).proteins) == 2
    with pytest.raises(ValueError, match="Conflicting protein identity"):
        direct.merge(exported("vaxrank", sequence="MCHANGED"))
    with pytest.raises(ValueError, match="Conflicting producer record identity"):
        direct.merge(exported("vaxrank", read_set="bam-2"))


def test_unknown_source_taxon_survives_and_requires_explicit_admission(reference, tmp_path):
    unknown = produced(taxon_id=None, origin="unknown", method="caller_supplied",
                       tool="lens", native_id="row-1", sample_id=None)
    standalone = Proteome.from_proteins(unknown, compatibility="mixed")
    assert standalone.proteins[0].taxon_id is None
    assert "|taxon=unknown|" in standalone.proteins[0].fasta_token()
    path = tmp_path / "unknown-taxon.json.gz"
    standalone.save(path)
    reloaded = Proteome.load(path)
    assert reloaded.proteins[0].taxon_id is None
    assert reloaded.reference is None

    with pytest.raises(ValueError, match="explicit compatibility"):
        Proteome(reference, ()).merge(standalone)
    admitted = Proteome(reference, ()).merge(standalone, compatibility="allow_unknown")
    assert admitted.reference_policy == "allow_unknown"
    assert admitted.proteins[0].taxon_id is None
    assert admitted._content()["reference_policy"] == "allow_unknown"

    dog = produced(taxon_id=9615, origin="unknown", method="caller_supplied",
                   tool="lens", native_id="dog-row", sample_id=None)
    dog_only = Proteome.from_proteins(dog, compatibility="mixed")
    with pytest.raises(ValueError, match="another species"):
        Proteome(reference, ()).merge(dog_only, compatibility="allow_unknown")


def test_v4_exporter_specific_snapshot_is_normalized(reference, tmp_path):
    entry = produced(reference=reference, exporter="vaxrank")
    record = ProducerRecord(
        "isovar", "v2", "translation", "T1",
        exporter="vaxrank", sample_id="tumor-1",
    )
    collection = Proteome.from_proteins(
        entry,
        reference=reference,
        records=(record,),
        links=(ProducerLink(record.key, entry.accession, "protein", "translates_to"),),
    )
    path = tmp_path / "legacy-v4.json.gz"
    collection.save(path)
    snapshot = json.loads(gzip.decompress(path.read_bytes()))
    content = snapshot["content"]
    old_accession = entry.accession.replace("|sample=", "|exporter=vaxrank|sample=")
    content["entries"][0]["accession"] = old_accession
    producer = content["entries"][0]["producer"]
    producer.pop("exporters")
    producer["exporter"] = "vaxrank"
    old_record_key = record.key.replace("|sample=", "|exporter=vaxrank|sample=")
    saved_record = content["producer_graph"]["records"][0]
    saved_record.pop("exporters")
    saved_record["exporter"] = "vaxrank"
    saved_link = content["producer_graph"]["links"][0]
    saved_link.update(source=old_record_key, target=old_accession, target_kind="entry")
    snapshot.update(schema="orfdata.proteome.v4", identity=_sha256(content))
    path.write_bytes(gzip.compress(json.dumps(snapshot).encode(), mtime=0))

    loaded = Proteome.load(path)
    assert loaded.proteins[0].accession == entry.accession
    assert loaded.proteins[0].producer.exporters == ("vaxrank",)
    assert loaded.records[0].key == record.key
    assert loaded.records[0].exporters == ("vaxrank",)
    assert loaded.links == (ProducerLink(record.key, entry.accession,
                                         "protein", "translates_to"),)
    current = tmp_path / "current.json.gz"
    loaded.save(current)
    assert Proteome.load(current).identity == loaded.identity
