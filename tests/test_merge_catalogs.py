"""Merging ORF catalogues for one reference without losing source records or evidence (#71)."""

import json
from dataclasses import replace

import pytest
from test_proteome import genome, make_catalog, orf  # noqa: F401

from orfdata import (
    MS,
    ORFCatalog,
    PeptideEvidence,
    Proteome,
    RiboSeq,
    extra_orfs,
    load_catalog,
    merge_catalogs,
)


def other_catalog():
    """A second source: one ORF matching ms-orf, one reusing the ID ribo-orf, one under another code."""
    peptide = PeptideEvidence("mhc", "KPG", ((2, 4),), "other", "fixture")
    records = [
        replace(orf("OP_1", peptides=(peptide,), ribosome_studies=("study-b",)), source="other"),
        replace(orf("ribo-orf", blocks=((206, 214),), protein="KP"), source="other"),
        replace(orf("OP_3", blocks=((218, 226),), protein="MK"), source="other", translation_table=2),
    ]
    return ORFCatalog(records, dict(genome_build="GRCm39", taxon_id=10090, species="mus_musculus",
                                    coordinates="1-based-inclusive", dataset_version="other-v2",
                                    sources=[dict(id="other", scope="test fixture")]))


def test_one_orf_from_two_sources_keeps_both_records_and_their_evidence(tmp_path):
    first, second = make_catalog(), other_catalog()
    merged = merge_catalogs([first, second])
    assert merged.source == "orfdata-catalog:fixture-v1+orfdata-catalog:other-v2"
    ms = merged.get("ms-orf")
    assert ms.source_ids == ("synthetic:ms-orf", "other:OP_1")
    assert [m.orf_id for m in ms.members] == ["ms-orf", "OP_1"]
    assert merged.get("OP_1") is ms and merged.get("other:OP_1") is ms
    assert ms.ribosome_studies == ("study-b",) and [p.source for p in ms.peptides] == ["synthetic", "other"]
    assert ms.has_evidence(RiboSeq() & MS(min_peptides=1)) and not ms.has_evidence(MS(min_peptides=2))
    # A source's own records stay as it gave them: absent-orf shares ms-orf's path but is a separate record.
    assert merged.get("absent-orf").members == ()
    # Different genetic codes on one path are different ORFs.
    assert merged.get("OP_3").translation_table == 2
    assert len(merged.records) == len(first.records) + len(second.records) - 1
    path = tmp_path / "merged.json"
    path.write_text(json.dumps(merged.to_dict()))
    assert load_catalog(path).records == merged.records


def test_one_id_for_different_orfs_is_qualified_and_never_guessed():
    merged = merge_catalogs([make_catalog(), other_catalog()])
    with pytest.raises(ValueError, match="ribo-orf names several ORFs: synthetic:ribo-orf, other:ribo-orf"):
        merged.get("ribo-orf")
    assert merged.get("synthetic:ribo-orf").protein_sequence == "MK"
    assert merged.get("other:ribo-orf").protein_sequence == "KP"
    with pytest.raises(KeyError):
        merged.get("missing")


def test_incompatible_catalogues_are_not_merged():
    first = make_catalog()
    assert merge_catalogs([first]) is first
    with pytest.raises(ValueError, match="merge one release of each source"):
        merge_catalogs([first, make_catalog()])
    with pytest.raises(ValueError, match="Cannot merge a homo_sapiens GRCh38 catalogue"):
        merge_catalogs([first, load_catalog()])
    with pytest.raises(TypeError):
        merge_catalogs([])


def test_merged_catalogues_place_each_orf_once(genome, tmp_path):  # noqa: F811
    catalogs = [make_catalog(), other_catalog()]
    found = extra_orfs(genome, catalog=catalogs)
    placed = [o for o in found if o.orf_id == "ms-orf"]
    assert placed and all(o.equivalent_ids == ("other:OP_1",) for o in placed)
    proteome = Proteome.from_genome(genome).with_extra_orfs(RiboSeq() & MS(), catalog=catalogs)
    entry, = [e for e in proteome if e.kind == "extra_orf"]
    assert (entry.accession, entry.equivalent_ids) == ("ms-orf", ("other:OP_1",))
    assert "equivalent=other:OP_1" in entry.fasta_header()
    assert proteome.layers[-1]["source"] == "orfdata-catalog:fixture-v1+orfdata-catalog:other-v2"
    snapshot = tmp_path / "proteome.json.gz"
    proteome.save(snapshot)
    loaded = Proteome.load(snapshot)
    assert loaded.identity == proteome.identity
    assert [e.equivalent_ids for e in loaded if e.kind == "extra_orf"] == [("other:OP_1",)]


def test_any_source_id_finds_the_merged_orf_and_the_cli_merges_files(genome, tmp_path, capsys, monkeypatch):  # noqa: F811
    from orfdata import proteome_cli
    from orfdata.cli import run

    found = extra_orfs(genome, catalog=[make_catalog(), other_catalog()])
    assert found.for_orf("OP_1") == found.for_orf("other:OP_1") == found.for_orf("ms-orf") != ()
    assert found.for_orf("missing") == ()
    with pytest.raises(ValueError, match="names several ORFs"):
        found.for_orf("ribo-orf")
    monkeypatch.setattr(proteome_cli, "_ensembl_genome", lambda release, species: genome)
    paths = []
    for name, catalog in (("first", make_catalog()), ("second", other_catalog())):
        paths += ["--catalog", str(tmp_path / f"{name}.json")]
        (tmp_path / f"{name}.json").write_text(json.dumps(catalog.to_dict()))
    run(["proteome", "--species", "mouse", "--release", "1", "--min-ribo-seq-studies", "1",
         "--min-ms-peptides", "1", *paths])
    assert "equivalent=other:OP_1" in capsys.readouterr().out
