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
    ReferenceContext,
    RiboSeq,
    extra_orfs,
    load_catalog,
    merge_catalogs,
)


def other_catalog():
    """A second source: one ORF matching ms-orf, one reusing the ID ribo-orf, one under another code."""
    peptide = PeptideEvidence("immunopeptidomics", "KPG", ((2, 4),), "other", "fixture")
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


def test_direct_catalogue_preserves_colliding_source_records_and_lookup(tmp_path):
    first = replace(orf("shared", ribosome_studies=("study-a",)), source="a")
    second = replace(orf("shared", ribosome_studies=("study-b",)), source="b")
    catalog = ORFCatalog([first, second], metadata("combined-v1", "combined"))
    assert [r.orf_id for r in catalog.records] == ["a:shared", "b:shared"]
    for original in (first, second):
        record = catalog.get("shared", source=original.source)
        assert record.members == (original,)
        assert catalog.get(f"{original.source}:shared") is record
        assert record.ribosome_studies == original.ribosome_studies
    with pytest.raises(ValueError, match="shared names several ORFs"):
        catalog.get("shared")
    with pytest.raises(KeyError):
        catalog.get("shared", source="absent")
    path = tmp_path / "collisions.json.gz"
    catalog.save(path)
    loaded = load_catalog(path)
    assert loaded.to_dict() == catalog.to_dict()
    assert loaded.get("shared", source="a").members == (first,)
    with pytest.raises(ValueError, match="names several ORFs"):
        loaded.get("shared")


def test_source_qualified_ids_cannot_overwrite_native_ids_or_delimiter_collisions():
    records = [replace(orf(identifier), source=source) for source, identifier in (
        ("a", "x"), ("b", "x"), ("a", "a:x"), ("a:b", "c"), ("a", "b:c"))]
    catalog = ORFCatalog(records, metadata("combined-v1", "combined"))
    assert len({r.orf_id for r in catalog.records}) == len(records)
    for original in records:
        resolved = catalog.get(original.orf_id, source=original.source)
        assert original in (resolved.members or (resolved,))
    for ambiguous in ("x", "a:x", "a:b:c"):
        with pytest.raises(ValueError, match="names several ORFs"):
            catalog.get(ambiguous)


@pytest.mark.parametrize("changed", [False, True])
def test_duplicate_ids_within_one_source_still_raise(changed):
    original = orf("same")
    other = replace(original, protein_sequence="MKPA") if changed else original
    with pytest.raises(ValueError, match="Duplicate ORF identifiers"):
        ORFCatalog([original, other], metadata("fixture-v1", "synthetic"))


def test_direct_collisions_keep_distinct_placements_and_proteome_records(genome):  # noqa: F811
    records = [replace(orf("shared"), source="a"),
               replace(orf("shared", blocks=((218, 226),), protein="MK"), source="b")]
    catalog = ORFCatalog(records, metadata("combined-v1", "combined"))
    found = extra_orfs(genome, catalog=catalog)
    assert {o.protein_sequence for o in found.for_orf("shared", source="a")} == {"MKPG"}
    assert {o.protein_sequence for o in found.for_orf("shared", source="b")} == {"MK"}
    with pytest.raises(ValueError, match="names several ORFs"):
        found.for_orf("shared")
    proteins = Proteome.from_genome(genome).with_extra_orfs(None, catalog=catalog)
    assert {(p.accession, p.sequence) for p in proteins if p.kind == "extra_orf"} == {
        ("a:shared", "MKPG"), ("b:shared", "MK")}


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
        paths += ["--nc-orf-sources", str(tmp_path / f"{name}.json")]
        (tmp_path / f"{name}.json").write_text(json.dumps(catalog.to_dict()))
    run(["proteome", "--species", "mouse", "--release", "1", *paths])
    assert "equivalent=other:OP_1" in capsys.readouterr().out


def single(source, *records, version="v1"):
    return ORFCatalog([replace(r, source=source) for r in records], dict(
        genome_build="GRCm39", taxon_id=10090, species="mus_musculus", coordinates="1-based-inclusive",
        dataset_version=f"{source}-{version}", sources=[dict(id=source, scope="test fixture")]))


def test_review_cases_for_ids_contexts_and_re_merging(genome):  # noqa: F811
    base = make_catalog()
    ms = base.get("ms-orf")
    # Two later sources call the same ORF IP_1: one ORF, not an ambiguity.
    b, c = single("b", replace(ms, orf_id="IP_1", initiation_context=None)), single("c", replace(ms, orf_id="IP_1"))
    merged = merge_catalogs([base, b, c])
    assert merged.get("IP_1") is merged.get("ms-orf") and merged.get("c:IP_1") is merged.get("ms-orf")
    assert merged.get("synthetic:absent-orf").orf_id == "absent-orf"  # Unmerged records answer to source:id too.
    # A merged catalogue merges again from its source records.
    again = merge_catalogs([merge_catalogs([base, b]), c])
    assert again.records == merged.records
    # Catalogue order does not decide which initiation context an ORF keeps.
    flipped = merge_catalogs([b, base])
    assert flipped.get("IP_1").initiation_context == ms.initiation_context
    # A member's evidence cannot be replaced by other evidence.
    with pytest.raises(ValueError, match="members' evidence"):
        replace(merged.get("ms-orf"), ribosome_studies=("invented",))
    assert merged.metadata["merged_from"][0]["metadata"]["dataset_version"] == "fixture-v1"


def test_pins_strains_conflicts_and_transcripts():
    base = make_catalog()
    ms = base.get("ms-orf")
    older = ReferenceContext(10090, "mus_musculus", "GRCm39", "synthetic", "0")
    pinned = single("b", replace(ms, orf_id="IP_1", reference_context=older))
    # Sources pinned to different annotations leave the merged ORF unpinned, in either order.
    for order in ([base, pinned], [pinned, base]):
        assert merge_catalogs(order).get("IP_1").reference_context is None
    strain = replace(older, strain="C57BL/6J")
    with pytest.raises(ValueError, match="different or unknown strains"):
        merge_catalogs([base, single("b", replace(ms, orf_id="IP_1", reference_context=strain))])
    merged = merge_catalogs([base, other_catalog()])
    # One path, two genetic codes; and the fixture's own changed-orf, whose protein differs from ms-orf's.
    assert merged.metadata["conflicts"] == [["OP_3", "synthetic:ribo-orf"], ["absent-orf", "changed-orf", "ms-orf"]]
    # absent-orf shares ms-orf's path on another transcript: a source record on that transcript joins it.
    late = single("late", replace(ms, orf_id="L_1", transcript_id="missing_rna", initiation_context=None))
    assert merge_catalogs([base, late]).get("L_1").orf_id == "absent-orf"


def metadata(version, source, **extra):
    return dict(genome_build="GRCm39", taxon_id=10090, species="mus_musculus", coordinates="1-based-inclusive",
                dataset_version=version, sources=[dict(id=source, scope="test fixture")], **extra)


def test_one_study_named_twice_is_counted_once_by_its_pmid():
    first = ORFCatalog([replace(orf("A1", ribosome_studies=("Chen_et_al_2020",)), source="a")],
                       metadata("a-v1", "a", ribosome_study_pmids={"Chen_et_al_2020": "32139545"}))
    second = ORFCatalog([replace(orf("B1", ribosome_studies=("Chen2020",)), source="b")],
                        metadata("b-v1", "b", ribosome_study_pmids={"Chen2020": "32139545"}))
    merged = merge_catalogs([first, second])
    assert merged.get("A1").ribosome_studies == ("Chen_et_al_2020",)
    assert not merged.get("A1").has_evidence(RiboSeq(min_studies=2))
    assert merged.metadata["ribosome_study_pmids"] == {"Chen_et_al_2020": "32139545", "Chen2020": "32139545"}
    # Without a PMID to tell, two labels stay two studies.
    unmapped = ORFCatalog([replace(orf("C1", ribosome_studies=("Chen2020",)), source="c")], metadata("c-v1", "c"))
    assert merge_catalogs([first, unmapped]).get("A1").has_evidence(RiboSeq(min_studies=2))
    # A merged record may list fewer studies than its members, never others.
    with pytest.raises(ValueError, match="must be its members' evidence"):
        replace(merged.get("A1"), ribosome_studies=("Elsewhere",))


def test_a_redundant_start_codon_and_a_shared_transcript_join_one_orf(genome):  # noqa: F811
    first = ORFCatalog([replace(orf("A1"), source="a")], metadata("a-v1", "a"))
    declared = ORFCatalog([replace(orf("C1"), source="c", start_codon="ATG")], metadata("c-v1", "c"))
    merged = merge_catalogs([first, declared])
    assert [r.orf_id for r in merged.records] == ["A1"] and merged.metadata["conflicts"] == []
    # B's record on A1's transcript joins A1, although B lists another transcript's record first.
    second = ORFCatalog([replace(orf("B_other", transcript_id="test_rna_long"), source="b"),
                         replace(orf("B_same"), source="b")], metadata("b-v1", "b"))
    merged = merge_catalogs([first, second])
    assert [m.orf_id for m in merged.get("A1").members] == ["A1", "B_same"]
    assert merged.get("B_other").members == ()
    # Every source's transcript and gene find the merged ORF.
    labelled = ORFCatalog([replace(orf("Y", transcript_id="test_rna_long"), source="y", gene_name="Other")],
                          metadata("y-v1", "y"))
    merged = merge_catalogs([first, labelled])
    assert [r.orf_id for r in merged.query(transcript_id="test_rna_long")] == ["A1"]
    assert [r.orf_id for r in merged.query(gene="Other")] == ["A1"]
    long = genome.transcript_by_id("test_rna_long")
    assert [o.orf_id for o in merged.for_transcript(long)] == ["A1"]
    found = extra_orfs(genome, None, catalog=merged)
    assert all(o.on_source_transcript for o in found.for_orf("A1")) and found.source_failures == ()
    assert [o.orf_id for o in found.for_gene("Other", source=True)] == ["A1", "A1"]


def test_records_of_one_catalogue_stay_as_given_whatever_else_is_merged():
    both = ORFCatalog([replace(orf("P"), source="s1"), replace(orf("Q"), source="s2")], metadata("x-v1", "s1"))
    unrelated = ORFCatalog([replace(orf("Z", blocks=((218, 226),), protein="MK"), source="z")], metadata("z-v1", "z"))
    assert sorted(r.orf_id for r in merge_catalogs([both, unrelated]).records) == ["P", "Q", "Z"]


def test_the_cli_finds_records_by_any_source_id_and_merges_repeated_catalogues(tmp_path, capsys):
    from orfdata.cli import run

    paths = []
    for name, catalog in (("a", make_catalog()), ("b", other_catalog())):
        paths += ["--catalog", str(tmp_path / f"{name}.json")]
        (tmp_path / f"{name}.json").write_text(json.dumps(catalog.to_dict()))
    run(paths + ["--orf-id", "other:OP_1", "--format", "json"])
    assert [r["orf_id"] for r in json.loads(capsys.readouterr().out)["records"]] == ["ms-orf"]


def test_a_biotype_query_finds_a_merged_orf_by_any_source_label():
    first, second = make_catalog(), other_catalog()
    second = ORFCatalog([replace(r, biotype="lncORF") if r.orf_id == "OP_1" else r for r in second.records],
                        second.metadata)
    merged = merge_catalogs([first, second])
    assert merged.get("ms-orf") in merged.query(biotype="lncORF")
    assert merged.get("ms-orf") in merged.query(biotype=merged.get("ms-orf").biotype)
