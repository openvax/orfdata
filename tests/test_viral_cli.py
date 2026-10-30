"""The viral CLI makes current and missing GDC coverage visible offline."""

import json

import pytest

from orfdata.cli import run


def test_available_table_and_missing_search(capsys):
    run(["viral", "coverage", "--available"])
    output = capsys.readouterr().out
    assert "1/200" in output and "HPV16" in output
    assert "HPV18" not in output

    run(["viral", "coverage", "--missing", "--virus", "HPV18"])
    output = capsys.readouterr().out
    assert "HPV18" in output and "not bundled" in output
    assert "HPV16" not in output


def test_json_and_tsv_are_machine_readable(capsys):
    run(["viral", "coverage", "--available", "--format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert (payload["total"], payload["available"], payload["shown"]) == (200, 1, 1)
    assert payload["entries"][0]["bundled_reference_id"] == "GRCh38.d1.vd1:HPV16"

    run(["viral", "coverage", "--missing", "--search", "NC_006273.2", "--format", "tsv"])
    output = capsys.readouterr().out.splitlines()
    assert len(output) == 2 and "abbreviation\tname\tstatus" in output[0]
    assert "CMV ( HHV-5)" in output[1] and "not_bundled" in output[1]


def test_inventory_download_routes_to_api(monkeypatch, capsys):
    import orfdata.viral_cli as viral_cli

    called = {}

    def fake_download(**options):
        called.update(options)
        return "/cache/gdc-inventory.txt"

    monkeypatch.setattr(viral_cli, "download_viral_inventory", fake_download)
    run(["viral", "inventory", "download", "--cache-root", "/cache", "--force"])
    assert capsys.readouterr().out.strip() == "/cache/gdc-inventory.txt"
    assert called == {"cache_root": "/cache", "force": True, "timeout": 60}

    with pytest.raises(SystemExit):
        run(["viral", "coverage", "--available", "--missing"])


def test_subsets_listing_reports_resolved_reference_counts(capsys):
    run(["viral", "subsets"])
    output = capsys.readouterr().out
    assert "SUBSET" in output and "cancer-causing" in output
    assert "20" in output and "IARC" in output

    run(["viral", "subsets", "--json"])
    rows = json.loads(capsys.readouterr().out)
    assert {row["name"]: row["references"] for row in rows} == {
        "cancer-associated": 37, "cancer-causing": 20, "cancer-integrated-driver": 15}


def test_contexts_listing_is_readable_and_machine_readable(capsys):
    run(["viral", "contexts"])
    assert "cancer-transcriptome-screening" in capsys.readouterr().out
    run(["viral", "contexts", "--json"])
    rows = json.loads(capsys.readouterr().out)
    assert rows[0]["references"] == 200 and "cancer" in rows[0]["aliases"]


def test_summary_describes_the_installed_layers(viral_panel, tmp_path, capsys):
    from orfdata import download_viral_references

    download_viral_references(species="human", cache_root=str(tmp_path))
    run(["viral", "summary", "--cache-root", str(tmp_path)])
    output = capsys.readouterr().out
    assert "gdc-viruses (GRCh38.d1.vd1) for homo_sapiens taxon 9606" in output
    assert "selection: whole collection" in output
    for layer in ("genome", "transcriptome", "proteome"):
        assert layer in output
    assert "annotation:" in output and "annotation.gtf" in output
    assert "coding transcripts have whole start and stop codons" in output

    run(["viral", "summary", "--context", "cancer", "--cache-root", str(tmp_path), "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["selection"]["contexts"] == ["cancer-transcriptome-screening"]
    assert (payload["genome_count"], payload["transcript_count"]) == (2, 3)
    assert [layer["name"] for layer in payload["layers"]] == ["genome", "transcriptome", "proteome"]
    assert payload["layers"][0]["records"] == 2


def test_summary_without_an_installed_collection_says_how_to_get_one(tmp_path, capsys):
    with pytest.raises(SystemExit):
        run(["viral", "summary", "--cache-root", str(tmp_path)])
    assert "orfdata download viral" in capsys.readouterr().err


def test_download_viral_builds_the_collection_then_a_selection(viral_panel, tmp_path, capsys):
    run(["download", "viral", "--cache-root", str(tmp_path), "--json"])
    collection = json.loads(capsys.readouterr().out)
    assert collection["kind"] == "viral" and collection["genome_count"] == 4
    assert [layer["name"] for layer in collection["layers"]] == [
        "genome", "transcriptome", "proteome"]

    run(["download", "viral", "--cache-root", str(tmp_path), "--context", "cancer", "--json"])
    selection = json.loads(capsys.readouterr().out)
    assert selection["kind"] == "viral-selection" and selection["genome_count"] == 2
    assert selection["source_collection"] == collection["path"]

    # Excluding a context without naming one starts from every context.
    run(["download", "viral", "--cache-root", str(tmp_path), "--exclude-context",
         "cancer-transcriptome-screening", "--json"])
    remainder = json.loads(capsys.readouterr().out)
    assert remainder["genome_count"] == 2 and remainder["contexts"] == []


def test_summary_reports_genomes_without_annotated_transcripts(viral_panel, tmp_path, capsys):
    from orfdata import download_viral_references

    download_viral_references(species="human", cache_root=str(tmp_path))
    run(["viral", "summary", "--cache-root", str(tmp_path)])
    output = capsys.readouterr().out
    assert "1 of 4 genomes have no annotated CDS feature" in output
    assert "unannotated genomes: TV4" in output
    assert "3 of 4 coding transcripts have whole start and stop codons" in output


def test_completeness_surveys_each_installed_genome(viral_panel, tmp_path, capsys):
    from orfdata import download_viral_references

    download_viral_references(species="human", cache_root=str(tmp_path))
    run(["viral", "completeness", "--cache-root", str(tmp_path)])
    output = capsys.readouterr().out
    assert "4 genomes" in output and "fully annotated:" in output
    assert "initiation: " in output and "atg " in output
    for column in ("VIRUS", "CDSS", "COMPLETE", "INCOMPLETE", "UNNAMED", "SPLICED"):
        assert column in output
    for abbreviation in ("TV1", "TV2", "TV3", "TV4"):
        assert abbreviation in output

    run(["viral", "completeness", "--cache-root", str(tmp_path), "--incomplete"])
    output = capsys.readouterr().out
    assert "TV2" in output and "TV4" in output  # One has no closed frame, one has no CDS.
    assert "TV1" not in output.split("VIRUS", 1)[1]

    run(["viral", "completeness", "--cache-root", str(tmp_path), "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["totals"]["genomes"] == 4
    assert {row["abbreviation"] for row in payload["genomes"]} == {"TV1", "TV2", "TV3", "TV4"}
    assert all("shortfall" in row for row in payload["genomes"])


def test_completeness_without_an_installed_collection_says_how_to_get_one(tmp_path, capsys):
    with pytest.raises(SystemExit):
        run(["viral", "completeness", "--cache-root", str(tmp_path)])
    assert "orfdata download viral" in capsys.readouterr().err
