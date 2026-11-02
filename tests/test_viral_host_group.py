"""What each source record says about its host and lineage reaches the reference."""

import csv

from conftest import serve_records, viral_documents

from orfdata import download_viral_references, viral_references
from orfdata.genbank import parse
from orfdata.reference import species_taxon
from orfdata.viral_annotation import genus, source_host


def test_the_record_parser_reads_lineage_and_isolation_host():
    from conftest import _genbank

    record = parse(_genbank("VTEST009.1", "ATG" + "GCT" * 8 + "TAA", [], name="VTEST009",
                            host="Homo sapiens", lineage="Viruses; Fooviridae; Barvirus"))
    assert record.lineage == ("Viruses", "Fooviridae", "Barvirus")
    assert genus(record.lineage) == "Barvirus"
    # A lineage that ends in a species group still gives its genus; none gives none.
    assert genus(("Viruses", "Papillomaviridae", "Alphapapillomavirus", "Alphapapillomavirus 9")) \
        == "Alphapapillomavirus"
    assert genus(("Viruses", "Papillomaviridae")) == ""
    assert source_host(record) == "Homo sapiens"
    bare = parse(_genbank("VTEST010.1", "ATG" + "GCT" * 8 + "TAA", [], name="VTEST010"))
    assert bare.lineage == () and source_host(bare) == ""


def test_the_genome_table_keeps_host_and_group_as_the_record_gave_them(
        viral_panel, tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents())
    bundle = download_viral_references(species="human", cache_root=str(tmp_path))
    rows = {row["accession_version"]: row for row in
            csv.DictReader(open(bundle.genome_table), delimiter="\t")}
    human, chicken, silent = rows["VTEST001.1"], rows["VTEST002.1"], rows["VTEST003.1"]
    assert (human["group"], human["source_host"], human["source_host_taxon_ids"]) == (
        "Alphatestvirus", "Homo sapiens", "9606")
    # A host ORFdata has no taxon for stays as text, with no taxon invented.
    # Notes after a semicolon, as GenBank allows, do not hide the species.
    assert species_taxon("Homo sapiens; child".split(";")[0]) == 9606
    assert (chicken["source_host"], chicken["source_host_taxon_ids"]) == ("Gallus gallus", "")
    # A record that says nothing gives nothing, whatever the collection's target is.
    assert (silent["group"], silent["source_host"], silent["source_host_taxon_ids"]) == ("", "", "")


def test_installed_references_select_by_host_and_group(viral_panel, tmp_path, monkeypatch):
    serve_records(monkeypatch, viral_documents())
    download_viral_references(species="human", cache_root=str(tmp_path))
    found = {item.source_accession: item for item in viral_references(cache_root=str(tmp_path))}
    human = found["VTEST001.1"]
    assert (human.group, human.host_taxon_ids) == ("Alphatestvirus", (9606,))
    assert "VTEST001.1" in human.host_source and "Homo sapiens" in human.host_source
    assert not found["VTEST003.1"].host_taxon_ids and found["VTEST003.1"].host_source == ""
    chosen = viral_references(host="human", cache_root=str(tmp_path))
    assert [item.source_accession for item in chosen if item.bundle] == ["VTEST001.1"]
    assert [item.source_accession for item in viral_references(
        group="betatestvirus", cache_root=str(tmp_path))] == ["VTEST002.1"]
