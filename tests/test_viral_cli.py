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
