"""The pinned source-record archive is listed, described and installed like any other data."""

import hashlib
import io
import json
import tarfile
from dataclasses import replace

import pytest

from orfdata import source_records
from orfdata.downloads import list_data
from orfdata.errors import NotInstalledError

NAME = "gdc-viruses-records"


def test_the_archive_is_listed_beside_other_data_and_missing_until_fetched(tmp_path):
    (row,) = list_data(kind="records", cache_root=str(tmp_path))
    assert (row.kind, row.name, row.version, row.status, row.acquire) == \
        ("records", NAME, "1", "missing", "download")
    assert row.name in {item.name for item in list_data(cache_root=str(tmp_path))}


def test_the_pin_names_a_release_asset_by_checksum_and_size():
    (pin,) = source_records.pinned_records()
    assert pin.url.startswith("https://github.com/openvax/orfdata/releases/download/data-")
    assert len(pin.sha256) == 64 and pin.size > 0 and pin.filename.endswith(".tar.gz")


def test_nothing_is_read_before_the_archive_is_installed(tmp_path):
    with pytest.raises(NotInstalledError, match="orfdata data fetch"):
        source_records.read_records(NAME, cache_root=str(tmp_path))


def test_an_installed_archive_is_checked_against_its_pin_and_read_offline(tmp_path, monkeypatch):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        data = b"LOCUS       X\n"
        member = tarfile.TarInfo("genbank__X1.1.gb.gz")
        member.size = len(data)
        archive.addfile(member, io.BytesIO(data))
    raw = buffer.getvalue()
    local = tmp_path / "records.tar.gz"
    local.write_bytes(raw)
    (pin,) = source_records.pinned_records()
    small = replace(pin, url=pin.url.rsplit("/", 1)[0] + "/records.tar.gz",
                    sha256=hashlib.sha256(raw).hexdigest(), size=len(raw))
    monkeypatch.setattr(source_records, "pinned_records", lambda: (small,))

    cache = str(tmp_path / "cache")
    source_records.install_records(NAME, cache_root=cache, source_path=local)
    assert source_records.read_records(NAME, cache_root=cache) == {"genbank__X1.1.gb.gz": data}
    assert list_data(kind="records", cache_root=cache)[0].status == "available"

    wrong = tmp_path / "wrong.tar.gz"
    wrong.write_bytes(raw + b"x")
    with pytest.raises(Exception, match="(?i)sha|size|checksum|match"):
        source_records.install_records(NAME, cache_root=str(tmp_path / "other"), source_path=wrong)


def test_the_cli_describes_the_archive_offline(tmp_path, capsys):
    from orfdata.cli import run

    assert run(["data", "info", NAME, "--cache-root", str(tmp_path), "--json"]) in (0, None)
    info = json.loads(capsys.readouterr().out)
    assert info["kind"] == "records" and info["source"]["sha256"] and info["installed"] is False
