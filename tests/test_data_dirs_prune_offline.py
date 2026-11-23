"""Tests for orfdata data dirs, prune (dry-run vs --yes), and --offline / OPENVAX_OFFLINE."""

import hashlib
import json
import os

import datacache
import pytest

from orfdata import cache, download_cli


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in (
        "ORFDATA_CACHE",
        "OPENVAX_DATA_CACHE",
        "DATACACHE_OFFLINE",
        "OPENVAX_OFFLINE",
    ):
        monkeypatch.delenv(key, raising=False)


def test_data_dirs_readonly_and_json(tmp_path, monkeypatch, capsys):
    shared = tmp_path / "shared"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")

    download_cli.run(["dirs", "--json"], prog="orfdata data")
    payload = json.loads(capsys.readouterr().out)
    assert payload["cache_root"] == str(shared / "orfdata")
    assert payload["origin"] == "shared_env"
    assert payload["detail"] == "OPENVAX_DATA_CACHE"
    assert payload["exists"] is False
    assert payload["offline"] is True
    assert payload["offline_source"] == "OPENVAX_OFFLINE=1"
    assert not shared.exists()

    # Populate an area and verify human-readable output
    area = shared / "orfdata" / "catalogs"
    area.mkdir(parents=True)
    (area / "sample.bin").write_bytes(b"z" * 4096)

    download_cli.run(["dirs"], prog="orfdata data")
    out = capsys.readouterr().out
    assert "OPENVAX_DATA_CACHE" in out
    assert "catalogs" in out
    assert "4.0 KiB" in out
    assert "enabled (OPENVAX_OFFLINE=1)" in out


def test_data_prune_dry_run_and_apply(tmp_path, capsys):
    root = tmp_path / "cache"
    store_dir = cache.get_data_dir(override=root) / "catalogs" / "demo_store"

    src1 = tmp_path / "v1.bin"
    b1 = b"old-generation" * 50
    src1.write_bytes(b1)
    h1 = hashlib.sha256(b1).hexdigest()
    datacache.install_bundle(store_dir, {"f.bin": {"url": src1.as_uri(), "sha256": h1, "size": len(b1)}})

    src2 = tmp_path / "v2.bin"
    b2 = b"new-generation" * 60
    src2.write_bytes(b2)
    h2 = hashlib.sha256(b2).hexdigest()
    datacache.install_bundle(
        store_dir,
        {"f.bin": {"url": src2.as_uri(), "sha256": h2, "size": len(b2)}},
        force=True,
    )

    # Dry-run default removes nothing
    download_cli.run(["prune", "--cache-root", str(root), "--json"], prog="orfdata data")
    dry = json.loads(capsys.readouterr().out)
    assert dry["dry_run"] is True
    assert dry["removed_count"] == 1
    assert dry["reclaimable_bytes"] > 0

    # --yes deletes stale generation while keeping current
    download_cli.run(["prune", "--cache-root", str(root), "--yes", "--json"], prog="orfdata data")
    applied = json.loads(capsys.readouterr().out)
    assert applied["dry_run"] is False
    assert applied["removed_count"] == 1

    # Second prune finds nothing left
    download_cli.run(["prune", "--cache-root", str(root), "--json"], prog="orfdata data")
    after = json.loads(capsys.readouterr().out)
    assert after["removed_count"] == 0


def test_offline_flag_restores_env_and_blocks_remote(tmp_path, monkeypatch, capsys):
    with pytest.raises(SystemExit) as exc:
        download_cli.run(
            ["fetch", "t2t-chm13", "--cache-root", str(tmp_path / "c"), "--offline"],
            prog="orfdata data",
        )
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert "Offline mode is enabled" in err
    assert "DATACACHE_OFFLINE" not in os.environ
