"""Inspect the pinned public dog result schema with explicit resource budgets."""

import argparse
import json
from pathlib import Path

from tools.ftps_sqlite import ChunkCache, FTPSRanges, Identity, connect, sha256

ASSETS = {
    "pdResultView": (1601896448, "20251024173126"),
    "msfView": (3546578944, "20251024172153"),
    "pdResult": (86023450624, "20251028233329"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("asset", choices=ASSETS)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--transfer-mib", type=int, default=64)
    parser.add_argument("--cache-mib", type=int, default=32)
    parser.add_argument("--seconds", type=int, default=600)
    args = parser.parse_args()
    size, modified = ASSETS[args.asset]
    path = "/v11/MSV000099744/search/search/IQP_1269_Dog_P1toP4." + args.asset
    source = FTPSRanges(Identity("massive-ftp.ucsd.edu", path, size, "213 " + modified),
                        transfer_limit=args.transfer_mib * 1024 * 1024,
                        seconds_limit=args.seconds)
    args.directory.mkdir(parents=True, exist_ok=True)
    cache = ChunkCache(source, args.directory / "chunks", cache_limit=args.cache_mib * 1024 * 1024)
    database, vfs = connect(cache)
    try:
        schema = [dict(zip(("type", "name", "table", "rootpage", "sql"), row))
                  for row in database.execute(
                      "select type,name,tbl_name,rootpage,sql from sqlite_master order by type,name")]
        schema_path = args.directory / "schema.json"
        schema_path.write_text(json.dumps(schema, sort_keys=True, indent=2) + "\n")
        print(json.dumps(dict(schema_objects=len(schema), schema_sha256=sha256(schema_path.read_bytes()),
                              tables=[r["name"] for r in schema if r["type"] == "table"])), flush=True)
    finally:
        database.close()
        vfs.unregister()
        receipt = cache.receipt()
        (args.directory / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
        print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
