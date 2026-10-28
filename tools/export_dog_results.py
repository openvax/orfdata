"""Explicitly export the selected dog vendor facts with resumable bounded reads."""

import argparse
import hashlib
import json
from pathlib import Path

from tools.export_sqlite_rows import CODE_SHA256 as ROW_EXPORT_SHA256
from tools.export_sqlite_rows import NativeExport
from tools.extraction_budget import ScratchBudget
from tools.ftps_sqlite import ChunkCache, FTPSRanges, Identity, connect

PLAN_PATH = Path(__file__).with_name("dog_result_plan.json")
CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def read_ahead_kib(value):
    value = int(value)
    if value < 4 or value > 16384 or value % 4:
        raise argparse.ArgumentTypeError("Read-ahead must be a multiple of 4 KiB between 4 and 16384 KiB")
    return value


def main():
    plan_data = PLAN_PATH.read_bytes()
    plan = json.loads(plan_data)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Existing audit/cache root, or a new owned directory")
    parser.add_argument("--tables", nargs="+", choices=[p["table"] for p in plan["tables"]])
    parser.add_argument("--transfer-mib", type=int, default=2048)
    parser.add_argument("--cache-mib", type=int, default=512)
    parser.add_argument("--scratch-mib", type=int, default=4096)
    parser.add_argument("--output-mib", type=int, default=1024)
    parser.add_argument("--read-ahead-kib", type=read_ahead_kib, default=256)
    parser.add_argument("--lookup-read-ahead-kib", type=read_ahead_kib, default=4,
                        help="Bounded window for indexed acquisitions; default 4 KiB, independent of scan read-ahead")
    parser.add_argument("--requests", type=int, default=30000)
    parser.add_argument("--seconds", type=int, default=1200)
    parser.add_argument("--batch-rows", type=int, default=500)
    parser.add_argument("--lookup-keys", type=int, default=32)
    parser.add_argument("--extend-selection", action="store_true",
                        help="Explicitly add source tables and expand only unstarted acquisition selections")
    args = parser.parse_args()
    root = args.directory.resolve()
    scratch = ScratchBudget(root, args.scratch_mib * 1024 * 1024)
    identity = Identity("massive-ftp.ucsd.edu",
                        "/v11/MSV000099744/search/search/IQP_1269_Dog_P1toP4.pdResult",
                        86023450624, "213 20251028233329")
    source = FTPSRanges(identity, transfer_limit=args.transfer_mib * 1024 * 1024,
                        request_limit=args.requests, seconds_limit=args.seconds, reuse_control=True)
    cache = ChunkCache(source, root / "chunks", cache_limit=args.cache_mib * 1024 * 1024,
                       chunk_size=4096, scratch=scratch, read_ahead=args.read_ahead_kib * 1024)
    database, vfs, export = None, None, None
    invocation = len(list(root.glob("export-receipt-[0-9]*.json"))) + 1
    failure = None
    try:
        database, vfs = connect(cache)
        source_identity = dict(identity=cache.identity, header_sha256=cache.header_sha256)
        export = NativeExport(database, root / "native-evidence.partial.sqlite", plan, source_identity,
                              scratch, output_limit=args.output_mib * 1024 * 1024,
                              batch_rows=args.batch_rows, lookup_keys=args.lookup_keys,
                              extend_selection=args.extend_selection)
        selected = args.tables or list(export.tables)
        for table in selected:
            # Index lookups default to individual pages. A separately bounded
            # window can exploit measured locality without changing selection.
            # Both windows stop before any existing cached segment.
            window = args.lookup_read_ahead_kib if export.tables[table].get("select_from") else args.read_ahead_kib
            cache.read_ahead = window * 1024
            reported_rows = export.progress(table)["rows"]

            def report(name, progress):
                nonlocal reported_rows
                if progress["complete"] or progress["rows"] - reported_rows >= args.batch_rows * 10:
                    print(json.dumps(dict(table=name, **progress)), flush=True)
                    reported_rows = progress["rows"]

            export.copy(table, on_batch=report)
    except BaseException as error:
        failure = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        try:
            receipt = dict(acquisition=cache.receipt(), result=export.summary() if export else None,
                           failure=failure, tool_sha256=CODE_SHA256, row_export_sha256=ROW_EXPORT_SHA256,
                           export_options=dict(batch_rows=args.batch_rows, lookup_keys=args.lookup_keys,
                                               read_ahead_kib=args.read_ahead_kib,
                                               lookup_read_ahead_kib=args.lookup_read_ahead_kib,
                                               batch_bytes=export.batch_bytes if export else None,
                                               requested_tables=args.tables or [p["table"] for p in plan["tables"]]),
                           selection_sha256=hashlib.sha256(plan_data).hexdigest())
            print(json.dumps(dict(failure=failure, result=receipt["result"])), flush=True)
            scratch.atomic_write(root / f"export-receipt-{invocation:04d}.json",
                                 (json.dumps(receipt, sort_keys=True, indent=2) + "\n").encode())
        finally:
            if export is not None:
                export.close()
            if database is not None:
                database.close()
            if vfs is not None:
                vfs.unregister()
            cache.close()


if __name__ == "__main__":
    main()
