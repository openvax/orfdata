"""Stream one member out of a remote ZIP64 archive (Deflate64 or Deflate) with parallel HTTP
range requests consumed in order. Verifies the member's CRC-32 and size, records its SHA-256,
counts MS2 spectra, and writes it gzip-compressed (level 1) to OUT, via OUT.part.

On success prints one JSON receipt line: member, zip_url, compressed_bytes, mzml_bytes,
mzml_crc32, crc_ok, mzml_sha256, gz_bytes, n_ms2, download_s, MBps. On a CRC or size mismatch
renames OUT.part to OUT.part.bad and exits non-zero.

usage: python -I fetch_member.py URL MEMBER OUT.mzML.gz [NCONN]
  NCONN  parallel range requests (default 16); each request is a 16 MiB chunk.
"""
import concurrent.futures as cf
import gzip
import hashlib
import json
import os
import struct
import sys
import time
import urllib.request
import zipfile
import zlib

import inflate64

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remotezip import HTTPRangeFile  # noqa: E402

CHUNK = 16 << 20
MS2 = b'name="ms level" value="2"'


def fetch(url, s, e):
    for attempt in range(10):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e - 1}"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) == e - s:
                return data
            raise IOError(f"short read {len(data)} != {e - s}")
        except Exception as ex:
            print(f"retry {s} {attempt} {ex}", file=sys.stderr, flush=True)
            time.sleep(5 * (attempt + 1))
    raise IOError("range fetch failed")


def main():
    if len(sys.argv) not in (4, 5):
        raise SystemExit(__doc__)
    url, member, out = sys.argv[1:4]
    nconn = int(sys.argv[4]) if len(sys.argv) > 4 else 16

    f = HTTPRangeFile(url)
    zi = zipfile.ZipFile(f).getinfo(member)
    f.seek(zi.header_offset)
    sig, _, _, _, _, _, _, _, _, nlen, elen = struct.unpack("<IHHHHHIIIHH", f.read(30))
    assert sig == 0x04034B50, "not a local file header"
    start = zi.header_offset + 30 + nlen + elen
    end = start + zi.compress_size
    if zi.compress_type == 9:
        decomp = inflate64.Inflater().inflate
    elif zi.compress_type == 8:
        decomp = zlib.decompressobj(-15).decompress
    else:
        raise SystemExit(f"unsupported compression method {zi.compress_type}")

    ranges = [(s, min(s + CHUNK, end)) for s in range(start, end, CHUNK)]
    t0 = time.time()
    crc32 = 0
    sha = hashlib.sha256()
    n_in = n_out = n_ms2 = 0
    tail = b""
    part = out + ".part"
    with cf.ThreadPoolExecutor(nconn) as ex, open(part, "wb") as raw, \
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0) as gz:
        pending = {}
        nxt = 0
        for i in range(len(ranges)):
            while nxt < len(ranges) and nxt < i + 2 * nconn:
                pending[nxt] = ex.submit(fetch, url, *ranges[nxt])
                nxt += 1
            chunk = pending.pop(i).result()
            n_in += len(chunk)
            data = decomp(chunk)
            crc32 = zlib.crc32(data, crc32)
            sha.update(data)
            buf = tail + data
            n_ms2 += buf.count(MS2)
            tail = buf[-(len(MS2) - 1):]
            gz.write(data)
            n_out += len(data)
    dt = time.time() - t0
    ok = (n_out == zi.file_size) and (crc32 == zi.CRC)
    if not ok:
        os.rename(part, part + ".bad")
        raise SystemExit(f"CRC/size mismatch for {member}: {n_out} vs {zi.file_size}")
    os.rename(part, out)
    print(json.dumps(dict(member=member, zip_url=url, compressed_bytes=zi.compress_size, mzml_bytes=n_out,
                          mzml_crc32=f"{zi.CRC:08x}", crc_ok=ok, mzml_sha256=sha.hexdigest(),
                          gz_bytes=os.path.getsize(out), n_ms2=n_ms2, download_s=round(dt, 1),
                          MBps=round(n_in / dt / 1e6, 2))), flush=True)


if __name__ == "__main__":
    main()
