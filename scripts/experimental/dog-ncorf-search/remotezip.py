"""Read members of a remote ZIP (including ZIP64) with HTTP range requests, without
downloading the archive.

usage:
  python -I remotezip.py list URL              print member, size, compressed size, method
  python -I remotezip.py get URL MEMBER OUTPATH

fetch_member.py streams large members faster, with parallel requests.
"""
import io
import shutil
import sys
import time
import urllib.request
import zipfile


class HTTPRangeFile(io.RawIOBase):
    """A seekable, read-only file over HTTP range requests, cached in BLOCK-byte blocks."""

    def __init__(self, url, block=1 << 20):
        self.url = url
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req) as r:
            self.size = int(r.headers["Content-Length"])
        self.pos = 0
        self.block = block
        self.cache = {}
        self.nreq = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        if whence == 0:
            self.pos = off
        elif whence == 1:
            self.pos += off
        else:
            self.pos = self.size + off
        return self.pos

    def _fetch(self, start, end):
        for attempt in range(5):
            try:
                req = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end - 1}"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    self.nreq += 1
                    return r.read()
            except Exception as e:  # Retry transient errors.
                print("retry", attempt, e, file=sys.stderr)
                time.sleep(2 * (attempt + 1))
        raise IOError("range fetch failed")

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        # Large reads go straight through.
        if n > 4 * self.block:
            data = self._fetch(self.pos, self.pos + n)
            self.pos += len(data)
            return data
        out = bytearray()
        while n > 0:
            b = self.pos // self.block
            if b not in self.cache:
                if len(self.cache) > 64:
                    self.cache.clear()
                s = b * self.block
                self.cache[b] = self._fetch(s, min(s + self.block, self.size))
            chunk = self.cache[b][self.pos - b * self.block:][:n]
            out += chunk
            self.pos += len(chunk)
            n -= len(chunk)
        return bytes(out)

    def readinto(self, buf):
        data = self.read(len(buf))
        buf[: len(data)] = data
        return len(data)


def listing(url):
    """(member, uncompressed bytes, compressed bytes, compression method) for each member."""
    z = zipfile.ZipFile(HTTPRangeFile(url))
    return [(zi.filename, zi.file_size, zi.compress_size, zi.compress_type) for zi in z.infolist()]


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    cmd, url = sys.argv[1], sys.argv[2]
    if cmd == "list":
        for row in listing(url):
            print("\t".join(map(str, row)))
    elif cmd == "get":
        member, outpath = sys.argv[3], sys.argv[4]
        f = HTTPRangeFile(url)
        z = zipfile.ZipFile(f)
        zi = z.getinfo(member)
        f.block = 16 << 20  # Bigger blocks for streaming.
        t = time.time()
        with z.open(zi) as src, open(outpath, "wb") as dst:
            shutil.copyfileobj(src, dst, 16 << 20)
        print(f"extracted {member} {zi.file_size} bytes (compressed {zi.compress_size}) "
              f"in {time.time() - t:.0f}s, {f.nreq} requests", file=sys.stderr)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
