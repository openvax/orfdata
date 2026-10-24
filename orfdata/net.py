"""One retrying HTTP read, shared by the source-record downloads."""

import time
from http.client import HTTPException
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

USER_AGENT = "orfdata/viral (https://github.com/openvax/orfdata)"
MAX_ATTEMPTS = 5
# A request the server rejects will be rejected again; only a timeout or a rate
# limit is worth repeating.
_RETRYABLE_STATUS = (408, 429)


class SourceUnavailable(RuntimeError):
    """The source answered, and its answer is that it does not publish this record.

    Telling this apart from every other failure is what lets a ledger record an
    entry as ``unavailable`` rather than guessing from an error message.
    """


def fetch(url, *, timeout, params=None, attempts=MAX_ATTEMPTS, description=None):
    """Read one URL, retrying transport failures with exponential backoff.

    A reset connection or a body that stops early is an ``OSError`` or an
    ``HTTPException``, not only a ``URLError``, and every one of them is worth
    retrying rather than losing the whole download. A response the server
    refuses outright, such as a missing record, is not retried.
    """
    if params:
        url = f"{url}?{urlencode(params)}"
    request = Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except HTTPError as error:
            if error.code < 500 and error.code not in _RETRYABLE_STATUS:
                raise SourceUnavailable(
                    f"{description or url} was refused: HTTP {error.code} {error.reason}"
                ) from error
            if attempt == attempts:
                raise RuntimeError(f"Could not read {description or url}") from error
            time.sleep(2 ** (attempt - 1))
        except (OSError, HTTPException) as error:
            if attempt == attempts:
                raise RuntimeError(f"Could not read {description or url}") from error
            time.sleep(2 ** (attempt - 1))
    raise AssertionError("unreachable")
