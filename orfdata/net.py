"""One retrying HTTP read, shared by the source-record downloads."""

import math
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from http.client import HTTPException
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

# Defined in orfdata.cache with the rest of ORFdata's download settings; still importable from here.
from .cache import DOWNLOAD_TIMEOUT  # noqa: F401
from .version import __version__

USER_AGENT = f"orfdata/{__version__} (https://github.com/openvax/orfdata)"
MAX_ATTEMPTS = 5
# A request the server rejects will be rejected again; only a timeout or a rate
# limit is worth repeating.
_RETRYABLE_STATUS = (408, 429)


class SourceUnavailable(RuntimeError):
    """The source answered, and its answer is that it does not publish this record.

    Telling this apart from every other failure is what lets a ledger record an
    entry as ``unavailable`` rather than guessing from an error message.
    """


def parse_retry_after(value, *, now=None):
    """Parse an HTTP ``Retry-After`` header value into a nonnegative delay in seconds.

    Supports both delta-seconds (e.g. ``"120"``) and HTTP-date values (e.g.
    ``"Wed, 21 Oct 2015 07:28:00 GMT"``) per RFC 9110.

    Parameters
    ----------
    value : str or None
        Raw ``Retry-After`` header value.
    now : datetime.datetime, optional
        Reference UTC timestamp used for HTTP-date parsing; defaults to the
        current UTC time.

    Returns
    -------
    float or None
        Delay in seconds (clamped to ``>= 0.0``), or ``None`` when ``value`` is
        absent or unparseable.
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        seconds = float(text)
    except ValueError:
        seconds = None
    if seconds is not None:
        return seconds if math.isfinite(seconds) and seconds >= 0 else None
    try:
        target = parsedate_to_datetime(text)
    except (TypeError, ValueError, IndexError):
        return None
    if target is None:
        return None
    if target.tzinfo is None:
        target = target.replace(tzinfo=timezone.utc)
    reference = datetime.now(timezone.utc) if now is None else now
    if reference.tzinfo is None:
        reference = reference.replace(tzinfo=timezone.utc)
    return max(0.0, (target - reference).total_seconds())


def fetch(url, *, timeout, params=None, attempts=MAX_ATTEMPTS, description=None):
    """Read one URL, retrying transport failures and rate limits with backoff.

    A reset connection or a body that stops early is an ``OSError`` or an
    ``HTTPException``, not only a ``URLError``, and every one of them is worth
    retrying rather than losing the whole download. A response the server
    refuses outright, such as a missing record, is not retried. Retryable HTTP
    responses (408, 429, and 5xx) honour a ``Retry-After`` header when present,
    falling back to exponential backoff.

    Parameters
    ----------
    url : str
        Base HTTP(S) URL to request.
    timeout : float or None
        Per-request socket timeout in seconds, or ``None`` to wait indefinitely.
    params : dict, optional
        Query parameters appended to ``url``.
    attempts : int, optional
        Maximum total request attempts before raising ``RuntimeError``.
    description : str, optional
        Human-readable label used in error messages (defaults to ``url``).

    Returns
    -------
    bytes
        Response body bytes.
    """
    import datacache

    datacache.require_online(description or url)
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
            retry_after = parse_retry_after(
                error.headers.get("Retry-After") if error.headers else None
            )
            time.sleep(retry_after if retry_after is not None else 2 ** (attempt - 1))
        except (OSError, HTTPException) as error:
            if attempt == attempts:
                raise RuntimeError(f"Could not read {description or url}") from error
            time.sleep(2 ** (attempt - 1))
    raise AssertionError("unreachable")
