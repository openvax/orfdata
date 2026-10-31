"""The one error type ORFdata raises for data that hasn't been installed."""


class NotInstalledError(FileNotFoundError, ValueError):
    """A resource a call needs is not installed locally; the message names the command that installs it.

    It is both a ``FileNotFoundError`` and a ``ValueError``, the types these
    places raised before it existed, so existing ``except`` clauses still catch
    it. The command line exits with status 3 for it.
    """
