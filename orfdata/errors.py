"""Shared exception types for missing installations and unaccepted resource terms."""


class NotInstalledError(FileNotFoundError, ValueError):
    """A resource a call needs is not installed locally; the message names the command that installs it.

    It is both a ``FileNotFoundError`` and a ``ValueError``, the types these
    places raised before it existed, so existing ``except`` clauses still catch
    it. The command line exits with status 3 for it.
    """


class TermsNotAccepted(ValueError):
    """A resource with usage restrictions was requested without accepting its terms.

    It is a ``ValueError`` so existing library callers still catch it, while
    :meth:`~orfdata.cli_parser.ArgumentParser.fail` maps it to command-line
    exit status 2 (missing ``--accept-terms``).
    """
