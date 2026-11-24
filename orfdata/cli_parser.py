"""Common help and error conventions for the command-line entry points.

Exit statuses: 0 on success; 1 when the command fails (network, checksums,
damaged or inconsistent data, values the data doesn't recognise); 2 when the
command line itself is invalid; 3 when a resource the command needs isn't
installed.
"""

import argparse
import os
import sys

from .cache import DOWNLOAD_TIMEOUT
from .cli_help import HelpArgumentParser
from .errors import NotInstalledError, TermsNotAccepted

FAILED, USAGE, NOT_INSTALLED = 1, 2, 3


class UsageError(ValueError):
    """Options that can't be combined, found after parsing; reported like an argparse error."""


def _apply_shared_option(parser, option_spec, *, below):
    """Attach or normalize one shared option specification on ``parser`` and its subparsers."""
    option_strings, kwargs = option_spec
    flag = option_strings[0]
    found = [action for action in parser._actions if flag in action.option_strings]
    for action in found:
        if "metavar" in kwargs:
            action.metavar = kwargs["metavar"]
        if below:
            action.default = argparse.SUPPRESS
    if not found:
        options = dict(kwargs)
        if below:
            options["default"] = argparse.SUPPRESS
        parser.add_argument(*option_strings, **options)
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                if hasattr(subparser, "_shared_options") and option_spec not in subparser._shared_options:
                    subparser._shared_options.append(option_spec)
                _apply_shared_option(subparser, option_spec, below=True)


class ArgumentParser(HelpArgumentParser):
    """Show a command guide for bare groups and a useful help pointer on errors."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._shared_options = []
        self.set_defaults(_cli_parser=self)

    def parse_args(self, args=None, namespace=None):
        result = super().parse_args(args, namespace)
        self._selected_parser = result._cli_parser
        return result

    def add_subparsers(self, **kwargs):
        kwargs.setdefault("title", "commands")
        kwargs.setdefault("metavar", "COMMAND")
        commands = super().add_subparsers(**kwargs)
        original_add_parser = commands.add_parser
        parent = self

        def add_parser(*parser_args, **parser_kwargs):
            subparser = original_add_parser(*parser_args, **parser_kwargs)
            for option_spec in parent._shared_options:
                if hasattr(subparser, "_shared_options") and option_spec not in subparser._shared_options:
                    subparser._shared_options.append(option_spec)
                _apply_shared_option(subparser, option_spec, below=True)
            return subparser

        commands.add_parser = add_parser
        return commands

    def parse_known_args(self, args=None, namespace=None):
        args = list(sys.argv[1:] if args is None else args)
        # A bare group shows its guide, except while the shell's Tab key asks argcomplete for its commands.
        if not args and self._subparsers is not None and "_ARGCOMPLETE" not in os.environ:
            self.print_help()
            self.exit()
        return super().parse_known_args(args, namespace)

    def error(self, message):
        prog = getattr(self, "_selected_parser", self).prog
        self.exit(USAGE, f"{prog}: error: {message}\nTry '{prog} --help' for usage.\n")

    def fail(self, error):
        """Report an exception from running a command, with its exit status.

        A :class:`UsageError` or :class:`~orfdata.errors.TermsNotAccepted` is
        a command-line mistake and points to help (exit 2); any other failure
        is reported without that pointer, exiting 3 for a missing installation
        and 1 otherwise.
        """
        if isinstance(error, (UsageError, TermsNotAccepted)):
            self.error(str(error))
        prog = getattr(self, "_selected_parser", self).prog
        status = NOT_INSTALLED if isinstance(error, NotInstalledError) else FAILED
        self.exit(status, f"{prog}: error: {error}\n")


CACHE_ROOT_HELP = "Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default"


def add_cache_root(parser, *, help=CACHE_ROOT_HELP, default=None):
    parser.add_argument("--cache-root", metavar="DIR", help=help, default=default)


def share_option(parser, *option_strings, below=False, **kwargs):
    """Accept an option before or after any subcommand of ``parser``, including ones added later.

    A subcommand's own default would overwrite a value parsed on its parent, so
    every descendant receives ``default=argparse.SUPPRESS`` while keeping any
    custom help text it already declared. Registering the option on
    :class:`ArgumentParser` ensures subparsers added after this call inherit
    the shared option automatically.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parent parser to share the option across.
    *option_strings : str
        Option flags such as ``"--cache-root"`` or ``"--registry"``.
    below : bool, optional
        True when applying to a child parser so its default is suppressed.
    **kwargs
        Keyword arguments forwarded to ``add_argument`` (e.g. ``metavar``,
        ``help``, ``default``).
    """
    option_spec = (tuple(option_strings), dict(kwargs))
    if hasattr(parser, "_shared_options") and option_spec not in parser._shared_options:
        parser._shared_options.append(option_spec)
    _apply_shared_option(parser, option_spec, below=below)


def share_cache_root(parser, *, below=False):
    """Accept ``--cache-root`` before or after any of this group's subcommands.

    A subcommand's own default would replace a value given before it, so
    below the group the option sets nothing unless it is given there. A
    subcommand's own help for the option is kept. Subcommands added to an
    :class:`ArgumentParser` after this call inherit ``--cache-root``
    automatically.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The group parser to share ``--cache-root`` across.
    below : bool, optional
        True when applying to a descendant parser.
    """
    share_option(parser, "--cache-root", below=below, metavar="DIR", help=CACHE_ROOT_HELP, default=None)


def seconds(value):
    """A positive, finite number of seconds, or "none" to wait indefinitely."""
    if value.lower() == "none":
        return None
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("expected a positive number of seconds, or none") from None
    if not 0 < number < float("inf"):
        raise argparse.ArgumentTypeError("expected a positive number of seconds, or none")
    return number


def add_timeout(parser, *, note=""):
    parser.add_argument("--timeout", type=seconds, default=DOWNLOAD_TIMEOUT, metavar="SECONDS",
                        help="Seconds a download may stall before it is retried, then fails "
                             f"(default %(default)s; 'none' waits indefinitely){note.replace('%', '%%')}")


def add_resume(parser):
    parser.add_argument("--resume", action="store_true",
                        help="Keep partial downloads so an interrupted one continues where it stopped; partials "
                             "and staging stay on disk after a failure (see the download guide)")


def add_json(parser, *, help="Write complete machine-readable JSON"):
    """Add ``--json``, which asks a command for machine-readable output.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser to add the option to.
    help : str, optional
        The option's help text.
    """
    parser.add_argument("--json", action="store_true", help=help)


def nonnegative(value):
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("expected a nonnegative integer") from None
    if number < 0:
        raise argparse.ArgumentTypeError("expected a nonnegative integer")
    return number


def add_limit(parser):
    parser.add_argument("--limit", type=nonnegative, default=25, metavar="N",
                        help="Human preview rows (default 25; 0 shows all); exports are always complete")


def parse_source_paths(items, *, metavar="NAME=PATH"):
    """Parse repeated ``NAME=PATH`` overrides into a dict, or ``None`` when none were given."""
    paths = {}
    for item in items:
        key, separator, path = item.partition("=")
        if not separator or not key or not path or key in paths:
            raise UsageError(f"Source overrides require unique {metavar} entries")
        paths[key] = path
    return paths or None
