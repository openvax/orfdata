"""Common help and error conventions for the command-line entry points."""

import argparse
import sys


class ArgumentParser(argparse.ArgumentParser):
    """Show a command guide for bare groups and a useful help pointer on errors."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", argparse.RawDescriptionHelpFormatter)
        super().__init__(*args, **kwargs)
        self.set_defaults(_cli_parser=self)

    def parse_args(self, args=None, namespace=None):
        result = super().parse_args(args, namespace)
        self._selected_parser = result._cli_parser
        return result

    def add_subparsers(self, **kwargs):
        kwargs.setdefault("title", "commands")
        kwargs.setdefault("metavar", "COMMAND")
        commands = super().add_subparsers(**kwargs)
        return commands

    def parse_known_args(self, args=None, namespace=None):
        args = list(sys.argv[1:] if args is None else args)
        if not args and self._subparsers is not None:
            self.print_help()
            self.exit()
        return super().parse_known_args(args, namespace)

    def error(self, message):
        prog = getattr(self, "_selected_parser", self).prog
        self.exit(2, f"{prog}: error: {message}\nTry '{prog} --help' for usage.\n")


CACHE_ROOT_HELP = "Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default"


def add_cache_root(parser):
    parser.add_argument("--cache-root", metavar="DIR", help=CACHE_ROOT_HELP)


def share_cache_root(parser):
    """Accept --cache-root before or after any of this group's subcommands.

    A subcommand's own default would replace a value given before it, so
    below the group the option sets nothing unless it is given there.
    """
    if not any("--cache-root" in action.option_strings for action in parser._actions):
        add_cache_root(parser)
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                _defer_cache_root(subparser)


def _defer_cache_root(parser):
    found = [action for action in parser._actions if "--cache-root" in action.option_strings]
    for action in found:
        action.default, action.metavar, action.help = argparse.SUPPRESS, "DIR", CACHE_ROOT_HELP
    if not found:
        parser.add_argument("--cache-root", metavar="DIR", default=argparse.SUPPRESS, help=CACHE_ROOT_HELP)
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for subparser in action.choices.values():
                _defer_cache_root(subparser)


def seconds(value):
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("expected a positive number of seconds") from None
    if not number > 0:
        raise argparse.ArgumentTypeError("expected a positive number of seconds")
    return number


def add_timeout(parser):
    parser.add_argument("--timeout", type=seconds, default=60, metavar="SECONDS",
                        help="Seconds a download may stall before it is retried, then fails (default 60)")


def add_json(parser):
    parser.add_argument("--json", action="store_true", help="Write complete machine-readable JSON")


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
