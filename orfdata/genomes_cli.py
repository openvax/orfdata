"""List and install the pinned annotations PyEnsembl's EnsemblRelease doesn't cover, such as T2T-CHM13."""

import sys
from contextlib import redirect_stdout
from dataclasses import asdict

from .cli_output import print_result
from .cli_parser import (
    ArgumentParser,
    UsageError,
    add_json,
    add_limit,
    add_timeout,
    share_cache_root,
)
from .genomes import install_genome, pinned_genomes


def run(argv=None, prog=None):
    parser = ArgumentParser(description=__doc__, prog=prog)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="Show pinned annotations, without downloading")
    listing.add_argument("--species", help="Only this species' annotations, e.g. human")
    add_json(listing)
    add_limit(listing)
    install = commands.add_parser("install", help="Download, check and index a pinned annotation")
    install.add_argument("name", help="Annotation name, as 'orfdata genomes list' shows it, e.g. t2t-chm13")
    install.add_argument("--force", action="store_true", help="Download again even when present")
    install.add_argument("--genome-fasta", help="Explicit chromosome FASTA path or URL to install through PyEnsembl")
    install.add_argument("--source-path", action="append", default=[], metavar="ROLE=PATH",
                         help="Install a local copy of the gtf, transcripts or proteins file instead of downloading")
    add_timeout(install, note="; PyEnsembl downloads --genome-fasta without one")
    add_json(install)
    share_cache_root(parser)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            output = [asdict(g) for g in pinned_genomes(species=args.species)]
        else:
            paths = {}
            for item in args.source_path:
                role, separator, path = item.partition("=")
                if not separator or not role or not path or role in paths:
                    raise UsageError("Source overrides require unique ROLE=PATH entries")
                paths[role] = path
            with redirect_stdout(sys.stderr):
                genome = install_genome(args.name, force=args.force, source_paths=paths or None,
                                        cache_root=args.cache_root, genome_fasta_path_or_url=args.genome_fasta,
                                        download_options=dict(timeout=args.timeout))
            output = dict(name=args.name, reference_name=genome.reference_name,
                          annotation=f"{genome.annotation_name} {genome.annotation_version}",
                          transcripts=len(genome.transcript_ids()))
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    print_result(output, json_output=args.json, columns=("name", "species", "reference_name", "annotation_version"),
                 limit=getattr(args, "limit", 25))
    if not args.json and args.command == "list":
        print("\nDetails: orfdata data info NAME\nInstall: orfdata genomes install NAME")
