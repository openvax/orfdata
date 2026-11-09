"""List and install the pinned annotations PyEnsembl's EnsemblRelease doesn't cover, such as T2T-CHM13, and the names
their assemblies go by."""

import json
import sys
from contextlib import redirect_stdout
from dataclasses import asdict

from .cli_output import print_result
from .cli_parser import UsageError
from .command_parsers import genomes_parser
from .genomes import contig_aliases, install_genome, pinned_genomes, reference_names


def run(argv=None, prog=None):
    parser = genomes_parser(prog)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            output = [asdict(g) for g in pinned_genomes(species=args.species)]
        elif args.command in ("names", "contigs"):
            frame = reference_names() if args.command == "names" else contig_aliases(args.name)
            output = json.loads(frame.to_json(orient="records"))
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
    columns = {"names": ("name", "genome", "mitochondrion", "note"),
               "contigs": ("name", "ucsc", "genbank", "refseq", "length")}.get(
        args.command, ("name", "species", "reference_name", "annotation_version"))
    print_result(output, json_output=args.json, columns=columns, limit=getattr(args, "limit", 25))
    if not args.json and args.command == "list":
        print("\nDetails: orfdata data info NAME\nInstall: orfdata genomes install NAME\n"
              "Other names, such as hs1: orfdata genomes names")
