"""Write a proteome from selected host protein sources and optional viruses."""

import sys
from contextlib import redirect_stdout
from pathlib import Path

from .catalog import load_catalog
from .cli_output import print_result
from .cli_parser import UsageError
from .command_parsers import proteome_parser
from .proteome import Proteome


def _ensembl_genome(release, species):
    from pyensembl import EnsemblRelease

    return EnsemblRelease(release, species=species)


def run(argv=None, prog=None, *, genome_loader=None):
    parser = proteome_parser(prog)
    args = parser.parse_args(argv)
    try:
        # The command line is checked before any file, catalogue or genome is read.
        if args.genome is not None and (args.species is not None or args.release is not None):
            raise UsageError("Pass --genome or --species and --release, not both")
        if args.genome is None and (args.species is None or args.release is None):
            raise UsageError("Pass --species and --release, or --genome")
        names = [name.strip() for value in args.nc_orf_sources or ()
                 for name in ([value] if Path(value).is_file() else value.split(","))]
        if not all(names):
            raise UsageError("--nc-orf-sources requires nonempty source names")
        if args.no_reference_cds and not names and not args.viral:
            raise UsageError("No protein layers selected: add --nc-orf-sources or --viral")
        if args.min_protein_length is not None and not names:
            raise UsageError("--min-protein-length applies to ncORFs; add --nc-orf-sources")
        if args.min_protein_length is not None and args.min_protein_length < 1:
            raise UsageError("--min-protein-length must be a positive number of amino acids")
        sources = [load_catalog(name) if "/" in name or name.endswith((".json", ".json.gz"))
                   or Path(name).is_file() else name for name in names]
        if args.genome is not None:
            from .genomes import load_genome

            genome = load_genome(args.genome, cache_root=args.cache_root)
        else:
            genome = (genome_loader or _ensembl_genome)(args.release, args.species)
        with redirect_stdout(sys.stderr):
            proteome = Proteome.build(
                genome, reference_cds=not args.no_reference_cds, isoforms=args.isoforms,
                nc_orf_sources=sources, viral=args.viral or (), min_protein_length=args.min_protein_length,
                cache=args.cache, cache_root=args.cache_root)
        proteome.write_fasta(args.output or sys.stdout)
        annotation = proteome.write_annotation(args.annotation) if args.annotation else None
        if args.snapshot:
            proteome.save(args.snapshot)
    except (ValueError, OSError, ImportError, RuntimeError) as error:
        parser.fail(error)
    summary = proteome.summary()
    summary["excluded_orfs"] = [dict(orf_id=e.orf_id, transcript_id=e.transcript_id, reason=e.reason)
                                for e in proteome.excluded]
    if annotation is not None:
        summary["annotation"] = annotation
    if not args.json:
        summary = dict(output=args.output or "stdout (FASTA)", **summary)
    print_result(summary, json_output=args.json, stream=sys.stderr)

