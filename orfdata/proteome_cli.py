"""Write a proteome from selected host protein sources and optional viruses."""

import argparse
import json
import sys
from pathlib import Path

from .catalog import load_catalog
from .proteome import Proteome


def _ensembl_genome(release, species):
    try:
        from pyensembl import EnsemblRelease
    except ImportError as error:
        raise ImportError("Proteomes are built from PyEnsembl genomes: pip install 'orfdata[pyensembl]'") from error

    return EnsemblRelease(release, species=species)


def run(argv=None, prog=None, *, genome_loader=None):
    parser = argparse.ArgumentParser(
        description=__doc__, prog=prog,
        epilog="The PyEnsembl annotation must already be installed, e.g. "
               "'pyensembl install --release 115 --species human', or, for --genome, with "
               "'orfdata genomes install NAME'.")
    parser.add_argument("--species", help="PyEnsembl species, e.g. human, mouse or dog")
    parser.add_argument("--release", type=int, help="Ensembl release")
    parser.add_argument("--genome", metavar="NAME",
                        help="A pinned annotation installed with 'orfdata genomes install', e.g. t2t-chm13, "
                             "instead of --species and --release")
    parser.add_argument("--isoforms", choices=("all", "longest"), default="all",
                        help="Reference CDS isoforms: all (default), or the longest protein per gene")
    parser.add_argument("--nc-orf-sources", action="append", metavar="SOURCE[,SOURCE...]",
                        help="Add ncORFs from 'all', evidence presets such as 'ms' or 'ribosome', "
                             "dataset names, or local JSON/JSON.gz files. Repeat or separate by commas "
                             "to combine sources. No ncORFs by default")
    parser.add_argument("--no-reference-cds", action="store_true",
                        help="Omit proteins translated from the reference annotation's CDS features")
    parser.add_argument("--viral", action="append", metavar="REFERENCE",
                        help="Add a bundled virus by type (e.g. HPV16) or exact reference ID; "
                             "use 'all' for every distinct bundled virus recorded for this host; repeat for several")
    parser.add_argument("--snapshot", help="Also save a checksummed snapshot that Proteome.load reads offline")
    parser.add_argument("--cache", action="store_true",
                        help="Reuse an identical earlier build from the OpenVax DataCache, or keep this one there")
    parser.add_argument("--cache-root", help="Cache directory holding catalogues, chains, pinned genomes and cached builds, instead of "
                                             "ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", help="FASTA path; standard output by default")
    args = parser.parse_args(argv)
    try:
        sources = []
        for value in args.nc_orf_sources or ():
            names = [value] if Path(value).is_file() else value.split(",")
            for name in names:
                name = name.strip()
                if not name:
                    raise ValueError("--nc-orf-sources requires nonempty source names")
                sources.append(load_catalog(name) if "/" in name or name.endswith((".json", ".json.gz"))
                               or Path(name).is_file() else name)
        if args.no_reference_cds and not sources and not args.viral:
            raise ValueError("No protein layers selected: add --nc-orf-sources or --viral")
        if args.genome is not None:
            if args.species is not None or args.release is not None:
                raise ValueError("Pass --genome or --species and --release, not both")
            from .genomes import load_genome

            genome = load_genome(args.genome, cache_root=args.cache_root)
        elif args.species is None or args.release is None:
            raise ValueError("Pass --species and --release, or --genome")
        else:
            genome = (genome_loader or _ensembl_genome)(args.release, args.species)
        proteome = Proteome.build(
            genome, reference_cds=not args.no_reference_cds, isoforms=args.isoforms,
            nc_orf_sources=sources, viral=args.viral or (), cache=args.cache, cache_root=args.cache_root)
        proteome.write_fasta(args.output or sys.stdout)
        if args.snapshot:
            proteome.save(args.snapshot)
    except (TypeError, ValueError, OSError, ImportError) as error:
        parser.error(str(error))
    summary = proteome.summary()
    summary["excluded_orfs"] = [dict(orf_id=e.orf_id, transcript_id=e.transcript_id, reason=e.reason)
                                for e in proteome.excluded]
    json.dump(summary, sys.stderr, indent=2)
    sys.stderr.write("\n")

