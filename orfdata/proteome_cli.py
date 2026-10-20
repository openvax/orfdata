"""Write a proteome from selected host protein sources and optional viruses."""

import argparse
import json
import sys

from .catalog import load_catalog
from .cli import _add_rule_arguments, _catalog_names, _evidence_name, _rule
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
    parser.add_argument("--nc-orfs", "--nc-orf-evidence", "--extra-orf-evidence", dest="extra_orf_evidence",
                        type=_evidence_name("none"), metavar="NAME",
                        help="Add catalogue ncORFs with this evidence: 'ms', 'ribosome' (or 'riboseq'), etc.; "
                             "'all' for every catalogue ORF. ncORFs are omitted by default")
    source_group = parser.add_mutually_exclusive_group()
    source_group.add_argument("--no-reference-cds", action="store_true",
                              help="Omit proteins translated from the reference annotation's CDS features")
    source_group.add_argument("--host-sources", choices=("all", "reference-cds", "nc-orfs", "none"),
                              help=argparse.SUPPRESS)
    source_group.add_argument("--no-annotated", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--no-nc-orfs", action="store_true", help=argparse.SUPPRESS)
    _add_rule_arguments(parser)
    parser.add_argument("--catalog", action="append",
                        help="Use a local JSON/JSON.gz ORF catalogue for this genome; repeat to merge several")
    parser.add_argument("--catalogs", metavar="NAMES", type=_catalog_names,
                        help="Which of the species' installed catalogues: 'all' to add its optional ones, such as "
                             "human's Ribo-seq catalogues, or comma-separated names; its default set by default")
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
        sources = {"all": ("reference_cds", "nc_orfs"), "reference-cds": ("reference_cds",),
                   "nc-orfs": ("nc_orfs",), "none": ()}.get(args.host_sources)
        if sources is not None and args.no_nc_orfs:
            raise ValueError("--no-nc-orfs conflicts with --host-sources")
        if sources is not None and args.extra_orf_evidence is not None:
            if ("nc_orfs" in sources) == (args.extra_orf_evidence == "none"):
                raise ValueError("--nc-orfs conflicts with --host-sources")
        if args.no_nc_orfs and args.extra_orf_evidence not in (None, "none"):
            raise ValueError("--nc-orfs conflicts with --no-nc-orfs")
        requirements = (args.either or args.require_resolved or args.require_unique or args.require_reviewed
                        or any(value is not None for value in (
                            args.min_ribo_seq_studies, args.min_ms_peptides, args.ms_methods)))
        extra = ("nc_orfs" in sources if sources is not None else
                 not args.no_nc_orfs and args.extra_orf_evidence != "none" and
                 (args.extra_orf_evidence is not None or requirements))
        reference_cds = "reference_cds" in sources if sources is not None else not (
            args.no_reference_cds or args.no_annotated)
        if not extra and (args.catalog or args.catalogs or requirements):
            raise ValueError("ncORF options require --nc-orfs NAME (for example, --nc-orfs ms)")
        if not reference_cds and not extra and not args.viral:
            raise ValueError("No protein layers selected: add --nc-orfs NAME or --viral")
        rule = _rule(args, args.extra_orf_evidence, "--nc-orfs", default="ms") if extra else None
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
            genome, reference_cds=reference_cds, isoforms=args.isoforms,
            nc_orfs=(rule if rule is not None else "all") if extra else False,
            catalog=[load_catalog(path) for path in args.catalog] if args.catalog else None,
            catalogs=args.catalogs, viral=args.viral or (), cache=args.cache, cache_root=args.cache_root)
        proteome.write_fasta(args.output or sys.stdout)
        if args.snapshot:
            proteome.save(args.snapshot)
    except (ValueError, OSError, ImportError) as error:
        parser.error(str(error))
    summary = proteome.summary()
    summary["excluded_orfs"] = [dict(orf_id=e.orf_id, transcript_id=e.transcript_id, reason=e.reason)
                                for e in proteome.excluded]
    json.dump(summary, sys.stderr, indent=2)
    sys.stderr.write("\n")

