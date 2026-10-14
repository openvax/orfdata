"""Write a reference proteome: annotated proteins plus selected extra ORFs."""

import argparse
import json
import sys

from .catalog import load_catalog
from .cli import _add_rule_arguments, _rule
from .evidence import EVIDENCE_METHODS
from .proteome import Proteome


def _ensembl_genome(release, species):
    from pyensembl import EnsemblRelease

    return EnsemblRelease(release, species=species)


def run(argv=None, prog=None, *, genome_loader=None):
    parser = argparse.ArgumentParser(
        description=__doc__, prog=prog,
        epilog="The PyEnsembl annotation must already be installed, e.g. "
               "'pyensembl install --release 115 --species human'.")
    parser.add_argument("--species", required=True, help="PyEnsembl species, e.g. human, mouse or dog")
    parser.add_argument("--release", required=True, type=int, help="Ensembl release")
    parser.add_argument("--extra-orf-evidence", choices=EVIDENCE_METHODS + ("all", "none"),
                        help="Evidence an extra ORF needs: 'ms' (default), 'ribosome', etc.; "
                             "'all' for every catalogue ORF; 'none' for annotated proteins only")
    _add_rule_arguments(parser)
    parser.add_argument("--catalog", action="append",
                        help="Use a local JSON/JSON.gz ORF catalogue for this genome; repeat to merge several")
    parser.add_argument("--snapshot", help="Also save a checksummed snapshot that Proteome.load reads offline")
    parser.add_argument("--cache", action="store_true",
                        help="Reuse an identical earlier build from the OpenVax DataCache, or keep this one there")
    parser.add_argument("--cache-root", help="Cache directory holding catalogues and cached builds, instead of "
                                             "ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", help="FASTA path; standard output by default")
    args = parser.parse_args(argv)
    try:
        extra = args.extra_orf_evidence != "none"
        if not extra and (args.catalog or args.either or args.ms_methods is not None or any(
                getattr(args, name) not in (None, False) for name in (
                    "min_ribo_seq_studies", "min_ms_peptides", "require_resolved", "require_unique",
                    "require_reviewed"))):
            raise ValueError("Extra-ORF options require --extra-orf-evidence other than 'none'")
        rule = _rule(args, args.extra_orf_evidence, "--extra-orf-evidence", default="ms") if extra else None
        genome = (genome_loader or _ensembl_genome)(args.release, args.species)
        proteome = Proteome.build(
            genome, rule, extra_orfs=extra,
            catalog=[load_catalog(path) for path in args.catalog] if args.catalog else None,
            cache=args.cache, cache_root=args.cache_root)
        proteome.write_fasta(args.output or sys.stdout)
        if args.snapshot:
            proteome.save(args.snapshot)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    summary = proteome.summary()
    summary["excluded_orfs"] = [dict(orf_id=e.orf_id, transcript_id=e.transcript_id, reason=e.reason)
                                for e in proteome.excluded]
    json.dump(summary, sys.stderr, indent=2)
    sys.stderr.write("\n")

