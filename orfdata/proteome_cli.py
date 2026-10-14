"""Write a reference proteome: annotated proteins plus selected extra ORFs."""

import argparse
import json
import sys

from .catalog import load_catalog
from .evidence import EVIDENCE_METHODS, MS, RiboSeq, evidence_rule
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
    parser.add_argument("--min-ribo-seq-studies", type=int, metavar="N",
                        help="Require Ribo-seq in at least N studies the source names")
    parser.add_argument("--min-ms-peptides", type=int, metavar="N",
                        help="Require at least N distinct supporting MS peptides")
    parser.add_argument("--ms-methods", metavar="KINDS",
                        help="Count only these MS peptide kinds, comma-separated: shotgun, mhc, hla, ms_unclassified")
    parser.add_argument("--either", action="store_true",
                        help="Accept ORFs meeting the Ribo-seq or the MS requirement, instead of both")
    parser.add_argument("--require-resolved", action="store_true",
                        help="Require peptide records or named ribosome studies, not aggregate counts")
    parser.add_argument("--require-unique", action="store_true",
                        help="Require a peptide with one protein mapping in its MS source")
    parser.add_argument("--require-reviewed", action="store_true",
                        help="Require convincing reviewed spectra (benchmark scores 4-5)")
    parser.add_argument("--catalog", action="append",
                        help="Use a local JSON/JSON.gz ORF catalogue for this genome; repeat to merge several")
    parser.add_argument("--snapshot", help="Also save a checksummed snapshot that Proteome.load reads offline")
    parser.add_argument("--cache", action="store_true",
                        help="Reuse an identical earlier build from the OpenVax DataCache, or keep this one there")
    parser.add_argument("--cache-root", help="Cache directory instead of ORFDATA_CACHE, OPENVAX_DATA_CACHE or the default")
    parser.add_argument("--output", help="FASTA path; standard output by default")
    args = parser.parse_args(argv)
    try:
        if args.cache_root and not args.cache:
            raise ValueError("--cache-root needs --cache")
        evidence, options = _evidence(args)
        extra = args.extra_orf_evidence != "none"
        if not extra and (args.catalog or any(options.values())):
            raise ValueError("Extra-ORF options require --extra-orf-evidence other than 'none'")
        genome = (genome_loader or _ensembl_genome)(args.release, args.species)
        proteome = Proteome.build(
            genome, evidence_rule(evidence, **options) if extra else None, extra_orfs=extra,
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


def _evidence(args):
    """The evidence an extra ORF needs, and any require_* options for a method name."""
    options = dict(require_resolved=args.require_resolved, require_unique=args.require_unique,
                   require_reviewed=args.require_reviewed)
    ribo = args.min_ribo_seq_studies is not None
    ms = args.min_ms_peptides is not None or args.ms_methods is not None
    if args.either and not (ribo and ms):
        raise ValueError("--either needs both --min-ribo-seq-studies and an MS requirement")
    if not (ribo or ms):
        evidence = args.extra_orf_evidence or "ms"
        return (None if evidence == "all" else evidence), options
    if args.extra_orf_evidence is not None:
        raise ValueError("Use --extra-orf-evidence or the Ribo-seq/MS requirements, not both")
    if not ms and (args.require_unique or args.require_reviewed):
        raise ValueError("--require-unique and --require-reviewed need an MS requirement")
    rules = [RiboSeq(args.min_ribo_seq_studies, resolved=args.require_resolved)] if ribo else []
    if ms:
        peptides = 1 if args.min_ms_peptides is None else args.min_ms_peptides
        rules.append(MS(peptides, unique=args.require_unique, reviewed=args.require_reviewed,
                        resolved=args.require_resolved,
                        methods=args.ms_methods.split(",") if args.ms_methods is not None else None))
    if len(rules) == 1:
        return rules[0], {}
    return (rules[0] | rules[1] if args.either else rules[0] & rules[1]), {}
