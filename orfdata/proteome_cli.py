"""Write a reference proteome: annotated proteins plus selected extra ORFs."""

import argparse
import json
import sys

from .catalog import EVIDENCE_METHODS, load_catalog
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
    parser.add_argument("--extra-orf-evidence", choices=EVIDENCE_METHODS + ("all", "none"), default="ms",
                        help="Evidence an extra ORF needs: 'ms' (default), 'ribosome', etc.; "
                             "'all' for every catalogue ORF; 'none' for annotated proteins only")
    parser.add_argument("--require-resolved", action="store_true",
                        help="Require peptide records or named ribosome studies, not aggregate counts")
    parser.add_argument("--require-unique", action="store_true",
                        help="Require a peptide with one protein mapping in its MS source")
    parser.add_argument("--require-reviewed", action="store_true",
                        help="Require convincing reviewed spectra (benchmark scores 4-5)")
    parser.add_argument("--catalog", help="Use a local JSON/JSON.gz ORF catalogue for this genome")
    parser.add_argument("--placements", choices=("source", "all"), default="source",
                        help="'all' places an ORF whose source transcript is gone on another transcript "
                             "with the identical spliced coding path")
    parser.add_argument("--output", help="FASTA path; standard output by default")
    args = parser.parse_args(argv)
    try:
        genome = (genome_loader or _ensembl_genome)(args.release, args.species)
        proteome = Proteome.from_genome(genome)
        if args.extra_orf_evidence != "none":
            evidence = None if args.extra_orf_evidence == "all" else args.extra_orf_evidence
            proteome = proteome.with_extra_orfs(
                evidence, catalog=load_catalog(args.catalog) if args.catalog else None,
                require_resolved=args.require_resolved, require_unique=args.require_unique,
                require_reviewed=args.require_reviewed, placements=args.placements)
        elif (args.catalog or args.require_resolved or args.require_unique or args.require_reviewed
              or args.placements != "source"):
            raise ValueError("Extra-ORF options require --extra-orf-evidence other than 'none'")
        proteome.write_fasta(args.output or sys.stdout)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    summary = proteome.summary()
    summary["excluded_orfs"] = [dict(orf_id=e.orf_id, transcript_id=e.transcript_id, reason=e.reason)
                                for e in proteome.excluded]
    json.dump(summary, sys.stderr, indent=2)
    sys.stderr.write("\n")
