"""Offline, sequence-qualified viral CDSs bound to the existing RNA/ORF models."""

import copy
import gzip
import hashlib
import json
from functools import lru_cache
from importlib.resources import files

from .model import TranscriptORF
from .reference import ReferenceContext
from .rna import SourceRNA

_CATALOG_SHA256 = "52ab8181b9d3f199dde3afabc8440ee9a8c912291c6425c99fb4225e29ad1a70"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@lru_cache(maxsize=1)
def _catalog():
    raw = files("orfdata").joinpath("data/viral/hpv16-refseq-CDS.json.gz").read_bytes()
    if hashlib.sha256(raw).hexdigest() != _CATALOG_SHA256:
        raise ValueError("Bundled viral CDS catalogue checksum mismatch")
    data = json.loads(gzip.decompress(raw))
    if data["schema"] != "orfdata.viral_cds.v1":
        raise ValueError("Unsupported viral CDS catalogue schema")
    return data


def viral_references():
    """List released viral coding-reference identities, without network access.

    Returns
    -------
    tuple of dict
        Exact sequence hashes, accession/source pins, coding-model counts and
        coverage scope. A listed reference is not a complete transcript or
        exhaustive ORF catalogue. Returned dictionaries are independent copies.
    """
    data = _catalog()
    return tuple(dict(**{key: copy.deepcopy(value) for key, value in row.items()
                         if key != "coding_blocks"},
                      coding_models=len(row["coding_blocks"]), scope=data["scope"],
                      catalogue_sha256=_CATALOG_SHA256,
                      products=[dict(protein_id=cds["protein_id"], gene=cds["gene"],
                                     product=cds["product"]) for cds in data["cds"]])
                 for row in data["references"])


def viral_orfs(reference_sequence, *, contig=None):
    """Load source-defined CDSs for one exactly matching viral DNA reference.

    Parameters
    ----------
    reference_sequence : str
        Complete forward reference DNA, not sample RNA. No annotation is
        transferred to mismatching strains or by contig name alone.
    contig : str, optional
        The caller's name for this exact DNA in an alignment reference. The
        original versioned source identity remains on every ORF's context.

    Returns
    -------
    tuple of TranscriptORF
        Source CDS fragments represented with SourceRNA and TranscriptORF.at.
        Their genomic blocks are 1-based inclusive and their RNA sequence is
        only the available CDS. Full transcript boundaries are not asserted.
        No ribosome, peptide or sample protein-expression evidence is inferred.

    Raises
    ------
    ValueError
        No released coding reference matches, or a source model fails to bind
        and translate exactly. Discovery is available through viral_references.
    """
    sequence = str(reference_sequence).upper()
    fingerprint = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    data = _catalog()
    matches = [row for row in data["references"] if row["sequence_sha256"] == fingerprint
               and row["length"] == len(sequence)]
    if len(matches) != 1:
        raise ValueError("No unique released viral coding reference matches this DNA; "
                         "inspect viral_references() for available sequence identities")
    reference = matches[0]
    source = reference["source"]
    label = contig if contig is not None else reference["id"]
    if not isinstance(label, str) or not label:
        raise ValueError("contig must be a nonempty name")
    context = ReferenceContext(
        taxon_id=reference["taxon_id"], species=reference["species"],
        genome_build=reference["id"], annotation_name="NCBI RefSeq",
        annotation_version=source["accession"], sequence_sha256=fingerprint,
        source_sha256=source["sha256"])
    result = []
    for cds in data["cds"]:
        protein_id = cds["protein_id"]
        blocks = reference["coding_blocks"][protein_id]
        coding = "".join(sequence[a - 1:b] if cds["strand"] == "+" else
                         sequence[a - 1:b].translate(_COMPLEMENT)[::-1] for a, b in blocks)
        if cds["phase"] != 0:
            raise ValueError(f"Unsupported partial source codon: {protein_id}")
        rna = SourceRNA(
            id=f"{source['accession']}:CDS:{protein_id}", sequence=coding, reference=context,
            source=f"NCBI RefSeq {source['accession']}", contig=label, strand=cds["strand"],
            exons=blocks, gene_name=cds["gene"], biotype="source_CDS_fragment")
        orf = TranscriptORF.at(rna, 0, orf_id=protein_id, source=rna.source,
                               protein=cds["protein_sequence"], translation_table=cds["translation_table"],
                               orf_type="reference_CDS", start_complete=cds["start_complete"])
        if orf.coding_sequence != coding or orf.stop_complete != cds["stop_complete"]:
            raise ValueError(f"Source coding boundaries do not translate exactly: {protein_id}")
        result.append(orf)
    return tuple(result)
