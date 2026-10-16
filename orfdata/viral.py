"""Offline, sequence-qualified viral CDSs bound to the existing RNA/ORF models."""

import copy
import gzip
import hashlib
import json
from functools import lru_cache
from importlib.resources import files

from .model import TranscriptORF
from .reference import ReferenceContext, species_taxon
from .rna import SourceRNA

_CATALOG_SHA256 = "aca0c72948d4cda896d14b32203981bce4759e0c70212c05dca6f3d0e16ed4b4"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@lru_cache(maxsize=1)
def _catalog():
    raw = files("orfdata").joinpath("data/viral/hpv16-refseq-CDS.json.gz").read_bytes()
    if hashlib.sha256(raw).hexdigest() != _CATALOG_SHA256:
        raise ValueError("Bundled viral CDS catalogue checksum mismatch")
    data = json.loads(gzip.decompress(raw))
    if data["schema"] != "orfdata.viral_cds.v2":
        raise ValueError("Unsupported viral CDS catalogue schema")
    return data


def viral_references(*, host=None, group=None, virus_type=None):
    """Find released viral references by documented host, group and type.

    Parameters
    ----------
    host : str or int, optional
        Known host name such as ``"human"``, or an NCBI host taxon ID. A host
        association is source-attributed and is not a complete host range.
    group : str, optional
        Convenience group label, such as ``"HPV"``; not a taxonomy rank.
    virus_type : str, optional
        Named virus type, such as ``"HPV16"``; distinct from ICTV species.

    Returns
    -------
    tuple of dict
        Exact sequence hashes, source/host pins, virus identity and coverage.
        A listed reference is not a complete transcript or exhaustive ORF
        catalogue. Returned dictionaries are independent copies.
    """
    if host is None:
        host_taxon_id = None
    elif type(host) is int and host > 0:
        host_taxon_id = host
    elif isinstance(host, str) and species_taxon(host) is not None:
        host_taxon_id = species_taxon(host)
    else:
        raise ValueError("host must be a known species name or positive NCBI taxon ID")
    if any(value is not None and (not isinstance(value, str) or not value)
           for value in (group, virus_type)):
        raise ValueError("group and virus_type must be nonempty names")
    data = _catalog()
    return tuple(dict(**{key: copy.deepcopy(value) for key, value in row.items()
                         if key not in ("coding_blocks", "reference_sequence")},
                      coding_models=len(row["coding_blocks"]), scope=data["scope"],
                      catalogue_sha256=_CATALOG_SHA256,
                      products=[dict(protein_id=cds["protein_id"], gene=cds["gene"],
                                     product=cds["product"]) for cds in data["cds"]])
                 for row in data["references"]
                 if (host_taxon_id is None or host_taxon_id in row["host_taxon_ids"])
                 and (group is None or row["virus"]["group"].casefold() == group.casefold())
                 and (virus_type is None or row["virus"]["type"].casefold() == virus_type.casefold()))


def viral_reference_sequence(reference_id):
    """Return bundled, checksum-verified reference bases for one exact ID.

    The sequence uses the forward FASTA alphabet (T rather than U) and may
    represent a virus whose biological genome is RNA. No network is accessed.
    """
    matches = [row for row in _catalog()["references"] if row["id"] == reference_id]
    if len(matches) != 1:
        raise ValueError("No unique released viral reference ID matches this name")
    row = matches[0]
    sequence = row["reference_sequence"]
    if (len(sequence) != row["length"] or
            hashlib.sha256(sequence.encode("ascii")).hexdigest() != row["sequence_sha256"]):
        raise ValueError("Bundled viral reference sequence checksum mismatch")
    return sequence


def viral_orfs(reference_sequence=None, *, reference_id=None, contig=None):
    """Load source-defined CDSs for one exactly matching viral DNA reference.

    Parameters
    ----------
    reference_sequence : str, optional
        Complete forward reference DNA from an alignment. This checks actual
        sequence identity rather than transferring annotation by contig name.
    reference_id : str, optional
        Bundled reference ID from ``viral_references()``. Use this to load the
        reference ORFs without obtaining DNA separately. Supply exactly one
        of ``reference_sequence`` or ``reference_id``.
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
    if (reference_sequence is None) == (reference_id is None):
        raise ValueError("Provide exactly one of reference_sequence or reference_id")
    sequence = (viral_reference_sequence(reference_id) if reference_id is not None
                else str(reference_sequence).upper())
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
