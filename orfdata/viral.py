"""Offline, sequence-qualified viral CDSs bound to the existing RNA/ORF models."""

import gzip
import hashlib
import json
from dataclasses import asdict, dataclass
from functools import cached_property, lru_cache

from .reference import ReferenceContext, species_taxon
from .reference_table import packaged_bytes
from .source_cds import source_cds_orf

_CATALOG_SHA256 = "aca0c72948d4cda896d14b32203981bce4759e0c70212c05dca6f3d0e16ed4b4"


@dataclass(frozen=True)
class ViralProduct:
    """One source-labelled protein product in a viral coding reference."""

    protein_id: str
    gene: str
    product: str


@dataclass(frozen=True)
class ViralReference:
    """One exact, bundled viral reference and its source-defined coding paths.

    ``dna`` and ``orfs`` load only when accessed. Host associations describe
    documented hosts, not an exhaustive host range.
    """

    id: str
    species: str
    taxon_id: int
    host_taxon_ids: tuple[int, ...]
    host_source: str
    group: str
    virus_type: str
    ictv_species: str
    taxonomy_source: str
    length: int
    sequence_md5: str
    sequence_sha256: str
    native_to_reference_rotation: int
    source_accession: str
    source_sha256: str
    source_url: str
    coding_models: int
    products: tuple[ViralProduct, ...]
    scope: str
    catalogue_sha256: str
    bundle: str | None = None  # The installed collection this came from; None when bundled.

    @cached_property
    def dna(self):
        """Complete, checksum-verified forward FASTA bases for this reference."""
        if self.bundle is not None:
            from .viral_installed import installed_sequence
            return installed_sequence(self)
        return viral_reference_sequence(self.id)

    @cached_property
    def orfs(self):
        """Source-defined CDS paths placed on this exact reference."""
        if self.bundle is not None:
            from .viral_installed import installed_orfs
            return installed_orfs(self)
        return viral_orfs(reference_id=self.id)

    def to_dict(self):
        """Return serializable metadata without loading DNA or ORF models."""
        row = asdict(self)
        row.pop("bundle")  # A local path is not part of the reference's identity.
        row["host_taxon_ids"] = list(self.host_taxon_ids)
        row["products"] = list(row["products"])
        row["virus"] = dict(group=row.pop("group"), type=row.pop("virus_type"),
                            ictv_species=row.pop("ictv_species"),
                            taxonomy_source=row.pop("taxonomy_source"))
        row["source"] = dict(accession=row.pop("source_accession"),
                             sha256=row.pop("source_sha256"), url=row.pop("source_url"))
        return row


@lru_cache(maxsize=1)
def _catalog():
    raw = packaged_bytes("data/viral/hpv16-refseq-CDS.json.gz", _CATALOG_SHA256)
    data = json.loads(gzip.decompress(raw))
    if data["schema"] != "orfdata.viral_cds.v2":
        raise ValueError("Unsupported viral CDS catalogue schema")
    return data


def _reference(row, data):
    virus = row["virus"]
    source = row["source"]
    return ViralReference(
        id=row["id"], species=row["species"], taxon_id=row["taxon_id"],
        host_taxon_ids=tuple(row["host_taxon_ids"]), host_source=row["host_source"],
        group=virus["group"], virus_type=virus["type"],
        ictv_species=virus["ictv_species"], taxonomy_source=virus["taxonomy_source"],
        length=row["length"], sequence_md5=row["sequence_md5"],
        sequence_sha256=row["sequence_sha256"],
        native_to_reference_rotation=row["native_to_reference_rotation"],
        source_accession=source["accession"], source_sha256=source["sha256"],
        source_url=source["url"], coding_models=len(row["coding_blocks"]),
        products=tuple(ViralProduct(cds["protein_id"], cds["gene"], cds["product"])
                       for cds in data["cds"]),
        scope=data["scope"], catalogue_sha256=_CATALOG_SHA256)


def _every_reference(installed, cache_root):
    """Bundled references, then those of installed collections the bundle does not already hold.

    A collection's genome with the sequence of a bundled reference is that
    reference, so it is listed once, as the bundled one. A different sequence
    under a bundled ID would make the ID ambiguous, so that installed genome is
    ignored with a warning, as is a collection that is damaged: bundled
    lookups never fail because of what is installed.
    """
    data = _catalog()
    found = [_reference(row, data) for row in data["references"]]
    if not installed:
        return found
    import warnings

    from .viral_installed import installed_references

    held = {item.sequence_sha256 for item in found}
    taken = {item.id: item for item in found}
    for item in installed_references(cache_root):
        if item.sequence_sha256 in held:
            continue
        if item.id in taken:
            warnings.warn(f"Ignoring installed viral reference {item.id!r}: the bundled catalogue "
                          f"already uses that ID for a different sequence ({item.bundle})",
                          stacklevel=3)
            continue
        held.add(item.sequence_sha256)
        taken[item.id] = item
        found.append(item)
    return found


def viral_references(*, host=None, group=None, virus_type=None, installed=True, cache_root=None):
    """Find released viral references by documented host, group and type.

    References of installed collections are listed after the bundled ones.
    
    Parameters
    ----------
    host : str or int, optional
        Known host name such as ``"human"``, or an NCBI host taxon ID. A host
        association is source-attributed and is not a complete host range.
    group : str, optional
        Convenience group label, such as ``"HPV"``; not a taxonomy rank.
    virus_type : str, optional
        Named virus type, such as ``"HPV16"``; distinct from ICTV species.
    installed : bool, default True
        Also list the references of collections installed under ``cache_root``.
        This reads local files only. An installed reference has a ``bundle``.
        Its ``group`` is the genus in its GenBank record's lineage (a bundled
        reference's is a curated label such as ``"HPV"``),
        and its host is the isolation host that record names, with the record
        as ``host_source``; ``ictv_species`` is empty until curated (#104).
        Pass ``False`` for the bundled catalogue alone.
    cache_root : path, optional
        Where collections are installed; the default ORFdata cache otherwise.

    Returns
    -------
    tuple of ViralReference
        Exact sequence hashes, source/host pins, virus identity and coverage.
        A listed reference is not a complete transcript or exhaustive ORF
        catalogue. DNA and ORFs load only when their properties are accessed.
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
    return tuple(item for item in _every_reference(installed, cache_root)
                 if (host_taxon_id is None or host_taxon_id in item.host_taxon_ids)
                 and (group is None or item.group.casefold() == group.casefold())
                 and (virus_type is None or item.virus_type.casefold() == virus_type.casefold()))


def viral_reference(reference_id, *, installed=True, cache_root=None):
    """Return one exact viral reference by stable ID, without positional selection.

    Installed collections are searched after the bundled catalogue; see
    :func:`viral_references`.
    """
    matches = [item for item in _every_reference(installed, cache_root) if item.id == reference_id]
    if len(matches) != 1:
        raise ValueError("No unique released viral reference ID matches this name")
    return matches[0]


def viral_reference_sequence(reference_id, *, installed=True, cache_root=None):
    """Return checksum-verified reference bases for one exact ID.

    The sequence uses the forward FASTA alphabet (T rather than U) and may
    represent a virus whose biological genome is RNA. No network is accessed.
    An ID not in the bundled catalogue is looked up in installed collections.
    """
    matches = [row for row in _catalog()["references"] if row["id"] == reference_id]
    if not matches and installed:
        return viral_reference(reference_id, cache_root=cache_root).dna
    if len(matches) != 1:
        raise ValueError("No unique released viral reference ID matches this name")
    row = matches[0]
    sequence = row["reference_sequence"]
    if (len(sequence) != row["length"] or
            hashlib.sha256(sequence.encode("ascii")).hexdigest() != row["sequence_sha256"]):
        raise ValueError("Bundled viral reference sequence checksum mismatch")
    return sequence


def viral_orfs(reference_sequence=None, *, reference_id=None, contig=None, installed=True,
               cache_root=None):
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
    installed : bool, default True
        Also match the references of collections installed under ``cache_root``
        (see :func:`viral_references`). An installed collection offers a source
        CDS as an ORF only when its annotation is complete, translation starts
        at the first base and the source claims no recoding; the others stay in
        its transcriptome layer. Local files only; nothing is downloaded.
    cache_root : path, optional
        Where collections are installed; the default ORFdata cache otherwise.

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
    sequence = (viral_reference_sequence(reference_id, installed=installed, cache_root=cache_root)
                if reference_id is not None else str(reference_sequence).upper())
    fingerprint = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    data = _catalog()
    matches = [row for row in data["references"] if row["sequence_sha256"] == fingerprint
               and row["length"] == len(sequence)]
    if not matches and installed:
        from .viral_installed import installed_orfs, installed_references

        found = [item for item in installed_references(cache_root)
                 if item.sequence_sha256 == fingerprint and item.length == len(sequence)]
        if len(found) == 1:
            return installed_orfs(found[0], sequence=sequence, contig=contig)
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
        if cds["phase"] != 0:
            raise ValueError(f"Unsupported partial source codon: {protein_id}")
        orf = source_cds_orf(
            sequence, blocks, cds["strand"], context=context, label=label,
            accession=source["accession"], protein_id=protein_id, gene=cds["gene"],
            source=f"NCBI RefSeq {source['accession']}", protein=cds["protein_sequence"],
            translation_table=cds["translation_table"], start_complete=cds["start_complete"])
        if orf.stop_complete != cds["stop_complete"]:
            raise ValueError(f"Source coding boundaries do not translate exactly: {protein_id}")
        result.append(orf)
    return tuple(result)
