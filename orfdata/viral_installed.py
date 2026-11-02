"""Viral coding references read from installed collections.

A collection installed with ``orfdata data fetch`` (or ``download viral``) is
the genome, transcriptome and proteome layers plus their tables. The bundled
catalogue answers ``viral_references()``, ``viral_orfs()`` and
``viral_panel_coverage()`` for HPV16; this module lets the same three
functions answer for every installed collection, from the bundle's own files
and with no second format. Lookup is offline: nothing here downloads, and a
collection that is not installed contributes nothing.

A source CDS is offered as an ORF only when it binds exactly: the annotation is
complete, translation starts at the first base, the source claims no recoding,
and the bases translate to the source protein from a start and stop codon the
source supports. The rest stay in the transcriptome layer, where PyEnsembl
still sees them, and ``coding_models`` counts only the ORFs offered. Today that
leaves out frameshifts and recoding (#65) and a CDS wrapping a circular origin,
such as HBV's polymerase (#66), which SourceRNA cannot yet hold.
"""

import hashlib
import re
import warnings
from functools import lru_cache
from pathlib import Path

import pandas as pd

from .fasta import read_fasta
from .gtf import read_exons
from .model import TranscriptORF
from .reference import ReferenceContext, canonical_species
from .rna import SourceRNA

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_DECLARED = ("alternative", "declared", "non_canonical")
_UNESTABLISHED = ("partial", "unknown")


def _bundles(cache_root, only=None):
    """Each installed collection (or just ``only``) with its opened files; unreadable ones are skipped."""
    from . import viral_panel
    from .viral_registry import viral_collections

    found = []
    for collection in (viral_collections() if only is None else (only,)):
        species = canonical_species(collection.target_species)
        try:
            path = viral_panel._installed_bundle(cache_root, species, collection)
            if path is None:
                continue
            files = viral_panel._existing_bundle(path, species=species, collection=collection)
        except (ValueError, OSError) as error:
            warnings.warn(f"Ignoring the installed {collection.name} viral references: {error}",
                          stacklevel=4)
            continue
        found.append((collection, files))
    return tuple(found)


@lru_cache(maxsize=8)
def _tables(files):
    """Genome rows, transcript rows, exon paths and genome sequences of one bundle.

    The bundle's checksums are verified once per process, the first time any of
    its tables is read. A bundle is a versioned directory that is never edited, so
    this is cached for the process; a replaced bundle is a new directory.
    """
    from . import viral_panel

    viral_panel._verify_checksums(files.path)
    read = {"dtype": str, "sep": "\t", "keep_default_na": False}
    genomes = pd.read_csv(files.genome_table, **read)
    transcripts = pd.read_csv(files.transcript_table, **read)
    exons = read_exons(files.annotation_gtf)
    sequences = read_fasta(files.genome_fasta, upper=True)
    for row in genomes.itertuples():
        found = sequences.get(row.genome_id, "")
        if (len(found) != int(row.length_nt)
                or hashlib.sha256(found.encode("ascii")).hexdigest() != row.genome_sha256):
            raise ValueError(f"Viral genome {row.genome_id} does not match its table: {files.path}")
    return genomes, transcripts, exons, sequences, read_fasta(files.protein_fasta, upper=True)


def _candidates(transcripts):
    """The transcripts worth binding: complete, starting at the first base, with no claimed recoding."""
    return transcripts[(transcripts["annotation_complete"] == "yes")
                       & (transcripts["codon_start"] == "1")
                       & (transcripts["recoding"] == "")]


def _species(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _bind(files, genome_id, sequence, label):
    """Bind a genome's candidate CDSs; returns ``(orfs, skipped)``.

    A CDS that does not translate exactly to its source protein is skipped and
    named in ``skipped`` with the reason, so one odd feature does not make the
    rest of the reference unusable and nothing is dropped without a record.
    """
    genomes, transcripts, exons, _, proteins = _tables(files)
    row = genomes[genomes["genome_id"] == genome_id].iloc[0]
    context = ReferenceContext(
        taxon_id=int(row.viral_taxon_id), species=_species(row["name"]), genome_build=genome_id,
        annotation_name=files.collection, annotation_version=row.accession_version,
        sequence_sha256=row.genome_sha256)
    source = f"{files.collection} inventory {row.accession_version}"
    orfs, skipped = [], []
    for item in _candidates(transcripts[transcripts["genome_id"] == genome_id]).itertuples():
        blocks, strand = exons[item.transcript_id]
        coding = "".join(sequence[a - 1:b] if strand == "+" else
                         sequence[a - 1:b].translate(_COMPLEMENT)[::-1] for a, b in blocks)
        try:
            rna = SourceRNA(
                id=f"{row.accession_version}:CDS:{item.protein_id}", sequence=coding,
                reference=context, source=source, contig=label, strand=strand, exons=blocks,
                gene_name=item.gene, biotype="source_CDS_fragment")
            orf = TranscriptORF.at(
                rna, 0, orf_id=item.protein_id, source=source,
                protein=proteins[item.protein_id].removesuffix("*"),
                translation_table=int(item.translation_table), orf_type="reference_CDS",
                # How the source starts the CDS: a non-ATG start it declares is passed
                # on, and a start it does not establish is read as an ordinary codon.
                start_codon=item.start_codon if item.initiation in _DECLARED else None,
                start_complete=item.initiation not in _UNESTABLISHED)
            if orf.coding_sequence != coding or not orf.stop_complete:
                raise ValueError("coding boundaries do not translate exactly")
        except (ValueError, KeyError) as error:
            skipped.append((item.protein_id, str(error)))
        else:
            orfs.append(orf)
    return tuple(orfs), tuple(skipped)


@lru_cache(maxsize=512)
def _placed(files, genome_id):
    """The ORFs and skipped CDSs of one genome under its own name; cached for the process."""
    return _bind(files, genome_id, _tables(files)[3][genome_id], genome_id)


@lru_cache(maxsize=8)
def _references_of(files, collection_name):
    """One bundle's references, built once per process."""
    from .viral import ViralProduct, ViralReference

    genomes, transcripts, _, sequences, _ = _tables(files)
    catalogue = hashlib.sha256((files.path / "SHA256SUMS").read_bytes()).hexdigest()
    references = []
    for row in genomes.itertuples():
        sequence = sequences[row.genome_id]
        orfs, _ = _placed(files, row.genome_id)
        held = int((transcripts["genome_id"] == row.genome_id).sum())
        left = f"; {held - len(orfs)} of {held} source CDS not offered as ORFs" if held > len(orfs) else ""
        references.append(ViralReference(
            id=row.genome_id, species=_species(row.name), taxon_id=int(row.viral_taxon_id),
            host_taxon_ids=tuple(int(value) for value in row.source_host_taxon_ids.split(",") if value),
            host_source="", group="", virus_type=row.abbreviation, ictv_species="",
            taxonomy_source=f"{collection_name} inventory", length=len(sequence),
            sequence_md5=hashlib.md5(sequence.encode("ascii"), usedforsecurity=False).hexdigest(),
            sequence_sha256=row.genome_sha256, native_to_reference_rotation=0,
            source_accession=row.accession_version, source_sha256="", source_url="",
            coding_models=len(orfs),
            products=tuple(ViralProduct(orf.orf_id, orf.transcript.gene_name or "", "")
                           for orf in orfs),
            scope=(f"Source CDS features of the installed {collection_name} collection "
                   f"({files.path.name}); not complete transcript models{left}"),
            catalogue_sha256=catalogue, bundle=str(files.path)))
    return tuple(references)


def installed_references(cache_root=None):
    """Every reference of every installed collection, as :class:`ViralReference` rows.

    A genome whose bases an earlier one already has is listed once, so an exact
    sequence names one reference. A collection that is damaged is ignored with a
    warning rather than failing the lookup.
    """
    references, seen = [], set()
    for collection, files in _bundles(cache_root):
        try:
            found = _references_of(files, collection.name)
        except (ValueError, OSError, KeyError) as error:
            warnings.warn(f"Ignoring the installed {collection.name} viral references: {error}",
                          stacklevel=4)
            continue
        for item in found:
            if item.sequence_sha256 not in seen:
                seen.add(item.sequence_sha256)
                references.append(item)
    return tuple(references)


def _files_at(bundle):
    from . import viral_panel
    from .viral_registry import viral_collections

    for collection in viral_collections():
        species = canonical_species(collection.target_species)
        try:
            return viral_panel._existing_bundle(
                Path(bundle), species=species, collection=collection)
        except (ValueError, OSError):
            continue
    raise ValueError(f"Not an installed viral reference collection: {bundle}")


def installed_sequence(reference):
    """Checksum-verified bases of one installed reference."""
    return _tables(_files_at(reference.bundle))[3][reference.id]


def installed_orfs(reference, *, sequence=None, contig=None):
    """Source CDS paths of one installed reference, placed with SourceRNA and TranscriptORF.at."""
    label = contig if contig is not None else reference.id
    if not isinstance(label, str) or not label:
        raise ValueError("contig must be a nonempty name")
    files = _files_at(reference.bundle)
    stored = _tables(files)[3][reference.id]
    if (sequence is None or sequence == stored) and label == reference.id:
        return _placed(files, reference.id)[0]
    sequence = stored if sequence is None else sequence
    return _bind(files, reference.id, sequence, label)[0]


def installed_coverage(collection, cache_root=None):
    """For one collection's installed genomes: abbreviation -> (genome ID, ORFs offered).

    Empty when the collection is not installed.
    """
    found = {}
    for _, files in _bundles(cache_root, only=collection):
        for row in _tables(files)[0].itertuples():
            found[row.abbreviation] = (row.genome_id, len(_placed(files, row.genome_id)[0]))
    return found
