"""Cached viral source records for one analysis target, in three explicit layers.

A downloaded collection keeps the three molecular layers apart and names them
the same way everywhere:

``genome``
    Viral genomic DNA: one record per source reference, the same molecules a
    DNA reference such as GDC's ``GRCh38.d1.vd1`` carries as extra contigs.
``transcriptome``
    One coding transcript per source-annotated CDS, with the spliced
    nucleotide sequence and a GTF placing it on the viral genome.
``proteome``
    The source translation of each of those transcripts.

The GTF, transcript FASTA and protein FASTA use matching identifiers, so the
set opens as a PyEnsembl genome (:func:`viral_genome`) that varcode, isovar,
topiary and vaxrank can annotate viral transcripts with.
"""

import csv
import hashlib
import json
import os
import shutil
import stat
import tempfile
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from os import PathLike
from pathlib import Path

from .cache import DOWNLOAD_TIMEOUT, file_sha256, get_data_dir, pinned_registry
from .errors import NotInstalledError
from .fasta import write_fasta
from .record_sources import REQUEST_INTERVAL
from .reference import canonical_species, species_taxon
from .reference_ledger import empty_ledger, ledger_row, published, same_molecule, verify_record
from .viral_annotation import (
    NON_CANONICAL_START,
    gtf_gene_line,
    gtf_transcript_lines,
    name_transcripts,
    safe_token,
    unique_source_proteins,
)
from .viral_coverage import normalize_viral_contexts, viral_panel_coverage
from .viral_registry import ViralCollection, viral_collection
from .viral_subsets import normalize_viral_subsets, viral_subset_abbreviations

BUNDLE_SCHEMA_VERSION = "7"
BUNDLE_SCHEMA = "orfdata.viral_reference_bundle.v7"
INDEX_SCHEMA = "orfdata.viral_reference_index.v7"
VIRAL_LAYERS = ("genome", "transcriptome", "proteome")


def bundle_version(collection):
    """A collection's cache identity: its schema and the records it was built from."""
    return f"v{BUNDLE_SCHEMA_VERSION}-{collection.source_files[collection.primary_source].sha256}"


_GENOME_FASTA = "genomes.fasta"
_TRANSCRIPT_FASTA = "transcripts.fasta"
_PROTEIN_FASTA = "proteins.fasta"
_ANNOTATION_GTF = "annotation.gtf"
_GENOME_TABLE = "genomes.tsv"
_TRANSCRIPT_TABLE = "transcripts.tsv"
_METADATA = "metadata.json"
_CHECKSUMS = "SHA256SUMS"
# One place for source documents, whichever source served them.
_RECORD_DIRECTORY = "records"

_GENOME_COLUMNS = (
    "index", "abbreviation", "name", "genome_id", "accession_version", "source_genbank",
    "source_refseq", "viral_taxon_id", "group", "source_host", "source_host_taxon_ids",
    "length_nt", "genome_sha256",
    "transcript_count", "genome_source", "annotation_source", "alternate_source",
    "alternate_id", "alternate_sequence", "inventory_notes",
)
_TRANSCRIPT_COLUMNS = (
    "abbreviation", "genome_id", "transcript_id", "gene_id", "gene", "product", "protein_id",
    "source_protein_id", "location", "exon_count", "nucleotides", "codon_start", "start_codon",
    "initiation", "recoding", "translation_table", "amino_acids", "annotation_complete",
    "partial", "annotation_source",
)


@dataclass(frozen=True)
class SequenceFiles:
    """Files and record count for one molecular layer."""

    name: str
    molecule: str
    description: str
    fasta: Path
    records: int
    table: Path | None = None
    annotation: Path | None = None

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-ready description of this layer."""
        return {key: str(value) if isinstance(value, Path) else value
                for key, value in asdict(self).items()}


@dataclass(frozen=True)
class ViralReferenceFiles:
    """The genome, transcript and protein files opened for a viral collection."""

    target_species: str
    target_taxon_id: int
    collection: str
    path: Path
    genome_fasta: Path
    transcript_fasta: Path
    protein_fasta: Path
    annotation_gtf: Path
    genome_table: Path
    transcript_table: Path
    genome_count: int
    transcript_count: int
    protein_count: int

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-ready description of these files."""
        result = {key: str(value) if isinstance(value, Path) else value
                  for key, value in asdict(self).items()}
        result["layers"] = [layer.to_dict() for layer in self.layers]
        return result

    @property
    def layers(self) -> tuple[SequenceFiles, ...]:
        """The genome, transcriptome and proteome layers, in that order."""
        return (
            SequenceFiles(
                name="genome", molecule="dna",
                description="Viral genomic DNA, one record per source reference",
                fasta=self.genome_fasta, records=self.genome_count, table=self.genome_table),
            SequenceFiles(
                name="transcriptome", molecule="rna",
                description="One coding transcript per source-annotated CDS, spliced",
                fasta=self.transcript_fasta, records=self.transcript_count,
                table=self.transcript_table, annotation=self.annotation_gtf),
            SequenceFiles(
                name="proteome", molecule="protein",
                description="Source translation of each coding transcript",
                fasta=self.protein_fasta, records=self.protein_count,
                table=self.transcript_table),
        )

    def layer(self, name: str) -> SequenceFiles:
        """One named layer: ``genome``, ``transcriptome`` or ``proteome``."""
        found = {layer.name: layer for layer in self.layers}
        if name not in found:
            raise ValueError(f"Unknown viral layer {name!r}; available: {', '.join(VIRAL_LAYERS)}")
        return found[name]

    def pyensembl_genome(self, *, cache_root=None, rebuild=False):
        """Open these files as a PyEnsembl genome; see :func:`viral_genome`."""
        return viral_genome(self, cache_root=cache_root, rebuild=rebuild)


@dataclass(frozen=True)
class DownloadedViralReferences(ViralReferenceFiles):
    """A complete, downloaded viral source-reference collection for a target species."""

    contexts: tuple[str, ...]
    panel: str
    record_directory: Path
    retrieved: str


@dataclass(frozen=True)
class SelectedViralReferences(ViralReferenceFiles):
    """Files selected by target-species context or named subset."""

    contexts: tuple[str, ...]
    excluded_contexts: tuple[str, ...]
    context_match: str
    subsets: tuple[str, ...]
    source_collection: Path


def _selection(species, collection):
    """Resolve a species and collection name to the published collection it names."""
    chosen = viral_collection(collection, species=species)
    return canonical_species(chosen.target_species), chosen


def _primary_ledger(collection, inventory):
    """What the primary source publishes for each inventory row, in inventory order.

    The inventory names an entry; the ledger says which record at the source
    that name resolved to, and pins it. Reading them together here is what
    keeps a collection from downloading a record its own inventory never named.
    """
    source = collection.source()
    ledger = collection.source_summary(collection.primary_source)
    for entry in inventory.itertuples(index=False):
        row = ledger_row(ledger, entry.abbreviation)
        named = source.identifier(entry)
        if not published(row):
            raise ValueError(
                f"The {source.name} ledger does not publish {entry.abbreviation}")
        if row["index"] != entry.index or row["identifier"].split(".")[0] != named.split(".")[0]:
            raise ValueError(
                f"The {collection.panel} inventory names {named} for {entry.abbreviation}, "
                f"but its {source.name} ledger resolves {row['identifier']}")
    return ledger.loc[inventory["abbreviation"]]


def _write_checksums(path):
    with (path / _CHECKSUMS).open("w") as stream:
        for file_path in sorted((item for item in path.rglob("*") if item.is_file()
                                 and item.name != _CHECKSUMS
                                 and not item.name.startswith(".datacache-")),
                                key=lambda item: item.relative_to(path).as_posix()):
            stream.write(f"{file_sha256(file_path)}  {file_path.relative_to(path).as_posix()}\n")


def _quote_empty_final_tsv_fields(path):
    """Keep empty trailing TSV values parseable without trailing whitespace."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    path.write_text("".join(
        f'{line[:-2]}\t""\n' if line.endswith("\t\n") else line
        for line in lines), encoding="utf-8")


def _write_table(path, columns, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    _quote_empty_final_tsv_fields(path)


def _read_table(path):
    with path.open(newline="") as stream:
        return tuple(csv.DictReader(stream, delimiter="\t"))


def _count_fasta(path):
    with path.open() as stream:
        return sum(line.startswith(">") for line in stream)


def _metadata(path, *, expected_schema):
    try:
        metadata = json.loads((path / _METADATA).read_text())
        if not isinstance(metadata, dict):
            raise ValueError("metadata is not an object")
        if metadata.get("schema") != expected_schema:
            raise ValueError(f"unexpected schema {metadata.get('schema')!r}")
        for key in ("genome_count", "transcript_count", "protein_count"):
            metadata[key] = int(metadata[key])
    except PermissionError:
        raise  # Unreadable, not invalid: replacing it would hide the problem.
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid viral reference metadata in {path}: {error}") from error
    return metadata


def _paths(path):
    return dict(
        path=path,
        genome_fasta=path / _GENOME_FASTA,
        transcript_fasta=path / _TRANSCRIPT_FASTA,
        protein_fasta=path / _PROTEIN_FASTA,
        annotation_gtf=path / _ANNOTATION_GTF,
        genome_table=path / _GENOME_TABLE,
        transcript_table=path / _TRANSCRIPT_TABLE,
    )


def _bundle_result(path, *, species, collection):
    metadata = _metadata(path, expected_schema=BUNDLE_SCHEMA)
    if metadata.get("bundle_version") != bundle_version(collection):
        raise ValueError(f"Invalid viral reference metadata in {path}: unexpected version")
    return DownloadedViralReferences(
        target_species=species,
        target_taxon_id=collection.target_taxon_id,
        collection=collection.name,
        genome_count=metadata["genome_count"],
        transcript_count=metadata["transcript_count"],
        protein_count=metadata["protein_count"],
        contexts=collection.contexts,
        panel=collection.panel,
        record_directory=path / _RECORD_DIRECTORY,
        retrieved=str(metadata.get("retrieved", "")),
        **_paths(path),
    )


def _existing_bundle(path, *, species, collection):
    """Check the bundle structure without hashing every installed file."""
    result = _bundle_result(path, species=species, collection=collection)
    required = (result.genome_fasta, result.transcript_fasta, result.protein_fasta,
                result.annotation_gtf, result.genome_table, result.transcript_table,
                path / _CHECKSUMS)
    if any(not file_path.is_file() for file_path in required):
        raise ValueError(f"Incomplete viral reference collection: {path}")
    if not result.record_directory.is_dir():
        raise ValueError(f"Incomplete viral source records: {result.record_directory}")
    # os.listdir, unlike glob, raises for a directory it can't read rather than finding nothing in it.
    if sum(name.endswith(".gz") for name in os.listdir(result.record_directory)) != result.genome_count:
        raise ValueError(f"Incomplete viral source records: {result.record_directory}")
    return result


def _installed_files(path):
    """Every file under ``path``; an unreadable directory raises, where rglob would skip it."""
    def refuse(error):
        raise error

    return {Path(root, name).relative_to(path).as_posix()
            for root, _, names in os.walk(path, onerror=refuse) for name in names}


def _verify_checksums(path):
    sums_path = path / _CHECKSUMS
    try:
        sums = tuple(line.rstrip("\n").split("  ", 1) for line in sums_path.open())
    except PermissionError:
        raise
    except OSError as error:
        raise ValueError(f"Missing viral reference checksums: {sums_path}") from error
    expected_files = {
        relative for relative in _installed_files(path)
        if relative != _CHECKSUMS and not Path(relative).name.startswith(".datacache-")
    }
    if (len(sums) != len(expected_files) or any(len(parts) != 2 for parts in sums)
            or {relative for _, relative in sums} != expected_files):
        raise ValueError(f"Invalid viral reference checksum manifest: {sums_path}")
    if any(file_sha256(path / relative) != digest for digest, relative in sums):
        raise ValueError(f"Viral reference checksum mismatch: {path}")


def _verify_counts(result):
    """Every layer must hold exactly the records its metadata claims."""
    for fasta, expected in ((result.genome_fasta, result.genome_count),
                            (result.transcript_fasta, result.transcript_count),
                            (result.protein_fasta, result.protein_count)):
        if _count_fasta(fasta) != expected:
            raise ValueError(f"FASTA record count does not match metadata: {fasta}")
    for table, expected in ((result.genome_table, result.genome_count),
                            (result.transcript_table, result.transcript_count)):
        if len(_read_table(table)) != expected:
            raise ValueError(f"Table row count does not match metadata: {table}")


def _verified_existing_bundle(path, *, species, collection):
    result = _existing_bundle(path, species=species, collection=collection)
    _verify_checksums(path)
    _verify_counts(result)
    return result


def _genome_row(entry, primary, molecule, transcripts, alternate, record):
    """One row of the genome layer's table, as the inventory and the record describe it."""
    sources = {transcript["annotation_source"] for transcript in transcripts}
    return {
        "index": entry.index,
        "abbreviation": entry.abbreviation,
        "name": entry.name,
        "genome_id": molecule.genome_id,
        "accession_version": primary["identifier"],
        "source_genbank": entry.genbank,
        "source_refseq": entry.refseq,
        "viral_taxon_id": primary["taxon_id"],
        # What the primary GenBank record says: its genus and the host it was isolated from. The target species is an analysis
        # choice and says nothing about either.
        "group": record.group,
        "source_host": record.host,
        "source_host_taxon_ids": str(species_taxon(record.host.split(";")[0]) or ""),
        "length_nt": len(molecule.sequence),
        "genome_sha256": hashlib.sha256(molecule.sequence.encode("ascii")).hexdigest(),
        "transcript_count": len(transcripts),
        # Which source this molecule came from, and which annotated it.
        "genome_source": molecule.source,
        "annotation_source": ",".join(sorted(sources)),
        # The other reference the inventory names for this row, and whether this
        # molecule is it; where they differ, the reference the collection follows
        # was built from the alternate.
        "alternate_source": alternate["source"] if alternate is not None else "",
        "alternate_id": alternate["identifier"] if alternate is not None else "",
        "alternate_sequence": molecule.alternate_state,
        "inventory_notes": entry.pave_notes or "",
    }


def _transcript_row(abbreviation, transcript):
    """One row of the transcriptome layer's table, which the proteome layer shares."""
    return {
        "abbreviation": abbreviation,
        "genome_id": transcript["genome_id"],
        "transcript_id": transcript["transcript_id"],
        "gene_id": transcript["gene_id"],
        "gene": transcript["gene"],
        "product": transcript["product"],
        "protein_id": transcript["protein_id"],
        "source_protein_id": transcript["source_protein_id"],
        "location": transcript["location"],
        "exon_count": len(transcript["exons"]),
        "nucleotides": len(transcript["nucleotides"]),
        "codon_start": transcript["codon_start"],
        # The transcriptome records how translation begins; the residue it
        # yields belongs to the proteome, which publishes it as the protein.
        "start_codon": transcript["start_codon"],
        "initiation": transcript["initiation"],
        "recoding": transcript["recoding"],
        "translation_table": transcript["translation_table"],
        "amino_acids": len(transcript["translation"]),
        "annotation_complete": "yes" if transcript["complete"] else "no",
        # Which end, if either, the source marks as running past the record.
        "partial": transcript["partial"],
        "annotation_source": transcript["annotation_source"],
    }


@dataclass(frozen=True)
class _Molecule:
    """The sequence a reference row is, the source it came from, and that source's record."""

    genome_id: str
    sequence: str
    source: str
    features: tuple
    alternate_state: str
    record: object = None

    @property
    def document(self):
        """The source record to store beside this molecule."""
        return self.record


def _resolve_molecule(entry, record, ledger, *, collection, prefer_reference_molecule, timeout,
                      notes):
    """Choose the molecule for one row, and say which source it came from.

    An inventory can name an entry at more than one source, and those are not
    always the same sequence. Where they differ, the molecule to carry is the
    one the reference the collection follows was built from: for GDC's viral
    decoys that is the alternate source's, which is checkable for HPV16 — its
    PaVE reference reproduces GDC's published DNA checksum exactly, while the
    GenBank accession beside it is two bases shorter.

    A collection with no alternate source, or a row the alternate does not
    publish, keeps the primary record untouched.
    """
    alternate = ledger_row(ledger, entry.abbreviation)
    if not published(alternate):
        state = "unavailable" if alternate is not None else ""
        return _Molecule(genome_id=record.identifier, sequence=record.sequence,
                         source=record.source, features=record.coding_features,
                         alternate_state=state, record=record)
    if same_molecule(alternate, record.sequence):
        return _Molecule(genome_id=record.identifier, sequence=record.sequence,
                         source=record.source, features=record.coding_features,
                         alternate_state="same", record=record)
    if not prefer_reference_molecule:
        return _Molecule(genome_id=record.identifier, sequence=record.sequence,
                         source=record.source, features=record.coding_features,
                         alternate_state="differs", record=record)
    try:
        fetched = verify_record(alternate, collection.source(alternate["source"]).record(
            alternate["identifier"], timeout=timeout))
    except (RuntimeError, ValueError, OSError) as error:
        # Keep the row on its primary record rather than lose it, and say so.
        notes[entry.abbreviation] = str(error)
        return _Molecule(genome_id=record.identifier, sequence=record.sequence,
                         source=record.source, features=record.coding_features,
                         alternate_state="differs", record=record)
    return _Molecule(genome_id=fetched.identifier, sequence=fetched.sequence,
                     source=fetched.source, features=fetched.coding_features,
                     alternate_state="same", record=fetched)


def _annotate_from_alternate(molecule, entry, ledger, *, collection, timeout, notes,
                             prefer_alternate=False):
    """Take an alternate source's coding features for a molecule it also publishes.

    A few inventory rows point at submissions that annotate nothing at all.
    Where the inventory also names an alternate reference that is the identical
    molecule, its coding features apply without any coordinate transfer. This
    adds to the primary record rather than standing in for it, so an alternate
    that is unreachable leaves the row without transcripts, with its reason
    recorded, instead of discarding the whole download.

    By default only an unannotated molecule is filled, because a submission
    that annotated its own genes is the record the inventory names. With
    ``prefer_alternate`` the alternate's annotation is taken wherever it
    publishes the identical molecule and annotates at least as much of it,
    which trades the submissions' own gene naming for one source's uniform
    models across the collection. Either way the sequence is untouched: this changes only how the
    same bases are read.
    """
    alternate = ledger_row(ledger, entry.abbreviation)
    if not published(alternate):
        return molecule
    # Asking for the alternate's models means asking for them wherever it has
    # them: a row where both annotate the same number of features is exactly
    # the row where the two namings and splicings differ, which is the point.
    # Only an alternate that annotates strictly less is declined.
    if molecule.features and not (prefer_alternate
                                  and alternate["coding_features"] >= len(molecule.features)
                                  and same_molecule(alternate, molecule.sequence)):
        return molecule
    try:
        fetched = collection.source(alternate["source"]).record(
            alternate["identifier"], timeout=timeout)
        if fetched.sequence != molecule.sequence:
            raise ValueError(
                f"{alternate['source']} {alternate['identifier']} is not the same sequence as the "
                f"source record, so its coordinates do not apply"
            )
    except (RuntimeError, ValueError, OSError) as error:
        notes[entry.abbreviation] = str(error)
        return molecule
    return replace(molecule, features=fetched.coding_features)


def _write_bundle(path, collection, inventory, primaries, *, timeout,
                  prefer_reference_molecule=True, annotate_from_alternate=True,
                  prefer_alternate_annotation=False):
    """Fetch every reference once and write the three layers, with their provenance."""
    records_directory = path / _RECORD_DIRECTORY
    records_directory.mkdir(exist_ok=True)
    source = collection.source()
    ledger = (collection.source_summary(collection.alternate_source)
              if collection.alternate_source else empty_ledger())
    entries = list(inventory.itertuples(index=False))
    primaries = [row for _, row in primaries.iterrows()]
    pinned = {primary["identifier"]: primary for primary in primaries}
    downloaded, batch_size = {}, getattr(source, "batch_size", 1)
    request_interval = getattr(source, "request_interval", REQUEST_INTERVAL)
    for start_index in range(0, len(primaries), batch_size):
        batch = primaries[start_index:start_index + batch_size]
        for record in source.records([item["identifier"] for item in batch], timeout=timeout):
            if record.identifier not in pinned:
                raise ValueError(f"{source.name} returned an unexpected record {record.identifier}")
            if record.identifier in downloaded:
                raise ValueError(f"{source.name} returned {record.identifier} twice")
            # Every source is checked against its own ledger, primary or not.
            downloaded[record.identifier] = verify_record(pinned[record.identifier], record)
        if start_index + batch_size < len(primaries):
            time.sleep(request_interval)  # Respect the source's request rate.
    if set(downloaded) != set(pinned):
        missing = sorted(set(pinned) - set(downloaded))
        extra = sorted(set(downloaded) - set(pinned))
        raise ValueError(f"{source.name} response mismatch; missing={missing}, extra={extra}")

    notes = {}
    molecules = {}
    for entry, primary in zip(entries, primaries, strict=True):
        molecule = _resolve_molecule(
            entry, downloaded[primary["identifier"]], ledger, collection=collection,
            prefer_reference_molecule=prefer_reference_molecule, timeout=timeout, notes=notes)
        if annotate_from_alternate:
            molecule = _annotate_from_alternate(
                molecule, entry, ledger, collection=collection, timeout=timeout,
                notes=notes, prefer_alternate=prefer_alternate_annotation)
        molecules[entry.index] = molecule

    # A source protein accession names the protein only when it names one protein.
    unique_protein_ids = unique_source_proteins(
        feature for molecule in molecules.values() for feature in molecule.features)

    genome_rows, transcript_rows, gtf_lines, complete = [], [], [], 0
    with (path / _GENOME_FASTA).open("w") as genome_stream, \
         (path / _TRANSCRIPT_FASTA).open("w") as transcript_stream, \
         (path / _PROTEIN_FASTA).open("w") as protein_stream:
        for entry, primary in zip(entries, primaries, strict=True):
            molecule = molecules[entry.index]
            genome_id = molecule.genome_id
            molecule.record.write_document(records_directory, entry.abbreviation)
            write_fasta(genome_stream, genome_id, molecule.sequence,
                        description=f"abbreviation={safe_token(entry.abbreviation)} "
                                    f"taxon={primary['taxon_id']} molecule=genome "
                                    f"source={molecule.source}")
            placed = name_transcripts(genome_id, molecule.features,
                                      unique_protein_ids=unique_protein_ids)
            for transcript in placed:
                shared = (f"genome={genome_id} gene={transcript['gene_name']} "
                          f"product={safe_token(transcript['product'])}")
                write_fasta(transcript_stream, transcript["transcript_id"],
                            transcript["nucleotides"],
                            description=f"{shared} protein_id={transcript['protein_id']}")
                write_fasta(protein_stream, transcript["protein_id"], transcript["translation"],
                            description=f"{shared} transcript={transcript['transcript_id']}")
                transcript_rows.append(_transcript_row(entry.abbreviation, transcript))
            complete += sum(row["complete"] for row in placed)
            for gene_id in dict.fromkeys(row["gene_id"] for row in placed):
                members = [row for row in placed if row["gene_id"] == gene_id]
                gtf_lines.append(gtf_gene_line(genome_id, gene_id, members))
                for transcript in members:
                    gtf_lines.extend(gtf_transcript_lines(transcript))
            genome_rows.append(_genome_row(entry, primary, molecule, placed,
                                           ledger_row(ledger, entry.abbreviation),
                                           downloaded[primary["identifier"]]))

    _write_table(path / _GENOME_TABLE, _GENOME_COLUMNS, genome_rows)
    _write_table(path / _TRANSCRIPT_TABLE, _TRANSCRIPT_COLUMNS, transcript_rows)
    (path / _ANNOTATION_GTF).write_text(
        f"#!genome-build {collection.panel} viral references ({collection.name})\n"
        f"#!annotation-source orfdata {bundle_version(collection)}; "
        f"{', '.join(sorted({row['annotation_source'] for row in transcript_rows}))}\n"
        + "".join(gtf_lines), encoding="utf-8")
    retrieved = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    (path / _METADATA).write_text(json.dumps({
        "schema": BUNDLE_SCHEMA,
        "bundle_version": bundle_version(collection),
        "collection": collection.name,
        "panel": collection.panel,
        "target_species": canonical_species(collection.target_species),
        "target_taxon_id": collection.target_taxon_id,
        "contexts": list(collection.contexts),
        "retrieved": retrieved,
        "sources": sorted({row["genome_source"] for row in genome_rows}),
        "genome_count": len(genome_rows),
        "transcript_count": len(transcript_rows),
        "protein_count": len(transcript_rows),
        "complete_transcript_count": complete,
        "genome_sources": dict(sorted(Counter(row["genome_source"]
                                              for row in genome_rows).items())),
        "annotation_sources": dict(sorted(Counter(row["annotation_source"]
                                                  for row in transcript_rows).items())),
        "source_notes": notes,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (path / "README.md").write_text(
        _bundle_readme(collection, genome_rows, transcript_rows, complete, retrieved),
        encoding="utf-8")
    _write_checksums(path)


def _bundle_readme(collection, genome_rows, transcript_rows, complete, retrieved):
    sources = dict(sorted(Counter(row["genome_source"] for row in genome_rows).items()))
    alternate = [row["abbreviation"] for row in genome_rows
                 if row["alternate_sequence"] == "differs"]
    return (
        f"# {collection.panel} viral reference collection ({collection.name})\n\n"
        f"{collection.title}, retrieved at {retrieved} for analysis of "
        f"{collection.target_species} (NCBI taxon {collection.target_taxon_id}).\n\n"
        "Three layers, kept separate:\n\n"
        f"- genome: `{_GENOME_FASTA}`, {len(genome_rows)} viral genomic DNA records, described by "
        f"`{_GENOME_TABLE}`. These are the molecules a DNA reference such as {collection.panel} "
        "carries as extra contigs.\n"
        f"- transcriptome: `{_TRANSCRIPT_FASTA}`, {len(transcript_rows)} spliced coding "
        f"transcripts, one per source CDS feature, placed on the genomes by `{_ANNOTATION_GTF}` "
        f"and described by `{_TRANSCRIPT_TABLE}`.\n"
        f"- proteome: `{_PROTEIN_FASTA}`, the source translation of each transcript.\n\n"
        f"{complete} of {len(transcript_rows)} transcripts have whole start and stop codons "
        "(`annotation_complete`); the rest are frameshifted, or marked by the source as running "
        "past the record, which the `partial` column names as `5_prime`, `3_prime` or `both`.\n\n"
        f"Molecules came from {', '.join(f'{name} ({count})' for name, count in sources.items())}. "
        f"`{_RECORD_DIRECTORY}/` holds each source's own document, exactly as it was served, so "
        f"`{_CHECKSUMS}` attests to the source record rather than to a re-serialization of it. "
        "`genome_source` and `annotation_source` say where every genome and transcript came "
        "from.\n\n"
        + (f"{len(alternate)} genomes are not the alternate reference the inventory names: "
           f"{', '.join(alternate)}.\n\n" if alternate else "")
        + "The GTF, transcript FASTA and protein FASTA share identifiers, so this set opens as a "
        "PyEnsembl genome for varcode, isovar, topiary and vaxrank.\n\n"
        f"The `{collection.name}` collection is for {collection.target_species} analysis. That "
        "target does not assert a source host. `source_host` is the host each primary GenBank "
        "record names for its sequence, an isolation host and not a host range, and "
        "`source_host_taxon_ids` its NCBI taxon where ORFdata knows the species; both are empty "
        "when the record says nothing. `group` is the genus in the record's taxonomic "
        "lineage.\n\n"
        "This is a source-record download, not an ORFdata curated viral-reference release. "
        "It does not alter `viral_references()` or `Proteome.build(..., viral=...)`.\n"
    )


# Collections and selections are published through DataCache's public registries
# and materialization stores (VersionedDatasetRegistry for published release
# artifacts; datacache.materialize for local builds and derived selections).
# Each version has a DataCache store holding immutable bundles under bundles/,
# the newest of which is current. Replaced bundles remain for running readers
# until prune_viral_references() or 'orfdata data prune' removes them. ORFdata
# 0.43 and earlier kept one plain directory per version under "references" and
# "indexes"; such a directory is still read while no DataCache bundle exists.
_COLLECTIONS, _LEGACY_COLLECTIONS = "collections", "references"
_DOWNLOADS = "downloads"
_SELECTIONS, _LEGACY_SELECTIONS = "selections", "indexes"
_BASE_OUTPUTS = (
    _GENOME_FASTA,
    _TRANSCRIPT_FASTA,
    _PROTEIN_FASTA,
    _ANNOTATION_GTF,
    _GENOME_TABLE,
    _TRANSCRIPT_TABLE,
    _METADATA,
    "README.md",
    _CHECKSUMS,
)


@dataclass(frozen=True)
class ViralInspection:
    """Read-only inspection result for a viral collection.

    Attributes
    ----------
    status : str
        Normalized local status: ``"available"``, ``"missing"``,
        ``"damaged"``, or ``"inaccessible"``.
    path : pathlib.Path or None
        Resolved bundle directory when available, or store/bundle path when
        damaged or inaccessible; ``None`` when missing.
    origin : str or None
        ``"download"`` for a published DataCache artifact, ``"build"`` for a
        local DataCache materialization (or 0.43 directory), or ``None`` when
        not available.
    datacache_status : str or None
        Raw underlying DataCache inspection status (``"available"``,
        ``"missing"``, ``"invalid"``, or ``"inaccessible"``), or ``None`` for a
        legacy 0.43 directory.
    error : Exception or None
        Underlying error when ``status`` is ``"damaged"`` or
        ``"inaccessible"``.
    """

    status: str
    path: Path | None
    origin: str | None = None
    datacache_status: str | None = None
    error: Exception | None = None


def viral_store_path(cache_root, species, collection, *, area=_COLLECTIONS):
    """Return the DataCache store directory for a viral collection version.

    Parameters
    ----------
    cache_root : str or os.PathLike, optional
        Explicit cache root override.
    species : str
        Canonical target species name.
    collection : ViralCollection or str
        Published collection object or name.
    area : str, optional
        Cache area under ``viral/`` (``"collections"`` for local builds,
        ``"downloads"`` for published artifacts, or ``"references"`` for legacy
        0.43 directories).

    Returns
    -------
    pathlib.Path
        Versioned store path.
    """
    collection = viral_collection(collection)
    return (get_data_dir(cache_root) / "viral" / area / species /
            collection.name / bundle_version(collection))


_bundle_path = viral_store_path


def viral_download_registry(collection, cache_root=None):
    """Return the DataCache :class:`~datacache.VersionedDatasetRegistry` for a published viral collection.

    Parameters
    ----------
    collection : ViralCollection or str
        Published collection declaring checksummed ``assets``.
    cache_root : str or os.PathLike, optional
        Explicit cache root override.

    Returns
    -------
    datacache.VersionedDatasetRegistry
        Registry installing the collection's published artifact into its
        ``viral/downloads/`` store.

    Raises
    ------
    ValueError
        If ``collection`` has no published artifact registered.
    """
    collection = viral_collection(collection)
    if not collection.published:
        raise ValueError(
            f"The {collection.name} viral collection has no published artifact; "
            "it is built locally from pinned source ledgers"
        )
    species = canonical_species(collection.target_species)
    assets = {name: dict(spec) for name, spec in collection.assets.items()}
    return pinned_registry(
        collection.name,
        bundle_version(collection),
        assets,
        viral_store_path(cache_root, species, collection, area=_DOWNLOADS),
    )


def _is_directory(path):
    """Whether ``path`` is a real directory, not a link to one; an unreadable one counts as absent."""
    try:
        return stat.S_ISDIR(Path(path).lstat().st_mode)
    except OSError:
        return False


def _legacy_directory(path):
    """Whether a 0.43 directory is there, following a link to it as 0.43 did; unreadable counts as absent."""
    try:
        return Path(path).is_dir()
    except OSError:
        return False


def _collection_repair(collection):
    name = viral_collection(collection).name
    return f"download_viral_references(force=True), or 'orfdata data fetch {name} --force', repairs it"


def _validate_bundle_path(bundle_dir, *, species, collection, origin, verify_files):
    """Validate an installed bundle directory and return a :class:`ViralInspection`."""
    check = _verified_existing_bundle if verify_files else _existing_bundle
    try:
        check(bundle_dir, species=species, collection=collection)
    except PermissionError as error:
        return ViralInspection(
            status="inaccessible", path=bundle_dir, origin=None,
            datacache_status="inaccessible", error=error,
        )
    except OSError as error:
        return ViralInspection(
            status="inaccessible", path=bundle_dir, origin=None,
            datacache_status="inaccessible", error=error,
        )
    except ValueError as error:
        return ViralInspection(
            status="damaged", path=bundle_dir, origin=None,
            datacache_status="invalid", error=error,
        )
    return ViralInspection(
        status="available", path=bundle_dir, origin=origin,
        datacache_status="available", error=None,
    )


def inspect_viral_collection(
    collection: ViralCollection | str | None = None,
    *,
    species: str | int | None = None,
    cache_root: str | PathLike[str] | None = None,
    verify_files: bool = True,
) -> ViralInspection:
    """Inspect a viral collection's local DataCache stores offline without writing or downloading.

    Checks the published artifact store (when the collection is published),
    then the local build materialization store, and finally any legacy ORFdata
    0.43 directory.

    Parameters
    ----------
    collection : ViralCollection or str, optional
        Collection object, name, alias, or panel; defaults to the sole
        collection for ``species``.
    species : str or int, optional
        Target species alias or NCBI taxon ID.
    cache_root : str or os.PathLike, optional
        Explicit cache root override.
    verify_files : bool, optional
        When ``True``, verify file hashes and bundle counts against the
        DataCache receipt and ``SHA256SUMS``.

    Returns
    -------
    ViralInspection
        Normalized status (``"available"``, ``"missing"``, ``"damaged"``, or
        ``"inaccessible"``), resolved path, origin, and raw DataCache status.
    """
    import datacache

    species, collection = _selection(species, collection)
    first_problem = None

    if collection.published:
        dl_store = viral_store_path(cache_root, species, collection, area=_DOWNLOADS)
        dl_inspected = viral_download_registry(collection, cache_root).inspect(
            collection.name, bundle_version(collection), verify_files=verify_files,
        )
        if dl_inspected.status == "available":
            bundle_dir = Path(dl_inspected.bundle)
            checked = _validate_bundle_path(
                bundle_dir, species=species, collection=collection,
                origin="download", verify_files=verify_files,
            )
            if checked.status == "available":
                return checked
            first_problem = checked
        elif dl_inspected.status != "missing":
            status = "inaccessible" if dl_inspected.status == "inaccessible" else "damaged"
            target = Path(dl_inspected.bundle) if dl_inspected.bundle else dl_store
            first_problem = ViralInspection(
                status=status, path=target, origin=None,
                datacache_status=dl_inspected.status, error=dl_inspected.error,
            )

    build_store = viral_store_path(cache_root, species, collection, area=_COLLECTIONS)
    build_inspected = datacache.inspect_materialization(build_store, verify_files=verify_files)
    if build_inspected.status == "available":
        bundle_dir = Path(build_inspected.bundle)
        checked = _validate_bundle_path(
            bundle_dir, species=species, collection=collection,
            origin="build", verify_files=verify_files,
        )
        if checked.status == "available":
            return checked
        first_problem = first_problem or checked
    elif build_inspected.status != "missing":
        status = "inaccessible" if build_inspected.status == "inaccessible" else "damaged"
        target = Path(build_inspected.bundle) if build_inspected.bundle else build_store
        first_problem = first_problem or ViralInspection(
            status=status, path=target, origin=None,
            datacache_status=build_inspected.status, error=build_inspected.error,
        )

    if first_problem is not None:
        return first_problem

    legacy = viral_store_path(cache_root, species, collection, area=_LEGACY_COLLECTIONS)
    if _legacy_directory(legacy):
        checked = _validate_bundle_path(
            legacy, species=species, collection=collection,
            origin="build", verify_files=verify_files,
        )
        return replace(checked, datacache_status=None)

    return ViralInspection(status="missing", path=None, origin=None, datacache_status="missing")


def _installed_bundle(cache_root, species, collection):
    """The directory holding the installed collection, or None when it isn't installed."""
    inspected = inspect_viral_collection(
        collection, species=species, cache_root=cache_root, verify_files=False,
    )
    if inspected.status == "available":
        return inspected.path
    if inspected.status == "missing":
        return None
    if inspected.status == "inaccessible":
        raise inspected.error or PermissionError(f"Cannot read viral collection at {inspected.path}")
    raise ValueError(
        f"The viral reference store {inspected.path} is damaged ({inspected.error}); "
        f"{_collection_repair(collection)}"
    ) from inspected.error


def _predict_record_outputs(collection, inventory, primaries, *, prefer_reference_molecule):
    """Return the expected relative ``records/*.gz`` paths for a collection build."""
    primary_source = collection.source()
    primary_suffix = getattr(primary_source, "document_suffix", ".gb")
    alt_source = collection.source(collection.alternate_source) if collection.alternate_source else None
    alt_suffix = getattr(alt_source, "document_suffix", ".json") if alt_source is not None else ""
    alt_ledger = (collection.source_summary(collection.alternate_source)
                  if collection.alternate_source else empty_ledger())
    outputs = []
    for entry, (_, primary) in zip(inventory.itertuples(index=False), primaries.iterrows(), strict=True):
        alternate = ledger_row(alt_ledger, entry.abbreviation)
        use_alternate = (
            prefer_reference_molecule
            and published(alternate)
            and alternate["sequence_sha256"] != primary["sequence_sha256"]
        )
        if use_alternate:
            filename = (
                f"{safe_token(entry.abbreviation)}__{alternate['source']}"
                f"__{safe_token(alternate['identifier'])}{alt_suffix}.gz"
            )
        else:
            filename = (
                f"{safe_token(entry.abbreviation)}__{collection.primary_source}"
                f"__{safe_token(primary['identifier'])}{primary_suffix}.gz"
            )
        outputs.append(f"{_RECORD_DIRECTORY}/{filename}")
    return outputs


class _DegradedOutputs(Exception):
    """Raised internally when an alternate source fell back to the primary document name."""

    def __init__(self, actual_files, fallback_dir):
        super().__init__("alternate degraded to primary record")
        self.actual_files = actual_files
        self.fallback_dir = fallback_dir


def _build_collection_bundle(
    store_path,
    collection,
    species,
    *,
    timeout,
    force,
    prefer_reference_molecule,
    annotate_from_alternate,
    prefer_alternate_annotation,
):
    """Build and atomically publish a viral collection via :func:`datacache.materialize`."""
    import datacache

    inventory = viral_panel_coverage(collection=collection.name, installed=False)
    primaries = _primary_ledger(collection, inventory)
    predicted_records = _predict_record_outputs(
        collection, inventory, primaries, prefer_reference_molecule=prefer_reference_molecule,
    )
    declared_outputs = {name: {} for name in (*_BASE_OUTPUTS, *predicted_records)}
    transform = {
        "version": BUNDLE_SCHEMA,
        "options": {
            "collection": collection.name,
            "bundle_version": bundle_version(collection),
            "prefer_reference_molecule": prefer_reference_molecule,
            "annotate_from_alternate": annotate_from_alternate,
            "prefer_alternate_annotation": prefer_alternate_annotation,
        },
    }
    spec_payload = json.dumps({
        "collection": collection.name,
        "bundle_version": bundle_version(collection),
        "panel": collection.panel,
        "primary_source": collection.primary_source,
        "alternate_source": collection.alternate_source,
    }, sort_keys=True, separators=(",", ":")).encode("utf-8")
    spec_sha256 = hashlib.sha256(spec_payload).hexdigest()

    with tempfile.TemporaryDirectory(prefix="orfdata-viral-spec-") as spec_dir:
        spec_file = Path(spec_dir) / "collection_spec.json"
        spec_file.write_bytes(spec_payload)
        sources = {
            "collection_spec.json": {
                "path": str(spec_file),
                "sha256": spec_sha256,
                "size": len(spec_payload),
                "source_identity": {
                    "collection": collection.name,
                    "bundle_version": bundle_version(collection),
                },
            }
        }

        def builder(_source_paths, output_paths):
            staged_dir = Path(output_paths[_METADATA]).parent
            _write_bundle(
                staged_dir, collection, inventory, primaries, timeout=timeout,
                prefer_reference_molecule=prefer_reference_molecule,
                annotate_from_alternate=annotate_from_alternate,
                prefer_alternate_annotation=prefer_alternate_annotation,
            )
            actual_files = _installed_files(staged_dir)
            if actual_files != set(output_paths):
                fallback_dir = Path(tempfile.mkdtemp(prefix="orfdata-viral-fallback-"))
                shutil.copytree(staged_dir, fallback_dir / "bundle")
                raise _DegradedOutputs(actual_files, fallback_dir)
            result = _existing_bundle(staged_dir, species=species, collection=collection)
            _verify_counts(result)

        try:
            paths = datacache.materialize(
                store_path,
                sources,
                transform=transform,
                outputs=declared_outputs,
                builder=builder,
                force=force,
                verify_files=False,
                download_options={"show_progress": False},
            )
        except _DegradedOutputs as degraded:
            fallback_dir = degraded.fallback_dir
            actual_outputs = {name: {} for name in sorted(degraded.actual_files)}
            try:
                def copy_degraded(_source_paths, output_paths):
                    staged_dir = Path(output_paths[_METADATA]).parent
                    src_root = fallback_dir / "bundle"
                    for rel_path, target in output_paths.items():
                        shutil.copy2(src_root / rel_path, target)
                    result = _existing_bundle(staged_dir, species=species, collection=collection)
                    _verify_counts(result)

                paths = datacache.materialize(
                    store_path,
                    sources,
                    transform=transform,
                    outputs=actual_outputs,
                    builder=copy_degraded,
                    force=True,
                    verify_files=False,
                    download_options={"show_progress": False},
                )
            finally:
                shutil.rmtree(fallback_dir, ignore_errors=True)

    published = Path(paths[_METADATA]).parent
    return _existing_bundle(published, species=species, collection=collection)


def _is_store_level_error(inspected, store_path):
    """Whether ``inspected.error`` reflects structural store damage rather than bundle tampering."""
    message = str(inspected.error or "")
    return any(token in message for token in (
        "not a datacache",
        "not a directory",
        "expected a directory",
        ".datacache-store.json",
    ))


def _discover_viral_stores(root, *, depth):
    """Discover candidate store directories at ``depth`` under ``root``, plus any readable DataCache stores."""
    from datacache.cli import find_bundle_stores

    stores = []
    skipped = {}
    if not _legacy_directory(root):
        return stores, skipped
    seen = set()
    for store_dir in find_bundle_stores(str(root)):
        path = Path(store_dir)
        seen.add(path)
        stores.append(path)

    # Also detect directories at the store depth that lost their marker or cannot be read.
    frontier = [root]
    for _ in range(depth):
        next_level = []
        for parent in frontier:
            if not _legacy_directory(parent):
                continue
            try:
                entries = sorted(parent.iterdir(), key=lambda entry: entry.name)
            except OSError as error:
                skipped[str(parent)] = f"can't be listed: {error}"
                continue
            for entry in entries:
                if entry.name.startswith("."):
                    continue
                next_level.append(entry)
        frontier = next_level
    for candidate in frontier:
        if candidate not in seen:
            stores.append(candidate)
    return sorted(stores), skipped


def _bundle_count(store_dir):
    """Count bundle directories currently present under a validated DataCache store."""
    bundles_dir = store_dir / "bundles"
    if not _is_directory(bundles_dir):
        return 0
    try:
        return sum(1 for entry in bundles_dir.iterdir() if _is_directory(entry))
    except OSError:
        return 0


def prune_viral_references(
    *,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    keep: int = 1,
    cache_root: str | PathLike[str] | None = None,
) -> dict[str, object]:
    """Delete all but the newest ``keep`` bundles of a collection and of its selections.

    A refresh keeps the bundle it replaces, as DataCache keeps its own,
    because a running program may still be reading it; nothing is pruned
    automatically. Paths into deleted bundles stop working, so prune only
    when no program still uses them. Reaches every bundle version and selection
    schema under the collection's DataCache stores, using DataCache's own
    locked pruning. A directory installed by ORFdata 0.43 or earlier counts as
    the oldest bundle. The PyEnsembl indexes :func:`viral_genome` built under
    ``cache_root`` for deleted bundles are removed too.

    Returns ``{"deleted": (bundle, ...), "skipped": {path: reason}}``. A store
    that can't be read or doesn't validate is skipped untouched, with the
    reason, and pruning goes on with the rest; so is a 0.43 directory that is
    a link, and a bundle that can't be removed.
    """
    import datacache
    from datacache.cli import prune_root

    if isinstance(keep, bool) or not isinstance(keep, int) or keep < 1:
        raise ValueError("keep must be a positive integer; the current bundle is always kept")
    species, collection = _selection(species, collection)
    viral_root = get_data_dir(cache_root) / "viral"
    deleted = []
    skipped = {}

    area_specs = (
        (viral_root / _COLLECTIONS / species / collection.name, 1, "materialization"),
        (viral_root / _DOWNLOADS / species / collection.name, 1, "bundle"),
        (viral_root / _SELECTIONS / species / collection.name, 2, "materialization"),
    )
    for area_root, depth, store_kind in area_specs:
        stores, area_skipped = _discover_viral_stores(area_root, depth=depth)
        skipped.update(area_skipped)
        for store_dir in stores:
            if store_kind == "bundle":
                inspected = datacache.inspect_bundle(store_dir, verify_files=False)
            else:
                inspected = datacache.inspect_materialization(store_dir, verify_files=False)
            if inspected.status == "inaccessible":
                skipped[str(store_dir)] = f"can't be read: {inspected.error}"
                continue
            if inspected.status == "invalid":
                skipped[str(store_dir)] = str(inspected.error)
                continue
            if inspected.status == "missing":
                continue
            try:
                results = prune_root(str(store_dir), keep=keep, apply=True)
            except (OSError, ValueError) as error:
                skipped[str(store_dir)] = f"couldn't be pruned: {error}"
                continue
            for entry in results:
                for removed in entry["removed"]:
                    bundle = Path(removed)
                    shutil.rmtree(_annotation_index(cache_root, collection, bundle), ignore_errors=True)
                    deleted.append(bundle)

    # Prune legacy 0.43 directories once the corresponding modern store has at least ``keep`` bundles.
    legacy_pairs = [
        (viral_store_path(cache_root, species, collection, area=_COLLECTIONS),
         viral_store_path(cache_root, species, collection, area=_LEGACY_COLLECTIONS)),
    ]
    legacy_selections_root = viral_store_path(cache_root, species, collection, area=_LEGACY_SELECTIONS)
    modern_selections_root = viral_store_path(cache_root, species, collection, area=_SELECTIONS)
    if _legacy_directory(legacy_selections_root):
        try:
            for entry in sorted(legacy_selections_root.iterdir(), key=lambda item: item.name):
                if not entry.name.startswith("."):
                    legacy_pairs.append((modern_selections_root / entry.name, entry))
        except OSError as error:
            skipped[str(legacy_selections_root)] = f"can't be listed: {error}"

    for modern_store, legacy in legacy_pairs:
        if not _legacy_directory(legacy):
            continue
        if _bundle_count(modern_store) < keep:
            continue
        if not _is_directory(legacy):
            skipped[str(legacy)] = "is a link to a 0.43 directory; remove it by hand"
            continue
        index = _annotation_index(cache_root, collection, legacy)
        try:
            shutil.rmtree(legacy)
        except OSError as error:
            skipped[str(legacy)] = f"couldn't be removed: {error}"
            continue
        shutil.rmtree(index, ignore_errors=True)
        deleted.append(legacy)

    return dict(deleted=tuple(deleted), skipped=skipped)


def _annotation_index(cache_root, collection, path):
    """Where :func:`viral_genome` keeps the PyEnsembl index of the files in ``path``."""
    collection = viral_collection(collection)
    digest = hashlib.sha256(str(Path(path).resolve()).encode("utf-8")).hexdigest()[:16]
    return (get_data_dir(cache_root) / "viral" / "annotations" /
            collection.name / bundle_version(collection) / digest).absolute()


def open_viral_references(
    *,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    download_if_missing: bool = False,
    cache_root: str | PathLike[str] | None = None,
    force: bool = False,
    timeout: float | None = DOWNLOAD_TIMEOUT,
) -> ViralReferenceFiles:
    """Open, or explicitly download and open, a viral source-reference collection.

    ``download_if_missing=True`` calls :func:`download_viral_references` when
    needed; otherwise this performs an offline, checksummed cache lookup.
    ``force=True`` explicitly refreshes or repairs the cache, which downloads
    a published artifact or rebuilds from pinned sources.
    ``species`` describes the target analysis species and ``collection`` names
    one of the published collections, which :func:`viral_collections` lists.
    Analysis contexts are nonexclusive metadata and never determine cache
    identity.
    """
    if type(download_if_missing) is not bool:
        raise ValueError("download_if_missing must be True or False")
    species, collection = _selection(species, collection)
    if download_if_missing or force:
        return download_viral_references(
            species=species, collection=collection.name, cache_root=cache_root,
            force=force, timeout=timeout)
    path = _installed_bundle(cache_root, species, collection)
    if path is None:
        raise NotInstalledError(
            f"The {collection.name} viral references are not installed; run "
            f"'orfdata data fetch {collection.name}', or pass "
            "download_if_missing=True."
        )
    return _verified_existing_bundle(path, species=species, collection=collection)


def download_viral_references(
    *,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    cache_root: str | PathLike[str] | None = None,
    force: bool = False,
    timeout: float | None = DOWNLOAD_TIMEOUT,
    prefer_reference_molecule: bool = True,
    annotate_from_alternate: bool = True,
    prefer_alternate_annotation: bool = False,
    source_paths: dict[str, str | PathLike[str]] | None = None,
    download_options: dict[str, object] | None = None,
) -> ViralReferenceFiles:
    """Install a viral genome, transcript and protein collection for an analysis target.

    Downloads a pinned release artifact through DataCache's
    :class:`~datacache.VersionedDatasetRegistry` when one is registered and
    default content options are used; otherwise builds the collection locally
    from pinned source ledgers and publishes it atomically through
    :func:`datacache.materialize`. Passing any option that changes content
    (``prefer_reference_molecule=False``, ``annotate_from_alternate=False``,
    or ``prefer_alternate_annotation=True``) selects the local build route.

    :func:`viral_collections` lists what can be installed; the one published
    today is ``gdc-viruses``, the 200-reference GDC screen for human
    cancer-transcriptome analysis. Target species does not assert that every
    source virus has that host, collection membership is not a per-virus
    causality classification, and analysis contexts are not exclusive. To
    create a filtered derived set, use :func:`select_viral_references`. The
    collection is stored under ``ORFDATA_CACHE`` (or the normal OpenVax cache);
    a verified cache is reused and ``force=True`` refreshes or repairs it. A
    refresh is published beside the installed collection as a new DataCache
    generation, so a reader finds the old files or the new ones, never a gap;
    the current directory therefore changes with each refresh. Replaced copies
    stay for programs still reading them until :func:`prune_viral_references`
    or ``orfdata data prune`` removes them, and installs of one collection take
    turns. A request that stalls for ``timeout`` seconds is retried, then
    fails; ``timeout=None`` waits.

    This source-record collection is separate from the curated
    :func:`viral_references` catalogue.
    """
    species, collection = _selection(species, collection)
    if type(force) is not bool:
        raise ValueError("force must be True or False")
    if timeout is not None and (not isinstance(timeout, (int, float))
                                or isinstance(timeout, bool) or timeout <= 0):
        raise ValueError("timeout must be a positive number of seconds, or None to wait")

    default_content = (
        prefer_reference_molecule
        and annotate_from_alternate
        and not prefer_alternate_annotation
    )
    use_download = collection.published and default_content

    if use_download:
        registry = viral_download_registry(collection, cache_root)
        version = bundle_version(collection)
        dl_store = viral_store_path(cache_root, species, collection, area=_DOWNLOADS)
        if not force:
            inspected = inspect_viral_collection(
                collection, species=species, cache_root=cache_root, verify_files=True,
            )
            if inspected.status == "available":
                return _existing_bundle(inspected.path, species=species, collection=collection)
            if inspected.status == "inaccessible":
                raise ValueError(
                    f"The {collection.name} viral references at {inspected.path} can't be read "
                    f"({inspected.error}); fix their permissions, or replace them with force=True "
                    "(--force on the command line)"
                ) from inspected.error
            if inspected.status == "damaged":
                if _is_store_level_error(inspected, dl_store):
                    raise ValueError(
                        f"The viral reference store {dl_store} is damaged ({inspected.error}); "
                        f"{_collection_repair(collection)}"
                    ) from inspected.error
                raise FileNotFoundError(
                    f"The {collection.name} viral references at {inspected.path} are incomplete or invalid; "
                    "rebuild them with force=True (--force on the command line)"
                ) from inspected.error
        opts = dict(download_options or {})
        opts.setdefault("timeout", timeout)
        overrides = None if source_paths is None else {k: str(v) for k, v in source_paths.items()}
        registry.download(collection.name, version, force=force, source_paths=overrides, **opts)
        after = registry.inspect(collection.name, version, verify_files=False)
        return _existing_bundle(Path(after.bundle), species=species, collection=collection)

    build_store = viral_store_path(cache_root, species, collection, area=_COLLECTIONS)
    materialize_force = force
    if not force:
        if default_content:
            inspected = inspect_viral_collection(
                collection, species=species, cache_root=cache_root, verify_files=True,
            )
            if inspected.status == "available":
                return _existing_bundle(inspected.path, species=species, collection=collection)
            if inspected.status == "inaccessible":
                raise ValueError(
                    f"The {collection.name} viral references at {inspected.path} can't be read "
                    f"({inspected.error}); fix their permissions, or replace them with force=True "
                    "(--force on the command line)"
                ) from inspected.error
            if inspected.status == "damaged":
                if _is_store_level_error(inspected, build_store):
                    raise ValueError(
                        f"The viral reference store {build_store} is damaged ({inspected.error}); "
                        f"{_collection_repair(collection)}"
                    ) from inspected.error
                raise FileNotFoundError(
                    f"The {collection.name} viral references at {inspected.path} are incomplete or invalid; "
                    "rebuild them with force=True (--force on the command line)"
                ) from inspected.error
        else:
            import datacache

            desired_options = {
                "collection": collection.name,
                "bundle_version": bundle_version(collection),
                "prefer_reference_molecule": prefer_reference_molecule,
                "annotate_from_alternate": annotate_from_alternate,
                "prefer_alternate_annotation": prefer_alternate_annotation,
            }
            existing = datacache.inspect_materialization(build_store, verify_files=True)
            if existing.status == "available":
                if (isinstance(existing.transform, dict)
                        and existing.transform.get("options") == desired_options):
                    return _verified_existing_bundle(
                        Path(existing.bundle), species=species, collection=collection,
                    )
                materialize_force = True
            elif existing.status == "inaccessible":
                raise ValueError(
                    f"The {collection.name} viral references at {build_store} can't be read "
                    f"({existing.error}); fix their permissions, or replace them with force=True "
                    "(--force on the command line)"
                ) from existing.error
            elif existing.status == "invalid":
                raise FileNotFoundError(
                    f"The {collection.name} viral references at {build_store} are incomplete or invalid; "
                    "rebuild them with force=True (--force on the command line)"
                ) from existing.error

    return _build_collection_bundle(
        build_store,
        collection,
        species,
        timeout=timeout,
        force=materialize_force,
        prefer_reference_molecule=prefer_reference_molecule,
        annotate_from_alternate=annotate_from_alternate,
        prefer_alternate_annotation=prefer_alternate_annotation,
    )


def _copy_fasta(source, destination, keys):
    """Copy the records whose identifier is selected, keeping their order."""
    keep, count = False, 0
    with source.open() as input_stream, destination.open("w") as output_stream:
        for line in input_stream:
            if line.startswith(">"):
                keep = line[1:].split(maxsplit=1)[0] in keys
                count += keep
            if keep:
                output_stream.write(line)
    return count


def _copy_table(source, destination, column, values):
    rows = [row for row in _read_table(source) if row[column] in values]
    with source.open(newline="") as stream:
        columns = tuple(csv.DictReader(stream, delimiter="\t").fieldnames or ())
    _write_table(destination, columns, rows)
    return len(rows)


def _copy_gtf(source, destination, genome_ids):
    """Keep the header and every row on a selected viral genome."""
    with source.open() as input_stream, destination.open("w") as output_stream:
        for line in input_stream:
            if line.startswith("#") or line.split("\t", 1)[0] in genome_ids:
                output_stream.write(line)


def _index_path(cache_root, species, collection, contexts, excluded_contexts,
                context_match, subsets, genome_ids, *, area=_SELECTIONS):
    """The store a selection is published in; ``area`` selects the 0.43 directory."""
    version = bundle_version(viral_collection(collection))
    identity = json.dumps({
        "schema": INDEX_SCHEMA,
        "bundle_version": version,
        "contexts": contexts,
        "excluded_contexts": excluded_contexts,
        "context_match": context_match,
        "subsets": subsets,
        "genomes": genome_ids,
    }, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    return (get_data_dir(cache_root) / "viral" / area / species / collection /
            version / digest)


def _index_result(path, *, bundle, contexts, excluded_contexts, context_match, subsets):
    metadata = _metadata(path, expected_schema=INDEX_SCHEMA)
    if metadata.get("bundle_version") != bundle_version(viral_collection(bundle.collection)):
        raise ValueError(f"Viral reference index identity mismatch: {path}")
    expected = {
        "target_species": bundle.target_species,
        "target_taxon_id": bundle.target_taxon_id,
        "collection": bundle.collection,
        "contexts": list(contexts),
        "excluded_contexts": list(excluded_contexts),
        "context_match": context_match,
        "subsets": list(subsets),
    }
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError(f"Viral reference index identity mismatch: {path}")
    source = Path(str(metadata.get("source_bundle")))
    if _resolved(source) != _resolved(bundle.path):
        # Built from a collection bundle a refresh has since replaced: never
        # reused under the current collection's name, so it is rebuilt.
        raise ValueError(f"Viral reference index {path} was built from {source}, "
                         f"not the current collection {bundle.path}")
    return SelectedViralReferences(
        target_species=bundle.target_species,
        target_taxon_id=bundle.target_taxon_id,
        collection=bundle.collection,
        genome_count=metadata["genome_count"],
        transcript_count=metadata["transcript_count"],
        protein_count=metadata["protein_count"],
        contexts=contexts,
        excluded_contexts=excluded_contexts,
        context_match=context_match,
        subsets=subsets,
        source_collection=source,
        **_paths(path),
    )


def _resolved(path):
    try:
        return Path(path).resolve()
    except (OSError, RuntimeError):
        return Path(path).absolute()


def _verified_existing_index(path, *, bundle, contexts, excluded_contexts, context_match, subsets):
    result = _index_result(
        path, bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
        context_match=context_match, subsets=subsets)
    required = (result.genome_fasta, result.transcript_fasta, result.protein_fasta,
                result.annotation_gtf, result.genome_table, result.transcript_table,
                path / _CHECKSUMS)
    if any(not item.is_file() for item in required):
        raise ValueError(f"Incomplete viral reference index: {path}")
    _verify_checksums(path)
    _verify_counts(result)
    return result


def _selected_abbreviations(contexts, excluded_contexts, context_match, subsets):
    """The inventory entries a context and subset selection resolves to."""
    inventory = viral_panel_coverage(
        contexts=contexts or None, exclude_contexts=excluded_contexts or None,
        context_match=context_match, installed=False)
    if subsets:
        selected = list(viral_subset_abbreviations(subsets))
        inventory = inventory[inventory["abbreviation"].isin(selected)]
    return tuple(inventory["abbreviation"])


def _selected_genomes(installed, abbreviations):
    """The genome identifiers an installed collection carries for those entries.

    A collection names its genomes after the molecules it carries, which is not
    always the identifier the inventory lists first, so the installed table is
    what a selection resolves against.
    """
    rows = {row["abbreviation"]: row["genome_id"] for row in _read_table(installed.genome_table)}
    missing = [name for name in abbreviations if name not in rows]
    if missing:
        raise ValueError(
            f"The installed collection has no genome for {', '.join(missing[:5])}")
    return tuple(rows[name] for name in abbreviations)


def select_viral_references(
    *,
    contexts: Iterable[str] | None = None,
    exclude_contexts: Iterable[str] | None = None,
    context_match: str = "any",
    subsets: Iterable[str] | None = None,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    download_if_missing: bool = False,
    cache_root: str | PathLike[str] | None = None,
    force: bool = False,
    timeout: float | None = DOWNLOAD_TIMEOUT,
) -> ViralReferenceFiles:
    """Select reusable reference files by viral context or named subset.

    The current collection targets human analysis; future collections may
    target other species. Context membership is an association and never
    changes a FASTA record's accession-based identity. Repeated ``contexts``
    select their union by default or their intersection with
    ``context_match='all'``; excluded contexts are then removed. Excluding a
    context without naming one starts from every context. The source
    GenBank collection is cached once; selections reuse it rather than downloading the
    same reference for each context.

    The derived set carries the same genome, transcriptome and proteome layers
    as the collection, including its GTF. ``force=True`` rebuilds the selection from
    the collection, and also repairs the collection itself when ``download_if_missing``
    allows the network. A damaged selection, or a selection store that lost its
    marker, is rebuilt either way, since it is derived from local files; a path
    in the store's place that ORFdata didn't write is reported, never removed.
    Each rebuild is a new bundle, so the selection's directory changes, and
    replaced ones stay until :func:`prune_viral_references` removes them.
    """
    import datacache

    contexts = normalize_viral_contexts(contexts)
    excluded_contexts = normalize_viral_contexts(exclude_contexts)
    subsets = normalize_viral_subsets(subsets)
    if not contexts and not subsets and not excluded_contexts:
        raise ValueError("At least one viral context, subset or excluded context is required")
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    bundle = open_viral_references(
        species=species, collection=collection, download_if_missing=download_if_missing,
        cache_root=cache_root, force=force and download_if_missing, timeout=timeout)
    abbreviations = _selected_abbreviations(contexts, excluded_contexts, context_match, subsets)
    if not abbreviations:
        raise ValueError("No viral references match the selected contexts and subsets")
    genome_ids = _selected_genomes(bundle, abbreviations)
    selected = set(genome_ids)
    location = (cache_root, bundle.target_species, bundle.collection, contexts,
                excluded_contexts, context_match, subsets, genome_ids)
    identity = dict(bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
                    context_match=context_match, subsets=subsets)
    store_path = _index_path(*location)
    legacy = _index_path(*location, area=_LEGACY_SELECTIONS)

    if not force:
        inspected = datacache.inspect_materialization(store_path, verify_files=True)
        if inspected.status == "available":
            try:
                return _index_result(Path(inspected.bundle), **identity)
            except (OSError, ValueError):
                pass
        elif inspected.status == "missing" and _legacy_directory(legacy):
            try:
                return _verified_existing_index(legacy, **identity)
            except (OSError, ValueError):
                pass

    sources = {
        _GENOME_FASTA: {"path": str(bundle.genome_fasta)},
        _TRANSCRIPT_FASTA: {"path": str(bundle.transcript_fasta)},
        _PROTEIN_FASTA: {"path": str(bundle.protein_fasta)},
        _ANNOTATION_GTF: {"path": str(bundle.annotation_gtf)},
        _GENOME_TABLE: {"path": str(bundle.genome_table)},
        _TRANSCRIPT_TABLE: {"path": str(bundle.transcript_table)},
    }
    version = bundle_version(viral_collection(bundle.collection))
    transform = {
        "version": INDEX_SCHEMA,
        "options": {
            "bundle_version": version,
            "collection": bundle.collection,
            "target_species": bundle.target_species,
            "target_taxon_id": bundle.target_taxon_id,
            "contexts": list(contexts),
            "excluded_contexts": list(excluded_contexts),
            "context_match": context_match,
            "subsets": list(subsets),
            "genome_ids": list(genome_ids),
            "source_bundle": str(bundle.path),
        },
    }
    declared_outputs = {name: {} for name in _BASE_OUTPUTS}

    def build(source_paths, output_paths):
        staging = Path(output_paths[_METADATA]).parent
        transcript_table = Path(source_paths[_TRANSCRIPT_TABLE])
        genome_table = Path(source_paths[_GENOME_TABLE])
        transcripts = [row for row in _read_table(transcript_table)
                       if row["genome_id"] in selected]
        genome_count = _copy_fasta(Path(source_paths[_GENOME_FASTA]), staging / _GENOME_FASTA, selected)
        transcript_count = _copy_fasta(
            Path(source_paths[_TRANSCRIPT_FASTA]), staging / _TRANSCRIPT_FASTA,
            {row["transcript_id"] for row in transcripts})
        protein_count = _copy_fasta(
            Path(source_paths[_PROTEIN_FASTA]), staging / _PROTEIN_FASTA,
            {row["protein_id"] for row in transcripts})
        genome_rows = _copy_table(genome_table, staging / _GENOME_TABLE,
                                  "genome_id", selected)
        transcript_rows = _copy_table(transcript_table, staging / _TRANSCRIPT_TABLE,
                                      "genome_id", selected)
        _copy_gtf(Path(source_paths[_ANNOTATION_GTF]), staging / _ANNOTATION_GTF, selected)
        if (genome_count != genome_rows or genome_count != len(selected)
                or transcript_count != transcript_rows or protein_count != transcript_rows):
            raise ValueError("Viral reference index counts disagree")
        (staging / _METADATA).write_text(json.dumps({
            "schema": INDEX_SCHEMA,
            "bundle_version": version,
            "collection": bundle.collection,
            "target_species": bundle.target_species,
            "target_taxon_id": bundle.target_taxon_id,
            "contexts": list(contexts),
            "excluded_contexts": list(excluded_contexts),
            "context_match": context_match,
            "subsets": list(subsets),
            "genome_ids": list(genome_ids),
            "source_bundle": str(bundle.path),
            "genome_count": genome_count,
            "transcript_count": transcript_count,
            "protein_count": protein_count,
            "annotation_sources": dict(sorted(
                Counter(row["annotation_source"] for row in transcripts).items())),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (staging / "README.md").write_text(
            "# ORFdata viral reference index\n\n"
            f"Contexts ({context_match}): {', '.join(contexts) or 'all'}\n\n"
            f"Excluded contexts: {', '.join(excluded_contexts) or 'none'}\n\n"
            f"Named subsets: {', '.join(subsets) or 'none'}\n\n"
            f"Genome, transcriptome and proteome layers for {genome_count} of the "
            f"{bundle.genome_count} references in the `{bundle.collection}` collection, "
            f"version `{version}`, with the same "
            "identifiers and the matching GTF rows. Context, subset and collection membership "
            "are evidence associations and do not alter identity.\n",
            encoding="utf-8",
        )
        _write_checksums(staging)
        result = _index_result(staging, **identity)
        _verify_counts(result)

    paths = datacache.materialize(
        store_path,
        sources,
        transform=transform,
        outputs=declared_outputs,
        builder=build,
        force=True,
        verify_files=False,
        download_options={"show_progress": False},
    )
    return _index_result(Path(paths[_METADATA]).parent, **identity)


def viral_reference_summary(*, contexts=None, exclude_contexts=None, context_match="any",
                            subsets=None, species=None, collection=None,
                            cache_root=None):
    """Describe an installed collection, or one selection of it, without writing anything.

    The result is a plain dictionary: the collection's identity and cache
    path, its genome, transcriptome and proteome layers with their record
    counts, how many coding transcripts carry whole start and stop codons,
    which source annotated them, and which selected genomes have no annotated
    CDS feature at all.
    A context or subset selection narrows the counts to the matching genomes;
    it does not build a derived selection, which :func:`select_viral_references`
    does.
    """
    contexts = normalize_viral_contexts(contexts)
    excluded_contexts = normalize_viral_contexts(exclude_contexts)
    subsets = normalize_viral_subsets(subsets)
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    bundle = open_viral_references(species=species, collection=collection, cache_root=cache_root)
    metadata = _metadata(bundle.path, expected_schema=BUNDLE_SCHEMA)
    genomes = _read_table(bundle.genome_table)
    transcripts = _read_table(bundle.transcript_table)
    selected = bool(contexts or excluded_contexts or subsets)
    if selected:
        chosen = set(_selected_abbreviations(contexts, excluded_contexts, context_match, subsets))
        genomes = [row for row in genomes if row["abbreviation"] in chosen]
        carried = {row["genome_id"] for row in genomes}
        transcripts = [row for row in transcripts if row["genome_id"] in carried]
    counts = (len(genomes), len(transcripts), len(transcripts))
    return {
        "collection": bundle.collection,
        "bundle_version": bundle_version(viral_collection(bundle.collection)),
        "panel": bundle.panel,
        "target_species": bundle.target_species,
        "target_taxon_id": bundle.target_taxon_id,
        "retrieved": bundle.retrieved,
        "path": str(bundle.path),
        "selection": {
            "selected": selected,
            "contexts": list(contexts),
            "excluded_contexts": list(excluded_contexts),
            "context_match": context_match,
            "subsets": list(subsets),
        },
        "genome_count": counts[0],
        "transcript_count": counts[1],
        "protein_count": counts[2],
        "complete_transcript_count": sum(row["annotation_complete"] == "yes"
                                         for row in transcripts),
        # How the selected transcripts begin, and which ones the source says
        # the ribosome recodes rather than reading straight through.
        "initiations": dict(sorted(Counter(row["initiation"] for row in transcripts).items())),
        "non_canonical_starts": [row["transcript_id"] for row in transcripts
                                 if row["initiation"] == NON_CANONICAL_START],
        "recoded_transcripts": {row["transcript_id"]: row["recoding"]
                                for row in transcripts if row["recoding"]},
        # Which source annotated the selected transcripts.
        "annotation_sources": dict(sorted(
            Counter(row["annotation_source"] for row in transcripts).items())),
        # A source record can carry no CDS feature at all: its genome is still
        # selectable, but it contributes no transcript and no protein.
        "genomes_without_transcripts": [row["abbreviation"] for row in genomes
                                        if int(row["transcript_count"]) == 0],
        # Where filling that gap was tried and did not check out.
        "annotation_fill_skipped": {name: reason for name, reason
                                    in metadata.get("annotation_fill_skipped", {}).items()
                                    if any(row["abbreviation"] == name for row in genomes)},
        # Where each molecule came from, and any entry still carrying a record
        # that is not the alternate reference its inventory row names.
        "genome_sources": dict(sorted(Counter(row["genome_source"] for row in genomes).items())),
        "genomes_unlike_alternate": [row["abbreviation"] for row in genomes
                                     if row.get("alternate_sequence") == "differs"],
        "layers": [dict(layer.to_dict(), records=count)
                   for layer, count in zip(bundle.layers, counts, strict=True)],
    }


def viral_annotation_survey(*, contexts=None, exclude_contexts=None, context_match="any",
                            subsets=None, species=None, collection=None,
                            cache_root=None):
    """Say, genome by genome, how complete an installed collection's annotation is.

    A reference inventory is only as good as the submissions behind it, and
    those differ in how much of a genome they annotate: some name every gene,
    some leave the name in free text, some annotate no coding feature at all,
    and few annotate spliced products. This walks the installed collection and
    reports what each genome actually carries.

    Where a collection names an alternate source, that source's packaged ledger
    records how many coding features it publishes for the same entry, so a
    genome can be measured against it without fetching anything: ``shortfall``
    is how many more coding features the alternate publishes, and is only
    counted where the alternate publishes the identical molecule, so it is a
    difference of annotation and never of sequence. No gene list is built in,
    which is what keeps this true of any collection, species and source.
    """
    summary = viral_reference_summary(
        contexts=contexts, exclude_contexts=exclude_contexts, context_match=context_match,
        subsets=subsets, species=species, collection=collection, cache_root=cache_root)
    bundle = open_viral_references(species=species, collection=collection, cache_root=cache_root)
    chosen = viral_collection(bundle.collection)
    ledger = (chosen.source_summary(chosen.alternate_source)
              if chosen.alternate_source else empty_ledger())
    genomes = _read_table(bundle.genome_table)
    if summary["selection"]["selected"]:
        wanted = {row["abbreviation"] for row in genomes} & set(_selected_abbreviations(
            normalize_viral_contexts(contexts), normalize_viral_contexts(exclude_contexts),
            context_match, normalize_viral_subsets(subsets)))
        genomes = [row for row in genomes if row["abbreviation"] in wanted]
    carried = {row["genome_id"] for row in genomes}
    transcripts = defaultdict(list)
    for row in _read_table(bundle.transcript_table):
        if row["genome_id"] in carried:
            transcripts[row["abbreviation"]].append(row)

    surveyed = []
    for genome in genomes:
        rows = transcripts.get(genome["abbreviation"], [])
        alternate = ledger_row(ledger, genome["abbreviation"])
        # Only an alternate publishing this very molecule can be compared: a
        # different molecule's feature count says nothing about this one.
        comparable = published(alternate) and genome["alternate_sequence"] == "same"
        surveyed.append({
            "abbreviation": genome["abbreviation"],
            "genome_id": genome["genome_id"],
            "genome_source": genome["genome_source"],
            "annotation_source": genome["annotation_source"],
            "transcripts": len(rows),
            "complete": sum(row["annotation_complete"] == "yes" for row in rows),
            # A reading frame the source does not close, whether or not it also
            # marks an end as running past the record.
            "incomplete": sum(row["annotation_complete"] != "yes" for row in rows),
            "partial_ends": sum(bool(row["partial"]) for row in rows),
            # A gene whose name is only in free-text prose cannot be looked up.
            "unnamed_genes": sum(not row["gene"] for row in rows),
            "spliced": sum(int(row["exon_count"]) > 1 for row in rows),
            "non_canonical_starts": sum(row["initiation"] == NON_CANONICAL_START for row in rows),
            "recoded": sum(bool(row["recoding"]) for row in rows),
            "alternate_source": genome["alternate_source"],
            "alternate_features": int(alternate["coding_features"]) if comparable else None,
            "shortfall": max(int(alternate["coding_features"]) - len(rows), 0)
            if comparable else None,
        })
    return dict(summary, genomes=surveyed, totals={
        "genomes": len(surveyed),
        "genomes_fully_annotated": sum(
            row["transcripts"] > 0 and row["incomplete"] == 0 and row["unnamed_genes"] == 0
            and not row["shortfall"] for row in surveyed),
        "genomes_with_unnamed_genes": sum(bool(row["unnamed_genes"]) for row in surveyed),
        "genomes_with_incomplete_transcripts": sum(bool(row["incomplete"])
                                                   for row in surveyed),
        "genomes_with_spliced_transcripts": sum(bool(row["spliced"]) for row in surveyed),
        "genomes_behind_alternate": sum(bool(row["shortfall"]) for row in surveyed),
        "coding_features_behind_alternate": sum(row["shortfall"] or 0 for row in surveyed),
    })


def viral_reference_context(reference_set: ViralReferenceFiles):
    """The :class:`ReferenceContext` that declares what :func:`viral_genome` opens.

    A PyEnsembl genome over viral contigs carries no species, so consumers that
    check a host reference against the analysis (vaxrank's ``reference_context``)
    need it stated. The context names the collection's *target* species, the
    analysis the set is for; it asserts no source host. It validates against the
    genome ``viral_genome(reference_set)`` returns.
    """
    from .reference import ReferenceContext

    collection = viral_collection(reference_set.collection)
    return ReferenceContext(
        taxon_id=reference_set.target_taxon_id,
        species=canonical_species(reference_set.target_species),
        genome_build=f"{collection.panel}-viral",
        annotation_name=f"orfdata-{collection.name}",
        annotation_version=bundle_version(collection))


def viral_genome(
    reference_set: ViralReferenceFiles,
    *,
    cache_root: str | PathLike[str] | None = None,
    rebuild: bool = False,
):
    """Open viral reference files as a PyEnsembl genome of coding transcripts.

    The GTF, transcript FASTA and protein FASTA of the set are indexed into a
    PyEnsembl database kept beside the cache, never inside the checksummed
    reference directory. varcode, isovar, topiary and vaxrank take the result
    wherever they take a PyEnsembl genome, so a variant called on a viral
    contig gets its coding effect from the source CDS annotation, including a
    non-standard genetic code declared by the source record.
    """
    if not isinstance(reference_set, ViralReferenceFiles):
        raise TypeError("viral_genome needs ViralReferenceFiles")
    from pyensembl import Genome

    collection = viral_collection(reference_set.collection)
    version = bundle_version(collection)
    # Each published bundle has its own directory, so its index is never one
    # built from a bundle it replaced.
    index = _annotation_index(cache_root, collection, reference_set.path)
    if rebuild:
        shutil.rmtree(index, ignore_errors=True)
    genome = Genome(
        reference_name=f"{collection.panel}-viral",
        annotation_name=f"orfdata-{collection.name}",
        annotation_version=version,
        gtf_path_or_url=str(Path(reference_set.annotation_gtf).absolute()),
        transcript_fasta_paths_or_urls=[str(Path(reference_set.transcript_fasta).absolute())],
        protein_fasta_paths_or_urls=[str(Path(reference_set.protein_fasta).absolute())],
        genome_fasta_path_or_url=str(Path(reference_set.genome_fasta).absolute()),
        cache_directory_path=str(index),
    )
    genome.index()  # Local files only; done once.
    return genome
