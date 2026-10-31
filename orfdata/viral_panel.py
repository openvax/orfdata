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
import shutil
import stat
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from os import PathLike
from pathlib import Path
from uuid import uuid4

from .fasta import write_fasta
from .record_sources import REQUEST_INTERVAL
from .reference import canonical_species
from .reference_ledger import empty_ledger, ledger_row, published, same_molecule, verify_record
from .viral_annotation import (
    NON_CANONICAL_START,
    gtf_gene_line,
    gtf_transcript_lines,
    name_transcripts,
    safe_token,
)
from .viral_coverage import normalize_viral_contexts, viral_panel_coverage
from .viral_registry import ViralCollection, viral_collection
from .viral_subsets import normalize_viral_subsets, viral_subset_abbreviations

BUNDLE_SCHEMA_VERSION = "6"
BUNDLE_SCHEMA = "orfdata.viral_reference_bundle.v6"
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
    "source_refseq", "viral_taxon_id", "source_host_taxon_ids", "length_nt", "genome_sha256",
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


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
                                 and item.name != _CHECKSUMS),
                                key=lambda item: item.relative_to(path).as_posix()):
            stream.write(f"{_sha256(file_path)}  {file_path.relative_to(path).as_posix()}\n")


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
    if (not result.record_directory.is_dir()
            or len(tuple(result.record_directory.glob("*.gz"))) != result.genome_count):
        raise ValueError(f"Incomplete viral source records: {result.record_directory}")
    return result


def _verify_checksums(path):
    sums_path = path / _CHECKSUMS
    try:
        sums = tuple(line.rstrip("\n").split("  ", 1) for line in sums_path.open())
    except OSError as error:
        raise ValueError(f"Missing viral reference checksums: {sums_path}") from error
    expected_files = {item.relative_to(path).as_posix() for item in path.rglob("*")
                      if item.is_file() and item.name != _CHECKSUMS}
    if (len(sums) != len(expected_files) or any(len(parts) != 2 for parts in sums)
            or {relative for _, relative in sums} != expected_files):
        raise ValueError(f"Invalid viral reference checksum manifest: {sums_path}")
    if any(_sha256(path / relative) != digest for digest, relative in sums):
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


def _genome_row(entry, primary, molecule, transcripts, alternate):
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
        # The target species is an analysis choice, not a claim about the source host.
        "source_host_taxon_ids": "",
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
    records_directory.mkdir()
    source = collection.source()
    ledger = (collection.source_summary(collection.alternate_source)
              if collection.alternate_source else empty_ledger())
    entries = list(inventory.itertuples(index=False))
    primaries = [row for _, row in primaries.iterrows()]
    pinned = {primary["identifier"]: primary for primary in primaries}
    downloaded, batch_size = {}, getattr(source, "batch_size", 1)
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
            time.sleep(REQUEST_INTERVAL)  # Respect the source's unauthenticated request rate.
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
    protein_uses = Counter(feature["source_protein_id"] for molecule in molecules.values()
                           for feature in molecule.features if feature["source_protein_id"])
    unique_protein_ids = {value for value, uses in protein_uses.items() if uses == 1}

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
                                           ledger_row(ledger, entry.abbreviation)))

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
        "target does not assert a source host; `source_host_taxon_ids` is empty when the source "
        "inventory does not establish one.\n\n"
        "This is a source-record download, not an ORFdata curated viral-reference release. "
        "It does not alter `viral_references()` or `Proteome.build(..., viral=...)`.\n"
    )


# Collections and selections are published the way DataCache publishes its own
# bundles. Each version has a store directory holding complete bundles under
# bundles/, the newest of which is current. An install builds a whole new bundle
# beside the current one, under the store's lock, and renames it into place, so a
# reader finds the old bundle or the new one, never a missing or partial one.
# Replaced bundles are kept, as DataCache keeps its own, because a running
# program may still be reading them; prune_viral_references removes them. ORFdata
# 0.44 and earlier kept one plain directory per version, under "references" and
# "indexes"; such a directory is still read while the store has no bundle, and
# counts as the oldest bundle.
_STORE_KIND = "orfdata-viral"
_COLLECTIONS, _LEGACY_COLLECTIONS = "collections", "references"
_SELECTIONS, _LEGACY_SELECTIONS = "selections", "indexes"


def _bundle_store():
    """DataCache's store primitives: its lock, bundle naming and private working directories.

    ``datacache.bundle_store`` is not part of DataCache's exported API, so a
    DataCache release could move it; this is the one place ORFdata imports it,
    and it fails naming the cause. Publishing viral collections through
    DataCache's public registries (#159) removes the dependency.
    """
    try:
        from datacache import bundle_store

        for name in ("BundleStore", "FORMAT", "MARKER", "private_directory",
                     "discard_private_directory", "user_key"):
            getattr(bundle_store, name)
    except (ImportError, AttributeError) as error:
        raise ImportError(
            "ORFdata's viral collections need datacache.bundle_store, as DataCache 1.29.3 "
            f"provides it, and the installed DataCache doesn't ({error})") from error
    return bundle_store


def _bundle_path(cache_root, species, collection, *, area=_COLLECTIONS):
    """The store a collection's version is published in; ``area`` selects the 0.44 directory."""
    from .datasets import get_data_dir

    collection = viral_collection(collection)
    return (get_data_dir(cache_root) / "viral" / area / species /
            collection.name / bundle_version(collection))


def _store(path):
    return _bundle_store().BundleStore(path, _STORE_KIND)


def _is_directory(path):
    """Whether ``path`` is a real directory; an unreadable one, as in a shared cache, counts as absent."""
    try:
        return stat.S_ISDIR(Path(path).lstat().st_mode)
    except OSError:
        return False


def _store_problem(store):
    """What keeps a store from being read, as ``(reason, target, repairable)``; None when it reads.

    A forced install repairs a store whose marker is lost but which holds only
    what a viral store keeps, and one whose newest bundle isn't a directory,
    by publishing a newer bundle. Anything else is never taken over: the
    ``target`` path has to be removed by hand. An unreadable store raises.
    """
    path, marker = store.path, store.path / _bundle_store().MARKER
    try:
        info = path.lstat()
    except FileNotFoundError:
        return None
    if not stat.S_ISDIR(info.st_mode):
        return "is not a directory", path, False
    entries = [entry.name for entry in path.iterdir()]
    try:
        found = json.loads(marker.read_text())
    except FileNotFoundError:
        if all(name.startswith(".datacache-") for name in entries):
            return None  # Empty: the first install creates the store here.
        found = None
    except IsADirectoryError:
        return "keeps its marker in something that isn't a file", marker, False
    except (ValueError, RecursionError):
        found = None
    if found !={"format": _bundle_store().FORMAT, "kind": _STORE_KIND}:
        ours = all(name in ("bundles", marker.name) or name.startswith((".staging-", ".datacache-"))
                   for name in entries)
        if found is None and ours:
            return f"has lost its {marker.name} marker", path, True
        return "is not an ORFdata viral store", path, False
    try:
        bundles = store.bundles.lstat()
    except FileNotFoundError:
        return None
    if not stat.S_ISDIR(bundles.st_mode):
        return "keeps its bundles in something that isn't a directory", store.bundles, False
    names = store.bundle_names()
    if names and not _is_directory(store.bundles / names[-1]):
        return f"has a newest bundle, {names[-1]}, that isn't a directory", path, True
    return None


def _damaged(store, problem):
    reason, target, repairable = problem
    remedy = ("replace it with force=True" if repairable else
              f"ORFdata never takes over or deletes such a path; remove {target} and install again")
    return ValueError(f"The viral reference store {store.path} {reason}; {remedy}")


def _current(path, legacy):
    """The directory holding a store's current bundle, else a 0.44 directory, else None.

    A damaged store raises ValueError naming what fixes it, an unreadable one
    its OSError. An unreadable 0.44 directory counts as absent.
    """
    store = _store(path)
    problem = _store_problem(store)
    if problem is not None:
        raise _damaged(store, problem)
    current = store.current_bundle()
    if current is None and _is_directory(legacy):
        return legacy
    return current


def _installed_bundle(cache_root, species, collection):
    """The directory holding the installed collection, or None when it isn't installed."""
    return _current(_bundle_path(cache_root, species, collection),
                    _bundle_path(cache_root, species, collection, area=_LEGACY_COLLECTIONS))


def _repair(store):
    """Make a damaged store installable again, if a forced install can; otherwise raise."""
    problem = _store_problem(store)
    if problem is None:
        return
    if not problem[2]:
        raise _damaged(store, problem)
    marker = store.path / _bundle_store().MARKER
    if problem[0].startswith("has lost"):
        # A lost marker is rewritten; the bundle published next replaces a broken newest one.
        temporary = store.path / f".datacache-json-{uuid4().hex}"
        temporary.write_text(json.dumps({"format": _bundle_store().FORMAT, "kind": _STORE_KIND}))
        temporary.replace(marker)


def _publish(path, legacy, build, verify, describe, reuse, *, force, rebuild_damaged=False):
    """Build a bundle with ``build(directory)`` and make it the store's current one.

    ``reuse(current)`` returns the installed bundle to keep, or None to build a
    new one. It is asked once without the lock, so a cache hit writes nothing,
    and again under it, since another install may have finished meanwhile.
    ``force`` builds a new bundle whatever is installed, repairing the store
    where it can; ``rebuild_damaged`` does that only for a damaged store.
    Installs into one store take turns. ``verify`` checks the staged bundle,
    hashing it, and ``describe`` returns the published one without hashing it
    again.
    """
    store = _store(path)
    if not force:
        try:
            current = _current(path, legacy)
        except ValueError:
            if not rebuild_damaged:
                raise
            force = True
        else:
            if current is not None and (found := reuse(current)) is not None:
                return found
    path.parent.mkdir(parents=True, exist_ok=True)
    tools = _bundle_store()
    with store.lock():
        if force:
            _repair(store)
        store.create()
        if not force:
            current = _current(path, legacy)
            if current is not None and (found := reuse(current)) is not None:
                return found
        # Installs hold the lock, so a working directory of this user's found
        # here was left by one that died; it is swept before this one starts.
        prefix = f".staging-{tools.user_key()}-"
        for stale in path.glob(prefix + "*"):
            tools.discard_private_directory(stale)
        working = tools.private_directory(path / f"{prefix}{uuid4().hex}")
        try:
            # The working directory is private; the bundle inside it takes the
            # cache's normal sharing mode, so others can read it once published.
            staged = working / "bundle"
            staged.mkdir()
            build(staged)
            verify(staged)
            published = store.publish(staged)
        finally:
            tools.discard_private_directory(working)
    return describe(published)


def prune_viral_references(
    *,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    keep: int = 1,
    cache_root: str | PathLike[str] | None = None,
) -> tuple[Path, ...]:
    """Delete all but the newest ``keep`` bundles of a collection and of its selections.

    A refresh keeps the bundle it replaces, as DataCache keeps its own,
    because a running program may still be reading it; nothing is pruned
    automatically. Paths into deleted bundles stop working, so prune only
    when no program still uses them. A directory installed by ORFdata 0.44 or
    earlier counts as the oldest bundle. The PyEnsembl indexes
    :func:`viral_genome` built under ``cache_root`` for deleted bundles are
    removed too. Returns the deleted bundle directories.
    """
    if isinstance(keep, bool) or not isinstance(keep, int) or keep < 1:
        raise ValueError("keep must be a positive integer; the current bundle is always kept")
    species, collection = _selection(species, collection)
    stores = [(_bundle_path(cache_root, species, collection),
               _bundle_path(cache_root, species, collection, area=_LEGACY_COLLECTIONS))]
    selections = _bundle_path(cache_root, species, collection, area=_SELECTIONS)
    legacy_selections = _bundle_path(cache_root, species, collection, area=_LEGACY_SELECTIONS)
    names = set()
    for parent in (selections, legacy_selections):
        if _is_directory(parent):
            names.update(entry.name for entry in parent.iterdir() if not entry.name.startswith("."))
    stores.extend((selections / name, legacy_selections / name) for name in sorted(names))
    deleted = []
    for path, legacy in stores:
        store = _store(path)
        problem = _store_problem(store)
        if problem is not None:
            raise _damaged(store, problem)
        if not _is_directory(path):
            continue  # At most a 0.44 directory: the current bundle, which is kept.
        with store.lock():
            bundles = [legacy] if _is_directory(legacy) else []
            if _is_directory(store.bundles):
                bundles += [store.bundles / name for name in store.bundle_names()]
            for bundle in bundles[:-keep]:
                if not _is_directory(bundle):
                    continue  # Never follow a link out of the store.
                index = _annotation_index(cache_root, collection, bundle)
                shutil.rmtree(bundle)
                shutil.rmtree(index, ignore_errors=True)
                deleted.append(bundle)
    return tuple(deleted)


def _annotation_index(cache_root, collection, path):
    """Where :func:`viral_genome` keeps the PyEnsembl index of the files in ``path``."""
    from .datasets import get_data_dir

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
    timeout: float | None = 120,
) -> ViralReferenceFiles:
    """Open, or explicitly download and open, a viral source-reference collection.

    ``download_if_missing=True`` calls :func:`download_viral_references` when
    needed; otherwise this performs an offline, checksummed cache lookup.
    ``force=True`` explicitly refreshes or repairs the cache, which downloads.
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
        raise FileNotFoundError(
            f"The {collection.name} viral references are not installed; run "
            f"'orfdata download viral --collection {collection.name}', or pass "
            "download_if_missing=True."
        )
    return _verified_existing_bundle(path, species=species, collection=collection)


def download_viral_references(
    *,
    species: str | int | None = None,
    collection: ViralCollection | str | None = None,
    cache_root: str | PathLike[str] | None = None,
    force: bool = False,
    timeout: float | None = 120,
    prefer_reference_molecule: bool = True,
    annotate_from_alternate: bool = True,
    prefer_alternate_annotation: bool = False,
) -> ViralReferenceFiles:
    """Install a viral genome, transcript and protein collection for an analysis target.

    :func:`viral_collections` lists what can be installed; the one published
    today is ``gdc-viruses``, the 200-reference GDC screen for human
    cancer-transcriptome analysis. Target species does not assert that every
    source virus has that host, collection membership is not a per-virus
    causality classification, and analysis contexts are not exclusive. To
    create a filtered derived set, use :func:`select_viral_references`. The
    collection is stored under ``ORFDATA_CACHE`` (or the normal OpenVax cache);
    a verified cache is reused and ``force=True`` refreshes or repairs it. A
    refresh is built beside the installed collection and switched in with one
    rename, so a reader finds the old files or the new ones, never a gap; the
    current directory therefore changes with each refresh. Replaced copies
    stay for programs still reading them until :func:`prune_viral_references`
    removes them, and installs of one collection take turns. A request that
    stalls for ``timeout`` seconds is retried, then fails; ``timeout=None``
    waits.

    An inventory can name an entry at two sources whose molecules differ. With
    ``prefer_reference_molecule=True`` (the default) the collection carries the
    one the reference it follows was built from — for GDC's viral decoys the
    alternate source's, which is checkable for HPV16 — and records
    ``genome_source`` on every genome. With
    ``annotate_from_alternate=True`` (the default) an entry whose own record
    annotates nothing takes its coding features from the alternate, when that
    is the identical molecule, and ``prefer_alternate_annotation=True`` widens
    that to every entry the alternate annotates at least as fully, trading
    each submission's own gene naming for one source's uniform models. All three
    verify against packaged checksums before anything is written, and all
    three degrade to the primary record, with the reason recorded, rather than
    discarding a download.

    This source-record collection is separate from the curated
    :func:`viral_references` catalogue.
    """
    species, collection = _selection(species, collection)
    if type(force) is not bool:
        raise ValueError("force must be True or False")
    if timeout is not None and (not isinstance(timeout, (int, float))
                                or isinstance(timeout, bool) or timeout <= 0):
        raise ValueError("timeout must be a positive number of seconds, or None to wait")

    def reuse(current):
        if force:
            return None
        try:
            return _verified_existing_bundle(current, species=species, collection=collection)
        except (OSError, ValueError) as error:
            raise FileNotFoundError(
                "Viral source-reference cache is incomplete or invalid; rerun with force=True."
            ) from error

    def build(directory):
        inventory = viral_panel_coverage(collection=collection.name)
        primaries = _primary_ledger(collection, inventory)
        _write_bundle(directory, collection, inventory, primaries, timeout=timeout,
                      prefer_reference_molecule=prefer_reference_molecule,
                      annotate_from_alternate=annotate_from_alternate,
                      prefer_alternate_annotation=prefer_alternate_annotation)

    return _publish(
        _bundle_path(cache_root, species, collection),
        _bundle_path(cache_root, species, collection, area=_LEGACY_COLLECTIONS), build,
        lambda path: _verified_existing_bundle(path, species=species, collection=collection),
        lambda path: _existing_bundle(path, species=species, collection=collection),
        reuse, force=force)


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
    """The store a selection is published in; ``area`` selects the 0.44 directory."""
    from .datasets import get_data_dir

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
        source_collection=bundle.path,
        **_paths(path),
    )


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
        context_match=context_match)
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
    timeout: float | None = 120,
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

    def reuse(current):
        if force:
            return None
        try:
            return _verified_existing_index(current, **identity)
        except (OSError, ValueError):
            return None  # A derived index is rebuilt from the local bundle, not repaired by hand.

    def build(staging):
        transcripts = [row for row in _read_table(bundle.transcript_table)
                       if row["genome_id"] in selected]
        genome_count = _copy_fasta(bundle.genome_fasta, staging / _GENOME_FASTA, selected)
        transcript_count = _copy_fasta(
            bundle.transcript_fasta, staging / _TRANSCRIPT_FASTA,
            {row["transcript_id"] for row in transcripts})
        protein_count = _copy_fasta(
            bundle.protein_fasta, staging / _PROTEIN_FASTA,
            {row["protein_id"] for row in transcripts})
        genome_rows = _copy_table(bundle.genome_table, staging / _GENOME_TABLE,
                                  "genome_id", selected)
        transcript_rows = _copy_table(bundle.transcript_table, staging / _TRANSCRIPT_TABLE,
                                      "genome_id", selected)
        _copy_gtf(bundle.annotation_gtf, staging / _ANNOTATION_GTF, selected)
        if (genome_count != genome_rows or genome_count != len(selected)
                or transcript_count != transcript_rows or protein_count != transcript_rows):
            raise ValueError("Viral reference index counts disagree")
        (staging / _METADATA).write_text(json.dumps({
            "schema": INDEX_SCHEMA,
            "bundle_version": bundle_version(viral_collection(bundle.collection)),
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
            f"version `{bundle_version(viral_collection(bundle.collection))}`, with the same "
            "identifiers and the matching GTF rows. Context, subset and collection membership "
            "are evidence associations and do not alter identity.\n",
            encoding="utf-8",
        )
        _write_checksums(staging)

    return _publish(
        _index_path(*location), _index_path(*location, area=_LEGACY_SELECTIONS), build,
        lambda path: _verified_existing_index(path, **identity),
        lambda path: _index_result(path, **identity), reuse, force=force, rebuild_damaged=True)


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
