"""Manage cached source records for ORFdata's GDC viral reference set."""

import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import tempfile
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .reference import canonical_species, species_taxon
from .viral_coverage import PANEL, normalize_viral_contexts, viral_panel_coverage


NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SUMMARY_RESOURCE = "data/viral/gdc-grch38-d1-vd1-genbank-summaries.tsv"
SUMMARY_SHA256 = "5d6c7cbc3f77ba5236acc2e94bc49cf78e3ff486c12535d2b51e80f3ab875498"
BUNDLE_SCHEMA_VERSION = "2"
BUNDLE_VERSION = f"v{BUNDLE_SCHEMA_VERSION}-{SUMMARY_SHA256}"
SUMMARY_HEADER = (
    "index", "abbreviation", "genbank_accession", "accession_version", "title", "length_nt",
    "taxid", "createdate", "updatedate",
)
USER_AGENT = "orfdata-gdc-viral-panel/1.0 (https://github.com/openvax/orfdata)"
_BATCH_SIZE = 25
_MAX_ATTEMPTS = 5
VIRAL_TARGET_SPECIES = "homo_sapiens"
VIRAL_TARGET_TAXON_ID = 9606
VIRAL_COLLECTION = "gdc-grch38-d1-vd1"
VIRAL_CONTEXTS = ("cancer-transcriptome-screening", "cancer-causing")
_COLLECTION_ALIASES = {
    "gdc": VIRAL_COLLECTION,
    "gdc-grch38-d1-vd1": VIRAL_COLLECTION,
    "grch38.d1.vd1": VIRAL_COLLECTION,
}


@dataclass(frozen=True)
class ViralReferenceBundle:
    """Paths and record counts for the installed viral source-reference bundle."""

    target_species: str
    target_taxon_id: int
    collection: str
    contexts: tuple[str, ...]
    panel: str
    path: Path
    genome_fasta: Path
    coding_sequence_fasta: Path
    protein_fasta: Path
    manifest: Path
    cds_table: Path
    genbank_directory: Path
    genome_count: int
    coding_sequence_count: int
    protein_count: int

    def to_dict(self):
        """Return a JSON-ready description of the installed bundle."""
        result = asdict(self)
        for key in ("path", "genome_fasta", "coding_sequence_fasta", "protein_fasta",
                    "manifest", "cds_table", "genbank_directory"):
            result[key] = str(result[key])
        return result


@dataclass(frozen=True)
class ViralReferenceIndex:
    """Cached FASTAs for a many-to-many context selection."""

    target_species: str
    target_taxon_id: int
    collection: str
    contexts: tuple[str, ...]
    excluded_contexts: tuple[str, ...]
    context_match: str
    path: Path
    source_bundle: Path
    genome_fasta: Path
    coding_sequence_fasta: Path
    protein_fasta: Path
    manifest: Path
    cds_table: Path
    genome_count: int
    coding_sequence_count: int
    protein_count: int

    def to_dict(self):
        """Return a JSON-ready description of the context index."""
        result = asdict(self)
        for key in ("path", "source_bundle", "genome_fasta", "coding_sequence_fasta",
                    "protein_fasta", "manifest", "cds_table"):
            result[key] = str(result[key])
        return result


@dataclass(frozen=True)
class _GenBankSummary:
    index: int
    abbreviation: str
    genbank_accession: str
    accession_version: str
    viral_taxon_id: int


def _selection(species, collection):
    """Normalize the one currently published viral collection's selectors."""
    taxon_id = species if type(species) is int else species_taxon(species)
    if taxon_id != VIRAL_TARGET_TAXON_ID:
        raise ValueError(
            "Viral reference downloads are currently available only for target species "
            "human (NCBI taxon 9606)."
        )
    normalized = str(collection).strip().lower().replace("_", "-")
    if normalized not in _COLLECTION_ALIASES:
        raise ValueError(
            f"Unknown viral collection {collection!r}; available: {VIRAL_COLLECTION} ({PANEL})."
        )
    return canonical_species(VIRAL_TARGET_SPECIES), _COLLECTION_ALIASES[normalized]


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _summaries():
    raw = files("orfdata").joinpath(SUMMARY_RESOURCE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != SUMMARY_SHA256:
        raise ValueError("Bundled GDC GenBank summary checksum mismatch")
    lines = raw.decode("utf-8").splitlines()
    if len(lines) < 5 or lines[0] != "# schema: orfdata.gdc_genbank_summaries.v1":
        raise ValueError("Unsupported GDC GenBank summary metadata")
    reader = csv.DictReader(io.StringIO("\n".join(line for line in lines if not line.startswith("#"))),
                            delimiter="\t")
    if tuple(reader.fieldnames or ()) != SUMMARY_HEADER:
        raise ValueError("Unsupported GDC GenBank summary columns")
    summaries = tuple(_GenBankSummary(
        index=int(row["index"]),
        abbreviation=row["abbreviation"],
        genbank_accession=row["genbank_accession"],
        accession_version=row["accession_version"],
        viral_taxon_id=int(row["taxid"]),
    ) for row in reader)
    if len(summaries) != 200 or tuple(summary.index for summary in summaries) != tuple(range(1, 201)):
        raise ValueError("Incomplete GDC GenBank summary")
    return summaries


def _normalise(value):
    """Make one FASTA header field safe and single-token."""
    return re.sub(r"[^A-Za-z0-9._:+-]+", "_", value.strip()).strip("_") or "unknown"


def _bio_seqio():
    try:
        from Bio import SeqIO
    except ImportError as error:
        raise ImportError(
            "Downloading viral source records needs Biopython, which is an ORFdata dependency; "
            "reinstall orfdata."
        ) from error
    return SeqIO


def _fetch_genbank(accessions, *, timeout):
    query = urlencode({
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "gbwithparts",
        "retmode": "text",
        "tool": "orfdata-gdc-viral-panel",
        "email": "support@openvax.org",
    })
    request = Request(f"{NCBI_EFETCH_URL}?{query}", headers={"User-Agent": USER_AGENT})
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        try:
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError) as error:
            if attempt == _MAX_ATTEMPTS:
                raise RuntimeError(f"NCBI EFetch failed for {', '.join(accessions)}") from error
            time.sleep(2 ** (attempt - 1))
    raise AssertionError("unreachable")


def _cds_products(record):
    """Yield source-annotated CDS products in their GenBank feature order."""
    for ordinal, feature in enumerate((feature for feature in record.features
                                       if feature.type == "CDS"), 1):
        qualifiers = feature.qualifiers
        coding_sequence = str(feature.extract(record.seq)).upper()
        codon_start = int(qualifiers.get("codon_start", ["1"])[0])
        translation_table = int(qualifiers.get("transl_table", ["1"])[0])
        translation = qualifiers.get("translation", [None])[0]
        if not translation:
            translated = feature.extract(record.seq)[codon_start - 1:].translate(
                table=translation_table)
            translation = str(translated).rstrip("*")
        yield {
            "ordinal": ordinal,
            "gene": qualifiers.get("gene", [""])[0],
            "product": qualifiers.get("product", [""])[0],
            "protein_id": qualifiers.get("protein_id", [""])[0],
            "location": str(feature.location),
            "coding_sequence": coding_sequence,
            "codon_start": codon_start,
            "translation_table": translation_table,
            "translation": translation.replace("\n", "").replace(" ", ""),
        }


def _write_fasta(stream, header, sequence):
    stream.write(f">{header}\n")
    for offset in range(0, len(sequence), 80):
        stream.write(f"{sequence[offset:offset + 80]}\n")


def _write_checksums(path):
    with (path / "SHA256SUMS").open("w") as stream:
        for file_path in sorted((item for item in path.rglob("*") if item.is_file()
                                 and item.name != "SHA256SUMS"),
                                key=lambda item: item.relative_to(path).as_posix()):
            stream.write(f"{_sha256(file_path)}  {file_path.relative_to(path).as_posix()}\n")


def _quote_empty_final_tsv_fields(path):
    """Keep empty trailing TSV values parseable without trailing whitespace."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    path.write_text("".join(
        f'{line[:-2]}\t""\n' if line.endswith("\t\n") else line
        for line in lines), encoding="utf-8")


def _result(path, *, species, collection):
    manifest_path = path / "manifest.tsv"
    with manifest_path.open(newline="") as stream:
        records = tuple(csv.DictReader(stream, delimiter="\t"))
    if len(records) != 200 or tuple(int(record["index"]) for record in records) != tuple(range(1, 201)):
        raise ValueError(f"Incomplete GDC viral bundle manifest: {manifest_path}")
    protein_count = sum(int(record["cds_count"]) for record in records)
    return ViralReferenceBundle(
        target_species=species,
        target_taxon_id=VIRAL_TARGET_TAXON_ID,
        collection=collection,
        contexts=VIRAL_CONTEXTS,
        panel=PANEL,
        path=path,
        genome_fasta=path / "genomes.fasta",
        coding_sequence_fasta=path / "coding_sequences.fasta",
        protein_fasta=path / "proteins.fasta",
        manifest=manifest_path,
        cds_table=path / "cds.tsv",
        genbank_directory=path / "genbank",
        genome_count=len(records),
        coding_sequence_count=protein_count,
        protein_count=protein_count,
    )


def _existing_bundle(path, *, species, collection):
    """Check the bundle structure without hashing every installed file."""
    result = _result(path, species=species, collection=collection)
    required_files = (
        result.genome_fasta,
        result.coding_sequence_fasta,
        result.protein_fasta,
        result.cds_table,
        path / "SHA256SUMS",
    )
    if any(not file_path.is_file() for file_path in required_files):
        raise ValueError(f"Incomplete viral bundle: {path}")
    if (not result.genbank_directory.is_dir()
            or len(tuple(result.genbank_directory.glob("*.gb.gz"))) != result.genome_count):
        raise ValueError(f"Incomplete viral bundle GenBank records: {result.genbank_directory}")
    return result


def _verify_checksums(path):
    sums_path = path / "SHA256SUMS"
    try:
        sums = tuple(line.rstrip("\n").split("  ", 1) for line in sums_path.open())
    except OSError as error:
        raise ValueError(f"Missing viral bundle checksums: {sums_path}") from error
    expected_files = {item.relative_to(path).as_posix() for item in path.rglob("*")
                      if item.is_file() and item.name != "SHA256SUMS"}
    if (len(sums) != len(expected_files) or any(len(parts) != 2 for parts in sums)
            or {relative for _, relative in sums} != expected_files):
        raise ValueError(f"Invalid viral bundle checksum manifest: {sums_path}")
    if any(_sha256(path / relative) != digest for digest, relative in sums):
        raise ValueError(f"Viral bundle checksum mismatch: {path}")


def _verified_existing_bundle(path, *, species, collection):
    result = _existing_bundle(path, species=species, collection=collection)
    _verify_checksums(path)
    for file_path, expected_count in ((result.genome_fasta, result.genome_count),
                                      (result.coding_sequence_fasta,
                                       result.coding_sequence_count),
                                      (result.protein_fasta, result.protein_count)):
        with file_path.open() as stream:
            actual_count = sum(line.startswith(">") for line in stream)
        if actual_count != expected_count:
            raise ValueError(f"FASTA record count does not match manifest: {file_path}")
    with result.cds_table.open(newline="") as stream:
        if sum(1 for _ in csv.DictReader(stream, delimiter="\t")) != result.protein_count:
            raise ValueError(f"CDS row count does not match manifest: {result.cds_table}")
    if len(tuple(result.genbank_directory.glob("*.gb.gz"))) != result.genome_count:
        raise ValueError(f"GenBank record count does not match manifest: {result.genbank_directory}")
    return result


def _write_bundle(path, entries, summaries, *, timeout):
    seqio = _bio_seqio()
    records_dir = path / "genbank"
    records_dir.mkdir()
    expected_versions = {summary.accession_version for summary in summaries}
    downloaded = {}
    for start in range(0, len(summaries), _BATCH_SIZE):
        batch = summaries[start:start + _BATCH_SIZE]
        payload = _fetch_genbank([summary.accession_version for summary in batch], timeout=timeout)
        records = tuple(seqio.parse(io.StringIO(payload.decode("utf-8")), "genbank"))
        if len(records) != len(batch):
            raise ValueError(f"NCBI returned {len(records)} records for {len(batch)} accessions")
        for record in records:
            if record.id not in expected_versions:
                raise ValueError(f"Unexpected NCBI record {record.id}")
            if record.id in downloaded:
                raise ValueError(f"Duplicate NCBI record {record.id}")
            downloaded[record.id] = record
        if start + _BATCH_SIZE < len(summaries):
            time.sleep(0.34)  # Respect NCBI's unauthenticated request rate.
    if set(downloaded) != expected_versions:
        missing = sorted(expected_versions - set(downloaded))
        extra = sorted(set(downloaded) - expected_versions)
        raise ValueError(f"NCBI response mismatch; missing={missing}, extra={extra}")

    manifest_columns = (
        "index", "abbreviation", "name", "source_genbank", "source_refseq",
        "accession_version", "viral_taxon_id", "source_host_taxon_ids", "length_nt",
        "genome_sha256", "cds_count", "pave_id", "pave_notes",
    )
    cds_columns = (
        "abbreviation", "accession_version", "cds_ordinal", "gene", "product", "protein_id",
        "location", "nucleotides", "codon_start", "translation_table", "amino_acids",
    )
    total_cds = 0
    with (path / "genomes.fasta").open("w") as genomes, \
         (path / "coding_sequences.fasta").open("w") as coding_sequences, \
         (path / "proteins.fasta").open("w") as proteins, \
         (path / "manifest.tsv").open("w", newline="") as manifest_stream, \
         (path / "cds.tsv").open("w", newline="") as cds_stream:
        manifest = csv.DictWriter(manifest_stream, fieldnames=manifest_columns, delimiter="\t",
                                  lineterminator="\n")
        cds_writer = csv.DictWriter(cds_stream, fieldnames=cds_columns, delimiter="\t",
                                    lineterminator="\n")
        manifest.writeheader()
        cds_writer.writeheader()
        for entry, summary in zip(entries, summaries, strict=True):
            record = downloaded[summary.accession_version]
            serialized = io.StringIO()
            seqio.write(record, serialized, "genbank")
            record_path = records_dir / f"{_normalise(entry.abbreviation)}__{record.id}.gb.gz"
            with gzip.GzipFile(record_path, "wb", mtime=0) as stream:
                stream.write(serialized.getvalue().encode("utf-8"))
            sequence = str(record.seq).upper()
            genome_id = f"viral:{record.id}"
            _write_fasta(genomes, genome_id, sequence)
            products = tuple(_cds_products(record))
            for product in products:
                protein_id = _normalise(product["protein_id"] or f"CDS{product['ordinal']}")
                header = (
                    f"{genome_id}|cds={product['ordinal']}|protein_id={protein_id}"
                    f"|gene={_normalise(product['gene'])}|product={_normalise(product['product'])}"
                )
                _write_fasta(coding_sequences, header, product["coding_sequence"])
                _write_fasta(proteins, header, product["translation"])
                cds_writer.writerow({
                    "abbreviation": entry.abbreviation,
                    "accession_version": record.id,
                    "cds_ordinal": product["ordinal"],
                    "gene": product["gene"],
                    "product": product["product"],
                    "protein_id": product["protein_id"],
                    "location": product["location"],
                    "nucleotides": len(product["coding_sequence"]),
                    "codon_start": product["codon_start"],
                    "translation_table": product["translation_table"],
                    "amino_acids": len(product["translation"]),
                })
            manifest.writerow({
                "index": entry.index,
                "abbreviation": entry.abbreviation,
                "name": entry.name,
                "source_genbank": entry.genbank,
                "source_refseq": entry.refseq,
                "accession_version": record.id,
                "viral_taxon_id": summary.viral_taxon_id,
                "source_host_taxon_ids": "",
                "length_nt": len(sequence),
                "genome_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "cds_count": len(products),
                "pave_id": entry.pave_id or "",
                "pave_notes": entry.pave_notes or "",
            })
            total_cds += len(products)

    _quote_empty_final_tsv_fields(path / "manifest.tsv")
    retrieved = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    (path / "README.md").write_text(
        f"# {PANEL} viral reference bundle\n\n"
        f"Retrieved from NCBI GenBank at {retrieved} using the exact GenBank accession in "
        "each row of the GDC GRCh38.d1.vd1 viral inventory. The accompanying manifest "
        "retains any RefSeq ID listed by GDC.\n\n"
        "Files:\n\n"
        "- `genomes.fasta`: one nucleotide reference per panel row.\n"
        "- `coding_sequences.fasta`: nucleotide sequences extracted from source CDS features.\n"
        "- `proteins.fasta`: source translations of those CDS features.\n"
        "- `cds.tsv`: CDS feature metadata, source locations and sequence lengths.\n"
        "- `manifest.tsv`: source identifiers, resolved accession versions, sequence hashes and CDS counts.\n"
        "- `genbank/`: gzipped GenBank records, one per source reference.\n\n"
        f"This bundle contains 200 genomes and {total_cds} paired coding and protein sequences. "
        "Genome FASTA IDs are `viral:<accession.version>`; paired CDS/protein IDs add "
        "the source CDS ordinal and labels. Collection and context provenance do not alter IDs.\n\n"
        "The collection is admitted for human analysis. That target does not assert a source host; "
        "`source_host_taxon_ids` is empty when the source inventory does not establish one.\n\n"
        "This is a source-record download, not an ORFdata curated viral-reference release. "
        "It does not alter `viral_references()` or `Proteome.build(..., viral=...)`.\n\n"
        "Source inventory: https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt\n"
        "NCBI retrieval: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\n",
        encoding="utf-8",
    )
    _write_checksums(path)


def _bundle_path(cache_root, species, collection):
    from .datasets import get_data_dir

    return (get_data_dir(cache_root) / "viral" / "references" / species / collection /
            BUNDLE_VERSION)


def viral_reference_bundle(*, species="human", collection=VIRAL_COLLECTION,
                           download_if_missing=False,
                           cache_root=None, force=False, timeout=120):
    """Open, or explicitly download and open, a viral source-reference bundle.

    ``download_if_missing=True`` calls :func:`download_viral_references` when
    needed; otherwise this performs an offline, checksummed cache lookup.
    ``force=True`` explicitly refreshes or repairs the cache. ``species``
    describes the target analysis species; ``collection`` selects a versioned reference
    set. Analysis contexts are nonexclusive metadata and never determine cache
    identity. Only the human GDC collection is currently published.
    """
    if type(download_if_missing) is not bool:
        raise ValueError("download_if_missing must be True or False")
    species, collection = _selection(species, collection)
    if download_if_missing or force:
        return download_viral_references(
            species=species, collection=collection, cache_root=cache_root,
            force=force, timeout=timeout)
    path = _bundle_path(cache_root, species, collection)
    if not path.is_dir():
        raise FileNotFoundError(
            "Viral source references are not installed; run "
            "'orfdata download viral --species human', or pass "
            "download_if_missing=True."
        )
    return _verified_existing_bundle(path, species=species, collection=collection)


def download_viral_references(*, species="human", collection=VIRAL_COLLECTION, cache_root=None,
                              force=False, timeout=120):
    """Install a viral genome and annotated-CDS collection for an analysis target.

    The currently published selection is the 200-reference GDC
    cancer-transcriptome-screening collection for a human analysis target.
    Collection membership is not a per-virus cancer-causality classification,
    and analysis contexts are not exclusive. To create a filtered derived
    index, use :func:`index_viral_references`. The bundle is stored under
    ``ORFDATA_CACHE`` (or the normal OpenVax cache); a verified cache is reused
    and ``force=True`` refreshes or repairs it.

    This source-record bundle is separate from the curated
    :func:`viral_references` catalogue.
    """
    species, collection = _selection(species, collection)
    if type(force) is not bool:
        raise ValueError("force must be True or False")
    if not isinstance(timeout, (int, float)) or isinstance(timeout, bool) or timeout <= 0:
        raise ValueError("timeout must be a positive number of seconds")
    output = _bundle_path(cache_root, species, collection)
    if output.exists() and not output.is_dir():
        raise ValueError(f"Viral reference cache path must be a directory: {output}")
    if output.is_dir():
        try:
            existing = _verified_existing_bundle(output, species=species, collection=collection)
        except (OSError, ValueError) as error:
            if not force:
                raise FileNotFoundError(
                    "Viral source-reference cache is incomplete or invalid; rerun with force=True."
                ) from error
        else:
            if not force:
                return existing

    entries = viral_panel_coverage()
    summaries = _summaries()
    if any(entry.genbank != summary.genbank_accession or entry.abbreviation != summary.abbreviation
           for entry, summary in zip(entries, summaries, strict=True)):
        raise ValueError("GDC coverage ledger and GenBank summary do not agree")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        _write_bundle(staging, entries, summaries, timeout=timeout)
        _verified_existing_bundle(staging, species=species, collection=collection)
        if output.exists():
            shutil.rmtree(output)
        staging.replace(output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return _verified_existing_bundle(output, species=species, collection=collection)


def _copy_selected_fasta(source, destination, accessions):
    keep = False
    count = 0
    with source.open() as input_stream, destination.open("w") as output_stream:
        for line in input_stream:
            if line.startswith(">"):
                token = line[1:].split(None, 1)[0]
                reference_id = token.split("|", 1)[0]
                accession = reference_id.removeprefix("viral:")
                keep = accession in accessions
                count += keep
            if keep:
                output_stream.write(line)
    return count


def _copy_selected_tsv(source, destination, accessions):
    with source.open(newline="") as input_stream, destination.open("w", newline="") as output_stream:
        reader = csv.DictReader(input_stream, delimiter="\t")
        writer = csv.DictWriter(output_stream, fieldnames=reader.fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        rows = 0
        for row in reader:
            if row["accession_version"] in accessions:
                writer.writerow(row)
                rows += 1
    return rows


def _index_path(cache_root, species, collection, contexts, excluded_contexts,
                context_match, accessions):
    from .datasets import get_data_dir

    identity = json.dumps({
        "schema": "orfdata.viral_reference_index.v1",
        "bundle_version": BUNDLE_VERSION,
        "contexts": contexts,
        "excluded_contexts": excluded_contexts,
        "context_match": context_match,
        "accessions": accessions,
    }, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    return (get_data_dir(cache_root) / "viral" / "indexes" / species / collection /
            BUNDLE_VERSION / digest)


def _index_result(path, *, bundle, contexts, excluded_contexts, context_match):
    try:
        metadata = json.loads((path / "index.json").read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid viral reference index metadata: {path}") from error
    expected = {
        "schema": "orfdata.viral_reference_index.v1",
        "bundle_version": BUNDLE_VERSION,
        "target_species": bundle.target_species,
        "target_taxon_id": bundle.target_taxon_id,
        "collection": bundle.collection,
        "contexts": list(contexts),
        "excluded_contexts": list(excluded_contexts),
        "context_match": context_match,
    }
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError(f"Viral reference index identity mismatch: {path}")
    return ViralReferenceIndex(
        target_species=bundle.target_species,
        target_taxon_id=bundle.target_taxon_id,
        collection=bundle.collection,
        contexts=contexts,
        excluded_contexts=excluded_contexts,
        context_match=context_match,
        path=path,
        source_bundle=bundle.path,
        genome_fasta=path / "genomes.fasta",
        coding_sequence_fasta=path / "coding_sequences.fasta",
        protein_fasta=path / "proteins.fasta",
        manifest=path / "manifest.tsv",
        cds_table=path / "cds.tsv",
        genome_count=int(metadata["genome_count"]),
        coding_sequence_count=int(metadata["coding_sequence_count"]),
        protein_count=int(metadata["protein_count"]),
    )


def _verified_existing_index(path, *, bundle, contexts, excluded_contexts, context_match):
    result = _index_result(
        path, bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
        context_match=context_match)
    required = (result.genome_fasta, result.coding_sequence_fasta, result.protein_fasta,
                result.manifest, result.cds_table, path / "SHA256SUMS")
    if any(not item.is_file() for item in required):
        raise ValueError(f"Incomplete viral reference index: {path}")
    _verify_checksums(path)
    for fasta, expected in ((result.genome_fasta, result.genome_count),
                            (result.coding_sequence_fasta, result.coding_sequence_count),
                            (result.protein_fasta, result.protein_count)):
        with fasta.open() as stream:
            if sum(line.startswith(">") for line in stream) != expected:
                raise ValueError(f"FASTA record count does not match viral index: {fasta}")
    return result


def index_viral_references(*, contexts, exclude_contexts=None, context_match="any", species="human",
                           collection=VIRAL_COLLECTION, download_if_missing=False,
                           cache_root=None, force=False, timeout=120):
    """Build reusable FASTAs for selected viral contexts minus excluded contexts.

    Context membership is an association and never changes a FASTA record's
    accession-based identity. Repeated ``contexts`` select their union by
    default or their intersection with ``context_match='all'``; excluded
    contexts are then removed. The source GenBank bundle is cached once;
    indexes reuse it rather than downloading the same reference for each
    context.
    """
    contexts = normalize_viral_contexts(contexts)
    excluded_contexts = normalize_viral_contexts(exclude_contexts)
    if not contexts:
        raise ValueError("At least one viral context is required")
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    bundle = viral_reference_bundle(
        species=species, collection=collection, download_if_missing=download_if_missing,
        cache_root=cache_root, timeout=timeout)
    selected_entries = viral_panel_coverage(
        contexts=contexts, exclude_contexts=excluded_contexts, context_match=context_match)
    all_entries = viral_panel_coverage()
    summaries = _summaries()
    versions = {entry.index: summary.accession_version
                for entry, summary in zip(all_entries, summaries, strict=True)}
    accessions = tuple(versions[entry.index] for entry in selected_entries)
    if not accessions:
        raise ValueError("No viral references match the selected contexts")
    selected = set(accessions)
    output = _index_path(
        cache_root, bundle.target_species, bundle.collection, contexts, excluded_contexts,
        context_match, accessions)
    if output.is_dir() and not force:
        return _verified_existing_index(
            output, bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
            context_match=context_match)
    if output.exists() and not output.is_dir():
        raise ValueError(f"Viral reference index path must be a directory: {output}")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        genome_count = _copy_selected_fasta(
            bundle.genome_fasta, staging / "genomes.fasta", selected)
        coding_count = _copy_selected_fasta(
            bundle.coding_sequence_fasta, staging / "coding_sequences.fasta", selected)
        protein_count = _copy_selected_fasta(
            bundle.protein_fasta, staging / "proteins.fasta", selected)
        manifest_count = _copy_selected_tsv(
            bundle.manifest, staging / "manifest.tsv", selected)
        cds_count = _copy_selected_tsv(bundle.cds_table, staging / "cds.tsv", selected)
        if (genome_count != manifest_count or coding_count != protein_count
                or protein_count != cds_count):
            raise ValueError("Viral context index counts disagree")
        metadata = {
            "schema": "orfdata.viral_reference_index.v1",
            "bundle_version": BUNDLE_VERSION,
            "target_species": bundle.target_species,
            "target_taxon_id": bundle.target_taxon_id,
            "collection": bundle.collection,
            "contexts": list(contexts),
            "excluded_contexts": list(excluded_contexts),
            "context_match": context_match,
            "accession_versions": list(accessions),
            "genome_count": genome_count,
            "coding_sequence_count": coding_count,
            "protein_count": protein_count,
        }
        (staging / "index.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (staging / "README.md").write_text(
            "# ORFdata viral reference index\n\n"
            f"Contexts ({context_match}): {', '.join(contexts)}\n\n"
            f"Excluded contexts: {', '.join(excluded_contexts) or 'none'}\n\n"
            "The paired coding/protein FASTAs use accession-based identifiers. "
            "Context and collection membership are provenance associations and do not alter identity.\n",
            encoding="utf-8",
        )
        _write_checksums(staging)
        _verified_existing_index(
            staging, bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
            context_match=context_match)
        if output.exists():
            shutil.rmtree(output)
        staging.replace(output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return _verified_existing_index(
        output, bundle=bundle, contexts=contexts, excluded_contexts=excluded_contexts,
        context_match=context_match)
