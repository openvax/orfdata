"""Native RibORF reconstruction and conservative, sequence-checked RNA placement."""

import gzip
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import replace

from ..catalog import ORFCatalog, ORFEvidence, ORFRecord, _placement_position
from ..evidence_summary import EvidenceSummary
from ..model import ORFModel, _contig
from ..reference import ScopedTarget
from ..reference_mapping import ReferenceMapping
from ..sequence import genetic_code
from .ncorf_mammals import _ATTRIBUTE, _blocks
from .riborf import PRESENCE, STOPS, _parsed, _rows


def _sha(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def _annotations(path):
    """Keep full transcript versions and every distinct annotated locus."""
    groups = {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = dict(_ATTRIBUTE.findall(fields[8]))
            name = attributes["transcript_id"]
            version = attributes.get("transcript_version")
            if version and "." not in name:
                name += "." + version
            key = name, _contig(fields[0]), fields[6], attributes["gene_id"]
            entry = groups.setdefault(
                key,
                dict(
                    transcript=name,
                    contig=key[1],
                    strand=key[2],
                    gene_id=key[3],
                    gene_name=attributes.get("gene_name", ""),
                    exons=set(),
                ),
            )
            entry["exons"].add((int(fields[3]), int(fields[4])))
    found = defaultdict(list)
    for entry in groups.values():
        entry["exons"] = tuple(sorted(entry["exons"], reverse=entry["strand"] == "-"))
        ascending = sorted(entry["exons"])
        if any(b >= c for (_, b), (c, _) in zip(ascending, ascending[1:])):
            raise ValueError("Overlapping/contradictory reference exons: " + entry["transcript"])
        found[entry["transcript"]].append(entry)
    return found


def _fasta(paths):
    """Full versioned IDs; identical duplicates coalesce, conflicting ones fail."""
    found = {}
    for path in paths:
        name, parts = None, []

        def save():
            if name is not None:
                sequence = "".join(parts).upper()
                if name in found and found[name] != sequence:
                    raise ValueError("Conflicting FASTA records for transcript " + name)
                found[name] = sequence

        with gzip.open(path, "rt") as stream:
            for line in stream:
                if line.startswith(">"):
                    save()
                    name, parts = line[1:].split()[0].split("|")[0], []
                elif name is not None:
                    parts.append(line.strip())
        save()
    return found


def _bed(row):
    if len(row) != 12:
        raise ValueError("Native mouse ORF must provide all twelve BED fields")
    first, last, count = int(row[1]), int(row[2]), int(row[9])
    sizes = [int(x) for x in str(row[10]).rstrip(",").split(",")]
    starts = [int(x) for x in str(row[11]).rstrip(",").split(",")]
    if count != len(sizes) or count != len(starts) or any(size <= 0 for size in sizes):
        raise ValueError("Native BED block counts/sizes disagree")
    blocks = tuple(
        (first + offset + 1, first + offset + size) for offset, size in zip(starts, sizes)
    )
    if (
        not blocks
        or blocks[0][0] != first + 1
        or blocks[-1][1] != last
        or any(b >= c for (_, b), (c, _) in zip(blocks, blocks[1:]))
        or (int(row[6]), int(row[7])) != (first, last)
    ):
        raise ValueError("Native BED span, thick span or block order disagree")
    return _contig(str(row[0])), row[5], blocks


def _coding(sequence, start, end, start_codon):
    if not 1 <= start < end <= len(sequence) + 1:
        raise ValueError("Native RNA interval is outside the full transcript")
    coding = sequence[start - 1 : end - 1]
    if (
        len(coding) < 6
        or len(coding) % 3
        or set(coding) - set("ACGT")
        or coding[:3] != start_codon
        or coding[-3:] not in STOPS
        or any(coding[i : i + 3] in STOPS for i in range(0, len(coding) - 3, 3))
    ):
        raise ValueError("Native RNA does not encode the declared start through its first stop")
    code = genetic_code(1)
    declared = None if start_codon in code.start_codons else start_codon
    return coding, code.translate(coding, start_codon=declared)[:-1], declared


def reconstruct(inputs, manifest):
    """Return checked target records and retained native/unresolved source facts."""
    if (
        manifest.processing.get("importer") != "riborf-native-v3"
        or len(manifest.source_references) != 1
        or len(manifest.references) != 1
    ):
        raise ValueError("Native RibORF placement requires explicit source and target references")
    assets = {asset.asset_id: asset for asset in manifest.sources}
    if set(inputs) != set(assets):
        raise ValueError("Native RibORF placement requires every declared source file")
    for name, path in inputs.items():
        with open(path, "rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != assets[name].sha256:
                raise ValueError("RibORF source checksum mismatch: " + name)
    source_reference, target_reference = manifest.source_references[0], manifest.references[0]
    identities = [
        (r.taxon_id, r.genome_build, r.annotation_name, r.annotation_version)
        for r in (source_reference, target_reference)
    ]
    if identities != [(10090, "GRCm38", "gencode", "M25"), (10090, "GRCm39", "ensembl", "107")]:
        raise ValueError("Native RibORF v3 requires mouse GENCODE M25 and Ensembl 107 reference identities")
    if (
        source_reference.sequence_sha256 != assets["source_rna"].sha256
        or source_reference.source_sha256 != assets["source_gtf"].sha256
    ):
        raise ValueError("Native reference identity disagrees with its RNA/annotation pins")
    source_annotations = _annotations(inputs["source_gtf"])
    source_sequences = _fasta([inputs["source_rna"]])
    target_annotations = _annotations(inputs["gtf"])
    target_sequences = _fasta([inputs["cdna"], inputs["ncrna"]])
    by_sequence = defaultdict(list)
    by_stable_id = defaultdict(list)
    for name, sequence in target_sequences.items():
        by_sequence[_sha(sequence)].append(name)
        by_stable_id[name.rsplit(".", 1)[0]].append(name)
    processing = manifest.processing
    parsed, source_rows, short = _parsed(processing, _rows(inputs["table"], processing["sheet"]))
    records, unresolved, seen = [], [], set()
    source_pins = {name: dict(url=a.url, sha256=a.sha256) for name, a in assets.items()}
    native_valid = 0
    for transcript_label, length, start, end, orf_type, start_codon, row in parsed:
        row_id = row[3]
        if row_id in seen:
            raise ValueError("Duplicate native RibORF row identifier: " + row_id)
        seen.add(row_id)
        transcript = transcript_label.split(":")[0]
        provenance = dict(
            source_row_id=row_id,
            source_row=list(row),
            source_transcript=transcript,
            source_rna_interval=dict(start=start, end=end, coordinates="1-based-half-open"),
            source_pins=source_pins,
        )
        try:
            entries = source_annotations.get(transcript, ())
            if len(entries) != 1:
                raise ValueError(
                    "Native transcript version is absent or has multiple annotated loci"
                )
            entry = entries[0]
            sequence = source_sequences.get(transcript)
            if sequence is None:
                raise ValueError("Native transcript version has no source RNA sequence")
            if len(sequence) != length or sum(b - a + 1 for a, b in entry["exons"]) != length:
                raise ValueError(
                    "Native transcript length disagrees with RNA, annotation or source row"
                )
            blocks = _blocks(entry, start, end - 1)
            if _bed(row) != (entry["contig"], entry["strand"], blocks):
                raise ValueError(
                    "Native RNA offsets and BED blocks disagree with the source annotation"
                )
            labels = transcript_label.split(":")
            if len(labels) != 3 or (_contig(labels[1]), labels[2]) != (
                entry["contig"],
                entry["strand"],
            ):
                raise ValueError("Native transcript label and BED locus disagree")
            coding, protein, declared = _coding(sequence, start, end, start_codon)
        except ValueError as error:
            unresolved.append(dict(provenance, reason=str(error), native_record=None))
            continue
        native_valid += 1
        orf_id = "riborf_" + row_id
        native_model = ORFModel(
            orf_id,
            transcript,
            entry["gene_id"],
            entry["gene_name"],
            entry["contig"],
            entry["strand"],
            blocks,
            protein,
            orf_type,
            source_reference.genome_build,
            source_reference.taxon_id,
            source_reference.species,
            source=processing["source"],
            reference_context=source_reference,
            start_codon=declared,
        )
        source_target = ScopedTarget(source_reference, "orf", orf_id, processing["source"])
        summary = EvidenceSummary(
            processing["source"],
            processing["release"],
            source_target,
            _sha(protein),
            "ribosome",
            "RibORF 2.0 translated ORF",
            1,
            processing.get("presence", PRESENCE),
            "orf",
            assets["table"].url,
            assets["table"].sha256,
        )
        native = ORFRecord.from_model(native_model, ORFEvidence(summaries=(summary,)))
        provenance.update(source_coding_sequence=coding, source_rna_sha256=_sha(sequence))
        candidates, failures = [], []
        # Search all target RNA sequences, not just the surviving stable ID.
        for name in sorted(by_sequence.get(_sha(sequence), ())):
            target_sequence = target_sequences[name]
            if target_sequence != sequence:  # Hash equality alone is not the sequence comparison.
                continue
            if not target_annotations.get(name):
                failures.append(name + ": exact RNA has no annotated exon path")
            for target_entry in target_annotations.get(name, ()):
                if sum(b - a + 1 for a, b in target_entry["exons"]) != len(target_sequence):
                    failures.append(name + ": RNA and annotated exon lengths disagree")
                    continue
                target_blocks = _blocks(target_entry, start, end - 1)
                target_coding, target_protein, _ = _coding(target_sequence, start, end, start_codon)
                if target_coding != coding or target_protein != protein:
                    raise ValueError(
                        "Exact RNA comparison failed its CDS/protein consistency check"
                    )
                model = replace(
                    native_model,
                    transcript_id=name,
                    gene_id=target_entry["gene_id"],
                    gene_name=target_entry["gene_name"],
                    contig=target_entry["contig"],
                    strand=target_entry["strand"],
                    coding_blocks=target_blocks,
                    genome_build=target_reference.genome_build,
                    reference_context=target_reference,
                )
                candidates.append(ORFRecord.from_model(model))
        candidates.sort(key=lambda r: (r.transcript_id, r.contig, r.strand, r.coding_blocks))
        validation = dict(
            frame=dict(translation_table=1, source_start=start, target_start=start),
            sequence=dict(
                source_sha256=_sha(protein),
                target_sha256=_sha(protein),
                all_candidates_checked=True,
            ),
            coding_sequence=dict(
                source_sha256=_sha(coding), target_sha256=_sha(coding), exact=True
            ),
            rna_sequence=dict(
                source_sha256=_sha(sequence),
                target_sha256=_sha(sequence),
                exact=True,
                search_scope="All transcripts in pinned target cDNA and ncRNA FASTAs",
            ),
        )
        provenance["candidate_positions"] = [_placement_position(c) for c in candidates]
        if len(candidates) == 1 and not failures:
            target = candidates[0]
            validation["position"] = dict(
                source=_placement_position(native), target=_placement_position(target)
            )
            mapping = ReferenceMapping(
                orf_id,
                manifest.catalog_id,
                manifest.version,
                source_target,
                target_reference,
                (ScopedTarget(target_reference, "orf", orf_id, processing["source"]),),
                "unique",
                "exact_full_rna_and_cds_with_native_bed_validation",
                validation,
                provenance=provenance,
            )
            records.append(
                replace(
                    target,
                    summaries=native.summaries,
                    lifted_from=native,
                    reference_mapping=mapping,
                )
            )
            continue
        reason = (
            "Multiple target transcript/locus placements have the exact native RNA"
            if len(candidates) > 1
            else "No target transcript preserves the complete native RNA"
            if not candidates and not failures
            else "Exact target RNA has incomplete/contradictory annotation"
        )
        validation["position"] = dict(
            source=_placement_position(native), candidates=provenance["candidate_positions"]
        )
        targets = (
            tuple(
                ScopedTarget(
                    target_reference,
                    "orf",
                    orf_id
                    + "@"
                    + c.transcript_id
                    + "@"
                    + hashlib.sha256(
                        json.dumps(_placement_position(c), sort_keys=True).encode()
                    ).hexdigest()[:16],
                    processing["source"],
                )
                for c in candidates
            )
            if len(candidates) > 1
            else ()
        )
        if not targets:
            # No validated placement is asserted for an absent/changed/incompletely annotated RNA.
            validation["sequence"]["all_candidates_checked"] = False
            validation["rna_sequence"]["exact"] = False
        provenance["annotation_failures"] = failures
        stable_candidates = by_stable_id.get(transcript.rsplit(".", 1)[0], ())
        provenance["stable_id_candidates"] = [
            dict(
                transcript=name,
                rna_sha256=_sha(target_sequences[name]),
                rna_length=len(target_sequences[name]),
                same_rna=target_sequences[name] == sequence,
                same_cds=target_sequences[name][start - 1 : end - 1] == coding,
            )
            for name in stable_candidates
        ]
        mapping = ReferenceMapping(
            orf_id,
            manifest.catalog_id,
            manifest.version,
            source_target,
            target_reference,
            targets,
            "ambiguous" if targets else "unresolved",
            "exact_full_rna_and_cds_with_native_bed_validation",
            validation,
            reason,
            provenance,
        )
        unresolved.append(
            dict(reason=reason, native_record=native.to_dict(), reference_mapping=mapping.to_dict())
        )
    counts = dict(
        records=len(records),
        distinct_protein_sequences=len({r.protein_sequence for r in records}),
        evidence_summaries=len(records),
        source_rows=source_rows,
        short_rows=short,
        native_valid_rows=native_valid,
        resolved_rows=len(records),
        unresolved_rows=len(unresolved),
        short_records=sum(
            len(r.protein_sequence) < processing["report_protein_lengths_below"] for r in records
        ),
    )
    return records, unresolved, counts


def import_orfs(inputs, output, manifest):
    records, unresolved, counts = reconstruct(inputs, manifest)
    if counts != dict(manifest.coverage):
        raise ValueError(
            "Native RibORF coverage differs from manifest: " + json.dumps(counts, sort_keys=True)
        )
    reference = manifest.references[0]
    metadata = dict(
        genome_build=reference.genome_build,
        taxon_id=reference.taxon_id,
        species=reference.species,
        coordinates="1-based-inclusive",
        dataset_version=manifest.processing["release"],
        sources=[
            dict(
                id=manifest.processing["source"],
                release=manifest.processing["release"],
                url=next(a.url for a in manifest.sources if a.asset_id == "table"),
                scope=manifest.processing["scope"],
            )
        ],
        unresolved=unresolved,
        unresolved_reasons=dict(Counter(item["reason"] for item in unresolved)),
        catalog_manifest=manifest.to_dict(),
    )
    catalog = ORFCatalog(records, metadata)
    manifest.validate_catalog(catalog)
    with (
        open(output, "wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream,
    ):
        stream.write(
            json.dumps(
                catalog.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        )
