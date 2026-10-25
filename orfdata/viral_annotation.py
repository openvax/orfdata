"""Turn source GenBank records into a viral genome, transcriptome and proteome.

One source-annotated CDS feature becomes one coding transcript: its spliced
nucleotide sequence, the source translation of it, and the GTF rows placing it
on the viral genome. Identifiers are derived here once so the three layers and
the annotation always agree, which is what lets PyEnsembl — and through it
varcode, isovar, topiary and vaxrank — read the result.

Flat files are read by :mod:`orfdata.genbank`, which keeps a location's
structure (``complement(join(...))``) and its partial-end markers, so a spliced
CDS stays in transcript order and a truncated one says so.
"""

import re
from itertools import product

from . import genbank
from .gtf import gene_line, transcript_blocks, transcript_lines
from .sequence import genetic_code

COMPLEMENT = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")
# What each IUPAC nucleotide code stands for, so an ambiguous codon can be read.
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT",
         "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT"}
GENBANK_SOURCE = "genbank"
# A gene split by the origin of a circular genome is still one gene; joining the
# two halves is ORFdata's inference, so it is named as a separate provenance.
ORIGIN_JOIN_SOURCE = "genbank:origin_join"
# Qualifiers whose continuation lines are prose, not sequence.
_TEXT_QUALIFIERS = ("gene", "product", "note")


class GenBankRecord:
    """One parsed GenBank record: its sequence, its CDS features and its own text."""

    __slots__ = ("accession_version", "sequence", "source_text", "circular", "coding_features",
                 "definition", "taxon_id")

    def __init__(self, *, accession_version, sequence, source_text, coding_features,
                 circular=False, definition="", taxon_id=0):
        self.accession_version = accession_version
        self.sequence = sequence
        self.source_text = source_text
        self.circular = circular
        self.coding_features = tuple(coding_features)
        self.definition = definition
        self.taxon_id = taxon_id


def safe_token(value):
    """Make one FASTA header or GTF attribute field safe and single-token."""
    return re.sub(r"[^A-Za-z0-9._:+-]+", "_", value.strip()).strip("_") or "unknown"


def reverse_complement(sequence):
    """The reverse complement of a nucleotide string, ambiguity codes included."""
    return sequence.translate(COMPLEMENT)[::-1]


def split_genbank(payload):
    """Split one EFetch flat-file response into the source text of each record.

    Each record comes back exactly as EFetch returns it on its own, including
    the blank line it writes after the ``//`` terminator, so a stored record
    can be compared with the source byte for byte.
    """
    records, current = [], []
    for line in payload.decode("utf-8").splitlines(keepends=True):
        if not current and not line.strip():
            if records:  # The blank line separating two records ends the first one.
                records[-1] += line
            continue
        current.append(line)
        if line.rstrip("\r\n") == "//":
            records.append("".join(current))
            current = []
    if "".join(current).strip():
        raise ValueError("Truncated GenBank response: the last record has no '//' terminator")
    return tuple(records)


def exon_blocks(location):
    """Flatten a source location into ``(start, end, strand)`` blocks in transcript order.

    Blocks are 1-based and inclusive. A ``complement`` wrapper reverses the
    order of what it contains, because the transcript reads the other way.
    """
    kind = type(location).__name__
    if kind == "Range":  # Already 1-based and inclusive, as GenBank writes it.
        return ((location.start, location.end, "+"),)
    if kind == "Complement":
        return tuple((start, end, "-" if strand == "+" else "+")
                     for start, end, strand in reversed(exon_blocks(location.location)))
    if kind in ("Join", "Order"):
        return tuple(block for part in location.locations for block in exon_blocks(part))
    raise ValueError(f"Unsupported CDS location: {location!r}")


def partial_ends(location):
    """Whether the source marks the transcript's 5' and 3' ends as running past the record.

    GenBank writes ``<`` and ``>`` against genomic coordinates, so a
    ``complement`` wrapper swaps which end of the transcript each one describes.
    """
    kind = type(location).__name__
    if kind == "Range":
        return bool(location.before), bool(location.after)
    if kind == "Complement":
        start, end = partial_ends(location.location)
        return end, start
    if kind in ("Join", "Order"):
        parts = location.locations
        return partial_ends(parts[0])[0], partial_ends(parts[-1])[1]
    return False, False


def extract(sequence, exons):
    """The spliced nucleotides of one exon path, in transcript order."""
    return "".join(sequence[start - 1:end] if strand == "+"
                   else reverse_complement(sequence[start - 1:end])
                   for start, end, strand in exons)


def location_text(location):
    """The source location, written back as GenBank writes it.

    A reconstruction from coordinates alone would drop the ``<`` and ``>``
    markers that say a feature runs past the record, so this renders the parsed
    location itself; an exon path from another source renders without markers.
    """
    kind = type(location).__name__
    if kind == "Range":
        if location.start == location.end and not (location.before and location.after):
            return f"{'<' if location.before else '>' if location.after else ''}{location.start}"
        return (f"{'<' if location.before else ''}{location.start}.."
                f"{'>' if location.after else ''}{location.end}")
    if kind == "Complement":
        return f"complement({location_text(location.location)})"
    if kind in ("Join", "Order"):
        inner = ",".join(location_text(part) for part in location.locations)
        return f"{kind.lower()}({inner})"
    if isinstance(location, tuple):  # An exon path, as PaVE and other sources give one.
        inner = ",".join(f"{start}..{end}" for start, end, _ in location)
        text = f"join({inner})" if len(location) > 1 else inner
        return f"complement({text})" if location[0][2] == "-" else text
    raise ValueError(f"Unsupported CDS location: {location!r}")


def translate(nucleotides, *, table):
    """Translate whole codons with ORFdata's genetic codes, dropping a final stop.

    This is a plain codon translation, for reading and checking a source's own
    annotation rather than modelling an ORF: it does not require an initiation
    codon and does not read the first codon as methionine.
    :meth:`GeneticCode.translate` remains the initiation-aware translation.

    A codon carrying an IUPAC ambiguity code still translates to a residue
    whenever every base it stands for gives the same one — ``ATY`` is
    isoleucine either way — and only to ``X`` when they genuinely differ.
    Sequence databases publish the determined residue, so reading such a codon
    as ``X`` would make a correct record look like a disagreement.
    """
    codons = genetic_code(table).codons
    protein = "".join(_residue(codons, nucleotides[offset:offset + 3])
                      for offset in range(0, len(nucleotides) - 2, 3))
    return protein.rstrip("*")


def _residue(codons, codon):
    """One codon's residue, resolving an ambiguity code it does not depend on."""
    plain = codons.get(codon)
    if plain is not None:
        return plain
    try:
        resolved = {codons[option] for option in _resolutions(codon)}
    except KeyError:  # A character that is not an IUPAC nucleotide code at all.
        return "X"
    return resolved.pop() if len(resolved) == 1 else "X"


def _resolutions(codon):
    """Every unambiguous codon an IUPAC-coded one stands for."""
    return ("".join(bases) for bases in product(*(IUPAC[base] for base in codon.upper())))


def agrees_with(translated, published):
    """Whether a translation matches a published one, allowing the initiator convention.

    A source writes the first residue of a coding sequence as methionine
    whatever its initiation codon is, so a leading difference there is the
    convention rather than a disagreement.
    """
    published = (published or "").upper().rstrip("*")
    if not published:
        return False
    return translated == published or (published[0] == "M"
                                       and translated[1:] == published[1:])


# How a coding sequence begins, as the source records it.
ATG_START = "atg"                  # The canonical initiation codon.
ALTERNATIVE_START = "alternative"  # Another initiation codon of the same genetic code.
DECLARED_START = "declared"        # Neither, but the source declares the initiator residue.
NON_CANONICAL_START = "non_canonical"  # Not an initiation codon, and nothing declares it.
PARTIAL_START = "partial"          # The 5' end runs past the record, so there is no start.
UNKNOWN_START = "unknown"          # The source names a genetic code ORFdata does not carry.
INITIATIONS = (ATG_START, ALTERNATIVE_START, DECLARED_START, NON_CANONICAL_START, PARTIAL_START,
               UNKNOWN_START)
# Qualifiers naming a translational recoding event, which a plain codon reading
# of the genome cannot reproduce on its own.
_RECODING_QUALIFIERS = ("ribosomal_slippage", "artificial_location")


def initiation(start_codon, *, table, partial_start=False, declared=False):
    """How a coding sequence begins: which of :data:`INITIATIONS` describes it."""
    if partial_start or len(start_codon) < 3:
        return PARTIAL_START
    if start_codon == "ATG":
        return ATG_START
    if declared:
        return DECLARED_START
    try:
        starts = genetic_code(table).start_codons
    except ValueError:
        # Which codons initiate is a property of the code; without the code
        # there is no answer to give, and inventing one would be a guess.
        return UNKNOWN_START
    return ALTERNATIVE_START if start_codon in starts else NON_CANONICAL_START


def declares_initiator(value, exons, codon_start):
    """Whether a ``transl_except`` qualifier declares the residue of the first codon.

    A source writes a non-canonical initiator as an exception over the codon's
    own coordinates — CMV's UL30A reads ``(pos:complement(38131..38133),aa:Met)``
    — and uses the same qualifier for read-through residues elsewhere in the
    sequence, so only an exception covering the first codon is an initiator.
    The residue itself needs no separate record: it is the first residue of the
    translation the source publishes.
    """
    match = re.fullmatch(r"\(pos:(.+),aa:([A-Za-z]{3})\)", " ".join((value or "").split()))
    if match is None:
        return False
    try:
        excepted = exon_blocks(genbank.parse_location(match.group(1)))
    except ValueError:
        return False
    return _positions(excepted) == _positions(
        transcript_blocks(exons, codon_start - 1, 3))


def _positions(blocks):
    """Every genomic position a set of blocks covers."""
    return {position for block in blocks for position in range(block[0], block[1] + 1)}


def coding_feature(ordinal, exons, nucleotides, qualifiers, *, partial=(False, False),
                   location=None, annotation_source=GENBANK_SOURCE):
    """One CDS as the three layers record it, whichever source annotated it."""
    codon_start = int(qualifiers.get("codon_start") or 1)
    if not 1 <= codon_start <= 3:
        raise ValueError(f"Unsupported codon_start {codon_start}")
    table = int(qualifiers.get("transl_table") or 1)
    translation = qualifiers.get("translation")
    if translation:
        translation = translation.replace("\n", "").replace(" ", "").upper()
    else:
        translation = translate(nucleotides[codon_start - 1:], table=table)
    coding_length = len(nucleotides) - (codon_start - 1)
    start_partial, stop_partial = partial
    start_codon = nucleotides[codon_start - 1:codon_start + 2].upper()
    return {
        "ordinal": ordinal,
        "gene": qualifiers.get("gene") or "",
        # Some records describe the product only in a note, as HBV's polymerase does.
        "product": qualifiers.get("product") or qualifiers.get("note") or "",
        # What to call the gene: never the note, which is prose rather than a label.
        "gene_label": qualifiers.get("gene") or qualifiers.get("product") or "",
        "source_protein_id": qualifiers.get("protein_id") or "",
        "location": location_text(exons if location is None else location),
        "partial": ("both" if start_partial and stop_partial else
                    "5_prime" if start_partial else "3_prime" if stop_partial else ""),
        "exons": exons,
        "nucleotides": nucleotides,
        "codon_start": codon_start,
        # Where translation starts, and what the ribosome had to do that reading
        # the genome codon by codon does not: both are the source's own claim.
        "start_codon": start_codon,
        "initiation": initiation(
            start_codon, table=table, partial_start=start_partial,
            declared=declares_initiator(qualifiers.get("transl_except"), exons, codon_start)),
        "recoding": ",".join(name for name in _RECODING_QUALIFIERS if name in qualifiers),
        "translation_table": table,
        "translation": translation,
        "annotation_source": annotation_source,
        # A complete CDS spans whole codons and ends with the stop codon the
        # translation leaves out, and the source marks neither end as running
        # past the record; a partial or frameshifted one does not.
        "complete": (not start_partial and not stop_partial
                     and coding_length == 3 * (len(translation) + 1)),
    }


def _qualifiers(feature):
    """The first value of each qualifier, with a valueless one kept as an empty string.

    Continuation lines of a prose qualifier join with a space.
    """
    values = {}
    for key, value in feature.qualifiers:
        if key in values:
            continue
        if value is None:  # A qualifier that is a flag; its presence is the claim.
            values[key] = ""
        else:
            values[key] = " ".join(value.split()) if key in _TEXT_QUALIFIERS else value
    return values


def coding_features(record):
    """Source-annotated CDS features of one parsed record, in file order."""
    sequence = record.sequence
    features = []
    for ordinal, feature in enumerate((item for item in record.features if item.kind == "CDS"), 1):
        exons = exon_blocks(feature.location)
        nucleotides = extract(sequence, exons)
        if sum(end - start + 1 for start, end, _ in exons) != len(nucleotides):
            raise ValueError(
                f"CDS {ordinal} of {record.accession_version} runs past the end of the record")
        features.append(coding_feature(
            ordinal, exons, nucleotides, _qualifiers(feature),
            partial=partial_ends(feature.location), location=feature.location))
    return tuple(features)


def join_across_origin(features, sequence):
    """Join a CDS the circular origin splits into a 3'- and a 5'-partial half.

    A gene that crosses the origin of a circular genome cannot be written as one
    span in a linear record, so GenBank writes two partial features: one running
    off the end, one resuming at base 1. They are only joined here when the join
    reads as one gene — whole codons, a start codon, a terminal stop, no stop in
    between, and the two published translations end to end — so a pair that does
    not clearly belong together is left exactly as the source wrote it.
    """
    tail = [f for f in features if f["partial"] == "3_prime"
            and f["exons"][-1][1] == len(sequence) and f["exons"][0][2] == "+"]
    head = [f for f in features if f["partial"] == "5_prime"
            and f["exons"][0][0] == 1 and f["exons"][0][2] == "+"]
    if len(tail) != 1 or len(head) != 1:
        return features
    first, second = tail[0], head[0]
    if first["translation_table"] != second["translation_table"]:
        return features
    exons = first["exons"] + second["exons"]
    nucleotides = extract(sequence, exons)
    table = first["translation_table"]
    try:
        protein = translate(nucleotides, table=table)
    except ValueError:
        return features
    if (len(nucleotides) % 3 or "*" in protein
            or not agrees_with(protein, first["translation"] + second["translation"])
            or genetic_code(table).codons.get(nucleotides[-3:]) != "*"
            or nucleotides[:3] not in genetic_code(table).start_codons):
        return features
    # Each half's description covers only that half — "(AA 1-303)" and the rest —
    # so the joined gene keeps both rather than inheriting one as if it were whole.
    descriptions = [text for text in dict.fromkeys((first["product"], second["product"])) if text]
    joined = coding_feature(
        first["ordinal"], exons, nucleotides,
        {"gene": first["gene"] or second["gene"],
         "note": "; ".join(descriptions),
         "protein_id": first["source_protein_id"] or second["source_protein_id"],
         "translation": protein, "transl_table": table},
        annotation_source=ORIGIN_JOIN_SOURCE)
    kept = [joined if f is first else f for f in features if f is not second]
    return tuple(dict(feature, ordinal=ordinal)
                 for ordinal, feature in enumerate(kept, 1))


def source_taxon(record):
    """The NCBI taxon the record's own ``source`` feature cross-references, if any."""
    for feature in record.features:
        if feature.kind != "source":
            continue
        for name, value in feature.qualifiers:
            if name == "db_xref" and str(value).startswith("taxon:"):
                return int(str(value).split(":", 1)[1])
    return 0


def parse_genbank(payload):
    """Parse one EFetch response into records that keep their own source text."""
    records = []
    for text in split_genbank(payload):
        parsed = genbank.parse(text)
        features = coding_features(parsed)
        if parsed.circular:
            features = join_across_origin(features, parsed.sequence)
        records.append(GenBankRecord(
            accession_version=parsed.accession_version,
            sequence=parsed.sequence,
            source_text=text,
            circular=parsed.circular,
            coding_features=features,
            definition=parsed.definition,
            taxon_id=source_taxon(parsed),
        ))
    return tuple(records)


def name_transcripts(genome_id, features, *, unique_protein_ids):
    """Name each CDS feature's transcript, gene and protein once, deterministically."""
    rows = []
    for feature in features:
        transcript_id = f"{genome_id}_cds{feature['ordinal']}"
        source_protein_id = feature["source_protein_id"]
        protein_id = (source_protein_id if source_protein_id in unique_protein_ids
                      else f"{transcript_id}.p")
        gene = feature.get("gene_label") or f"cds{feature['ordinal']}"
        # Two CDSs sharing a gene label on opposite strands are two genes, not one.
        strand = "plus" if feature["exons"][0][2] == "+" else "minus"
        rows.append(dict(
            feature,
            genome_id=genome_id,
            transcript_id=transcript_id,
            protein_id=safe_token(protein_id),
            gene_id=f"{genome_id}_{safe_token(gene)}_{strand}",
            gene_name=safe_token(gene),
        ))
    return tuple(rows)


def gtf_transcript_lines(transcript):
    """GTF rows placing one coding transcript on its viral genome."""
    return transcript_lines(
        transcript["genome_id"], transcript["exons"],
        gene_id=transcript["gene_id"],
        gene_name=transcript["gene_name"],
        transcript_id=transcript["transcript_id"],
        protein_id=transcript["protein_id"],
        translation_table=transcript["translation_table"],
        coding_offset=transcript["codon_start"] - 1,
        has_stop_codon=transcript["complete"],
        extra=(("annotation_source", transcript["annotation_source"]),))


def gtf_gene_line(genome_id, gene_id, transcripts):
    """One gene row spanning every transcript that shares a source gene label."""
    return gene_line(genome_id, gene_id, transcripts[0]["gene_name"],
                     [exon for transcript in transcripts for exon in transcript["exons"]])
