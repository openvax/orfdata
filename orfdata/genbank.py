"""Read the GenBank flat files NCBI EFetch returns, with no third-party parser.

ORFdata needs one producer's output and a small part of the format: the
accession version, whether the molecule is circular, the sequence, and the
feature table with its locations and qualifiers. Keeping that here avoids a
dependency on a library whose data model loses what the file says — a
location's order matters for a spliced CDS, and whether it runs past the end of
the record matters for a gene that crosses a circular genome's origin.

Anything outside that subset — a remote accession, a ``123^124`` site or a
``(1.5)`` uncertain position — raises rather than being guessed at.
"""

import re
from dataclasses import dataclass

_FEATURE_KEY = re.compile(r"^ {5}(\S+)\s*(.*)$")
_LOCUS = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+(?:bp|aa)\b", re.IGNORECASE)
_QUALIFIER = re.compile(r'^/([^=]+)(?:=(.*))?$')
_RANGE = re.compile(r"^([<>]?)(\d+)(?:\.\.([<>]?)(\d+))?$")
_SEQUENCE_LINE = re.compile(r"^\s*\d+((?:\s+[A-Za-z]+)+)\s*$")
# Qualifiers whose wrapped lines are sequence, so they close up without a break.
_SEQUENCE_QUALIFIERS = ("translation", "peptide", "anticodon")


@dataclass(frozen=True)
class Range:
    """One span, 1-based and inclusive; ``before``/``after`` are GenBank's ``<`` and ``>``."""

    start: int
    end: int
    before: bool = False
    after: bool = False


@dataclass(frozen=True)
class Complement:
    """The reverse strand of what it contains."""

    location: object


@dataclass(frozen=True)
class Join:
    """Spans read as one sequence, in the order the file gives them."""

    locations: tuple

    def __init__(self, locations):
        object.__setattr__(self, "locations", tuple(locations))


@dataclass(frozen=True)
class Order(Join):
    """Spans the file does not promise are joined; read in the given order."""


@dataclass(frozen=True)
class Feature:
    """One feature-table entry: its key, its location and its qualifiers in order."""

    kind: str
    location: object
    qualifiers: tuple


@dataclass(frozen=True)
class Record:
    """One GenBank record: what ORFdata reads from a flat file."""

    accession_version: str
    sequence: str
    circular: bool
    features: tuple
    name: str = ""
    definition: str = ""


def positions(location):
    """Every coordinate a location names, for checking it against its record."""
    kind = type(location).__name__
    if kind == "Range":
        return (location.start, location.end)
    if kind == "Complement":
        return positions(location.location)
    return tuple(value for part in location.locations for value in positions(part))


def parse_location(text):
    """Parse a GenBank location descriptor into its structure."""
    text = "".join(text.split())
    for name, kind in (("complement", Complement), ("join", Join), ("order", Order)):
        if text.startswith(name + "(") and text.endswith(")"):
            inner = text[len(name) + 1:-1]
            if kind is Complement:
                return Complement(parse_location(inner))
            return kind([parse_location(part) for part in _split_parts(inner)])
    match = _RANGE.fullmatch(text)
    if match is None:
        raise ValueError(f"Unsupported GenBank location: {text!r}")
    opening, start, closing, end = match.groups()
    if end is None:  # A single position, which either marker can qualify.
        return Range(int(start), int(start), before=opening == "<", after=opening == ">")
    if opening == ">" or closing == "<":
        raise ValueError(f"Unsupported GenBank location markers: {text!r}")
    return Range(int(start), int(end), before=opening == "<", after=closing == ">")


def _split_parts(text):
    """Split a comma-separated location list, respecting nested parentheses."""
    parts, depth, current = [], 0, []
    for character in text:
        if character == "," and depth == 0:
            parts.append("".join(current))
            current = []
            continue
        depth += (character == "(") - (character == ")")
        current.append(character)
    parts.append("".join(current))
    if depth or not all(parts):
        raise ValueError(f"Unbalanced GenBank location list: {text!r}")
    return parts


def _is_open(raw, continuations):
    """Whether a quoted qualifier value is still waiting for its closing quote."""
    if raw is None or not raw.startswith('"'):
        return False
    text = "".join([raw, *continuations])
    # Pairs of quotes are literal; the value closes on the first unpaired one.
    return len(text) < 2 or (len(text) - len(text.rstrip('"'))) % 2 == 0


def _qualifier_value(name, raw, continuations):
    """One qualifier's value, joining continuation lines the way its kind needs.

    A wrapped sequence closes up with nothing between the lines; wrapped prose
    keeps the line break, which readers normalize to a space.
    """
    if raw is None:
        return None
    separator = "" if name in _SEQUENCE_QUALIFIERS else "\n"
    if not raw.startswith('"'):
        return separator.join([raw, *continuations]).strip()
    # A quoted value may run over lines and may contain "" for a literal quote.
    text = separator.join([raw, *continuations])
    if not text.endswith('"') or len(text) < 2:
        raise ValueError(f"Unterminated GenBank qualifier value: {raw!r}")
    return text[1:-1].replace('""', '"')


def _features(lines):
    """Every feature in the table, with multi-line locations and qualifiers rejoined."""
    features, key, location, qualifiers, pending = [], None, None, [], None

    def close_qualifier():
        if pending is not None:
            name, raw, continuations = pending
            qualifiers.append((name, _qualifier_value(name, raw, continuations)))

    def close_feature():
        close_qualifier()
        if key is not None:
            features.append(Feature(key, parse_location(location), tuple(qualifiers)))

    for line in lines:
        match = _FEATURE_KEY.match(line)
        if match:  # A new feature starts at column 6.
            close_feature()
            key, location, qualifiers, pending = match.group(1), match.group(2), [], None
            continue
        if key is None:
            continue
        stripped = line.strip()
        # A line inside an unterminated quoted value belongs to it, even when it
        # happens to begin with the character that starts a qualifier.
        if stripped.startswith("/") and not (pending is not None and _is_open(*pending[1:])):
            close_qualifier()
            qualifier = _QUALIFIER.match(stripped)
            if qualifier is None:
                raise ValueError(f"Unreadable GenBank qualifier: {stripped!r}")
            pending = (qualifier.group(1), qualifier.group(2), [])
        elif pending is not None:
            pending[2].append(stripped)
        else:  # A location continued onto the next line.
            location += stripped
    close_feature()
    return tuple(features)


def parse(text):
    """Read one GenBank flat file.

    The record is checked against what it says about itself: its LOCUS line
    declares a length, and a sequence that does not match it is a truncated or
    damaged record rather than a shorter molecule.
    """
    lines = text.splitlines()
    accession_version, circular, name, declared = None, False, "", None
    feature_lines, sequence, definition, section = [], [], [], None
    for line in lines:
        if line.startswith("LOCUS"):
            circular = " circular " in f" {' '.join(line.split()[1:])} "
            locus = _LOCUS.match(line)
            if locus is not None:
                name, declared = locus.group(1), int(locus.group(2))
        elif line.startswith("VERSION"):
            accession_version = line.split()[1]
        elif line.startswith("DEFINITION"):
            section = "definition"
            definition.append(line[len("DEFINITION"):].strip())
            continue
        elif line.startswith("FEATURES"):
            section = "features"
            continue
        elif line.startswith(("ORIGIN", "CONTIG", "BASE COUNT")):
            section = "sequence" if line.startswith("ORIGIN") else None
            continue
        elif line.startswith("//"):
            break
        elif line and not line.startswith(" "):
            section = None
        if section == "features":
            feature_lines.append(line)
        elif section == "definition":
            definition.append(line.strip())
        elif section == "sequence":
            match = _SEQUENCE_LINE.match(line)
            if match is None:
                raise ValueError(f"Unreadable GenBank sequence line: {line!r}")
            sequence.append("".join(match.group(1).split()))
    if accession_version is None:
        raise ValueError("GenBank record has no VERSION line")
    bases = "".join(sequence).upper()
    if bases and declared is not None and len(bases) != declared:
        raise ValueError(
            f"GenBank record {accession_version} declares {declared} bases but carries "
            f"{len(bases)}"
        )
    features = _features(feature_lines)
    if bases:
        for feature in features:
            outside = [value for value in positions(feature.location)
                       if not 1 <= value <= len(bases)]
            if outside:
                raise ValueError(
                    f"GenBank record {accession_version} places its {feature.kind} feature at "
                    f"{outside[0]}, outside its {len(bases)} bases"
                )
    return Record(accession_version=accession_version, sequence=bases, circular=circular,
                  features=features, name=name,
                  # NCBI writes a record's one-line description as a sentence; its
                  # own titles elsewhere drop the full stop, and so does this.
                  definition=" ".join(definition).strip().rstrip("."))
