"""Explicit NCBI genetic codes for validating source ORFs against references.

Tables follow NCBI gc.prt version 4.6:
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
Only codes with native binding fixtures are supported; other IDs fail closed.
A code is never inferred from species or a contig name.
Not a patient variant annotation engine.
"""

from dataclasses import dataclass, field
from itertools import product
from types import MappingProxyType
from typing import Mapping

_TRIPLETS = tuple("".join(c) for c in product("TCAG", repeat=3))


@dataclass(frozen=True)
class GeneticCode:
    """An NCBI translation table and its alternative initiation codons.

    An initiation codon reads as methionine only in the first position;
    elsewhere it keeps its elongation meaning (CTG is leucine in table 1).
    """

    table_id: int
    name: str
    start_codons: frozenset
    codons: Mapping = field(compare=False, repr=False)

    def translate(self, coding, *, start_codon=None):
        """Translate complete codons, including any stop, from a supported start.

        ``start_codon`` declares a source-supported start outside this code's
        initiation codons; like them, it reads as methionine.
        """
        if start_codon is not None:
            if coding[:3] != start_codon:
                raise ValueError(f"ORF declared start codon {start_codon} does not match {coding[:3]}")
        elif coding[:3] not in self.start_codons:
            raise ValueError(
                "Unsupported ORF initiation codon %s for NCBI table %d; "
                "an explicit initiation model is required" % (coding[:3], self.table_id)
            )
        residues = "".join(self.codons[coding[i : i + 3]] for i in range(0, len(coding) - 2, 3))
        return "M" + residues[1:]


def _code(table_id, name, amino_acids, starts):
    return GeneticCode(
        table_id, name,
        frozenset(c for c, s in zip(_TRIPLETS, starts) if s == "M"),
        MappingProxyType(dict(zip(_TRIPLETS, amino_acids))),
    )


GENETIC_CODES = MappingProxyType({code.table_id: code for code in (
    _code(1, "Standard",
          "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "---M------**--*----M---------------M----------------------------"),
    _code(2, "Vertebrate Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
          "----------**--------------------MMMM----------**---M------------"),
    _code(3, "Yeast Mitochondrial",
          "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "----------**----------------------MM---------------M------------"),
)})

# Standard-code codon mapping retained for existing callers.
CODONS = dict(GENETIC_CODES[1].codons)


def genetic_code(table_id):
    """Return a supported NCBI table; unsupported IDs never fall back to table 1."""
    try:
        return GENETIC_CODES[table_id]
    except (KeyError, TypeError):
        raise ValueError(
            "Unsupported ORF translation table %s; no implicit standard-code fallback" % (table_id,)
        ) from None


def read_orf(sequence, start, *, code, start_codon=None, recoded=None, circular=False, initiation=True):
    """Read codons from ``start`` to the first stop; return ``(coding, protein, stopped)``.

    The first codon must be one of the code's initiation codons, or exactly
    ``start_codon`` when the source declares another one; either way it reads
    as methionine. ``recoded`` maps 0-based codon indexes to residues for stop
    codons the source says are read through, such as selenocysteine ``U``.
    A circular sequence is read across its junction at most once around; a
    linear one that ends before a stop gives ``stopped=False``.
    """
    sequence = str(sequence).upper()
    length = len(sequence)
    if not 0 <= start < length:
        raise ValueError(f"Start {start} is outside the RNA sequence")
    recoded = dict(recoded or {})
    first = (sequence[start:] + sequence[:start])[:3] if circular else sequence[start : start + 3]
    if not initiation:
        if start_codon is not None:
            raise ValueError("A fragment without a known start cannot declare a start codon")
    elif start_codon is not None:
        if first != start_codon:
            raise ValueError(f"Expected start codon {start_codon} at {start}, found {first}")
    elif first not in code.start_codons:
        raise ValueError(
            "Unsupported ORF initiation codon %s for NCBI table %d; declare start_codon "
            "for a source-supported non-standard start" % (first, code.table_id))
    residues, coding, offset = [], [], start
    while True:
        if circular:
            if offset - start + 3 > length:
                raise ValueError("No stop codon within one traversal of the circular RNA")
            codon = "".join(sequence[(offset + i) % length] for i in range(3))
        else:
            codon = sequence[offset : offset + 3]
            if len(codon) < 3:
                break
        if any(b not in "ACGT" for b in codon):
            raise ValueError(f"Ambiguous base in codon at {offset % length}")
        index = len(residues)
        residue = "M" if index == 0 and initiation else code.codons[codon]
        coding.append(codon)
        offset += 3
        if residue == "*":
            if index in recoded:
                residues.append(recoded.pop(index))
                continue
            break
        if index in recoded:
            raise ValueError(f"Codon {index} ({codon}) is not a stop codon; it cannot be recoded")
        residues.append(residue)
    if recoded:
        raise ValueError(f"Recoded codons {sorted(recoded)} lie beyond the ORF's stop")
    stopped = bool(coding) and code.codons[coding[-1]] == "*" and len(coding) - 1 == len(residues)
    return "".join(coding), "".join(residues), stopped


def reverse_complement_dna(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate_reference(sequence):
    """Match the original catalogue builder's standard-code initiation rules."""
    residues = []
    for offset in range(0, len(sequence) - 2, 3):
        codon = sequence[offset : offset + 3]
        aa = "M" if offset == 0 and codon in ("ATG", "CTG", "TTG") else CODONS[codon]
        if aa == "*":
            return "".join(residues), True
        residues.append(aa)
    return "".join(residues), False
