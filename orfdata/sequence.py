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

    def translate(self, coding):
        """Translate complete codons, including any stop, from a supported start."""
        if coding[:3] not in self.start_codons:
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
