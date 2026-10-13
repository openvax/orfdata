"""Native Ensembl search entries retain versioned IDs, alternate loci and stops."""

import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from orfdata import (
    ReferenceContext,
    SearchReference,
    SearchReferenceSpec,
    SearchSequence,
    SourceAsset,
)

ROOT = Path(__file__).parent / "data/search_reference"


def spec_for(path, format, count=3):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    source = SourceAsset("search.fasta.gz" if format.endswith("gzip") else "search.fasta",
                         "https://example.invalid/ensembl-fixture", digest, path.stat().st_size,
                         "Ensembl source sequence slice; original provenance in fixture notice")
    reference = ReferenceContext(9606, "human", "GRCh38", "ensembl", "83",
                                 assembly_accession="GCA_000001405.20", sequence_sha256=digest)
    return SearchReferenceSpec("human-search-fixture", "1", reference, source,
                               "ensembl-protein-versioned", format, count)


@pytest.mark.parametrize("compressed", (False, True))
def test_exact_real_ensembl_entries(compressed, tmp_path):
    native = ROOT / "ensembl83-human.fasta.gz"
    expected = json.loads((ROOT / "ensembl83-expected.json").read_text())
    assert hashlib.sha256(native.read_bytes()).hexdigest() == expected["slice_sha256"]
    path = native
    format = "ensembl-pep-fasta-gzip" if compressed else "ensembl-pep-fasta"
    if not compressed:
        path = tmp_path / "native.fasta"
        path.write_bytes(gzip.decompress(native.read_bytes()))
    spec = spec_for(path, format)
    search = SearchReference(path, spec)
    assert SearchReferenceSpec.from_dict(spec.to_dict()) == spec
    assert len(search) == 3
    for row, entry in zip(expected["records"], search.entries()):
        assert entry.header == row["header"] and entry.sequence == row["sequence"]
        assert entry.fasta_id == row["header"].split()[0]
        assert entry.target.identifier == entry.fasta_id
        assert entry.database == "ensembl" and entry.source_taxon_id is None
        assert entry.target.reference.taxon_id == 9606
        assert entry.decoy is None and entry.control is None
        assert SearchSequence.from_dict(entry.to_dict()) == entry
        assert search.membership(entry.target).status == "present"
        assert search.get(entry.fasta_id.split(".")[0]) is None  # No silent version removal.
    assert search.get("ENSP00000474693.2").sequence == "*LWC*L"
    assert "CHR_HSCHR" in search.get("ENSP00000484462.1").header
    protein = search.get("ENSP00000000233.5")
    assert protein.gene_label == "ENSG00000004059.10"
    other = replace(protein.target.reference, taxon_id=10090, species="mouse")
    with pytest.raises(ValueError, match="different reference"):
        search.membership(replace(protein.target, reference=other))
    with pytest.raises(ValueError, match="assembly"):
        SearchReference(path, replace(spec, reference=replace(spec.reference, genome_build="GRCh37")))
    with pytest.raises(ValueError, match="sequence count"):
        SearchReference(path, replace(spec, expected_sequences=4))


@pytest.mark.parametrize("change,error", (("gene", "header"), ("transcript", "header"),
                                          ("sequence", "invalid"), ("duplicate", "Duplicate")))
def test_invalid_ensembl_record_fails(change, error, tmp_path):
    text = gzip.decompress((ROOT / "ensembl83-human.fasta.gz").read_bytes()).decode()
    if change in ("gene", "transcript"):
        text = text.replace(" " + change + ":", " missing:")
    elif change == "sequence":
        lines = text.splitlines()
        lines[1] += "?"
        text = "\n".join(lines) + "\n"
    else:
        text += text
    path = tmp_path / "source.fasta"
    path.write_text(text)
    with pytest.raises(ValueError, match=error):
        SearchReference(path, spec_for(path, "ensembl-pep-fasta"))
