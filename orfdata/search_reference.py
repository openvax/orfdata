"""Immutable native search sequences, separate from abundance and existence evidence."""

import hashlib
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Optional

from .measurements import SourceAsset, _identifier, _Record
from .reference import ReferenceContext, ScopedTarget


@dataclass(frozen=True)
class SearchReferenceSpec(_Record):
    """A source adapter's explicit declaration of a native protein search set."""

    dataset_id: str
    dataset_version: str
    reference: ReferenceContext
    source: SourceAsset
    namespace: str
    format: str
    expected_sequences: Optional[int] = None

    def __post_init__(self):
        for name in ("dataset_id", "dataset_version", "namespace"):
            _identifier(getattr(self, name), name)
        if not isinstance(self.reference, ReferenceContext) or not isinstance(self.source, SourceAsset):
            raise TypeError("Search references require a context and pinned SourceAsset")
        if self.reference.sequence_sha256 != self.source.sha256:
            raise ValueError("Native search context must identify the exact source FASTA bytes")
        if self.format != "uniprot-fasta":
            raise ValueError("Unsupported declared search FASTA format")
        if self.expected_sequences is not None and (type(self.expected_sequences) is not int or self.expected_sequences < 1):
            raise ValueError("Expected sequence count must be positive or unknown")

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        values["reference"] = ReferenceContext(**values["reference"])
        values["source"] = SourceAsset.from_dict(values["source"])
        return cls(**values)


@dataclass(frozen=True)
class SearchSequence(_Record):
    """One original entry; source annotation does not establish protein detection."""

    target: ScopedTarget
    header: str
    fasta_id: str
    sequence: str
    sequence_sha256: str
    source_taxon_id: int
    gene_label: Optional[str]
    database: str
    source_line: int
    decoy: Optional[bool] = None
    control: Optional[bool] = None

    def __post_init__(self):
        if not isinstance(self.target, ScopedTarget) or self.target.kind != "protein":
            raise TypeError("Search sequence requires a native protein target")
        for name in ("header", "fasta_id", "sequence"):
            _identifier(getattr(self, name), name)
        if not re.fullmatch("[A-Z]+", self.sequence):
            raise ValueError("Invalid search sequence characters")
        if hashlib.sha256(self.sequence.encode("ascii")).hexdigest() != self.sequence_sha256:
            raise ValueError("Search sequence hash disagrees with its residues")
        if self.source_taxon_id != self.target.reference.taxon_id:
            raise ValueError("Native search sequence taxon disagrees with its reference")
        if self.database not in ("sp", "tr") or type(self.source_line) is not int or self.source_line < 1:
            raise ValueError("Invalid search sequence source metadata")
        if self.gene_label is not None:
            _identifier(self.gene_label, "gene_label")
        if any(value is not None and type(value) is not bool for value in (self.decoy, self.control)):
            raise TypeError("Decoy/control roles must be explicit booleans or unknown")


@dataclass(frozen=True)
class SearchMembership(_Record):
    """A scoped search-space result, never an abundance or detectability claim."""

    target: ScopedTarget
    status: str
    reason: str
    sequence: Optional[SearchSequence] = None

    def __post_init__(self):
        if not isinstance(self.target, ScopedTarget):
            raise TypeError("Search membership requires a scoped target")
        _identifier(self.reason, "membership reason")
        if self.status not in ("present", "not_in_search", "unknown"):
            raise ValueError("Unknown search membership status")
        if self.status == "present":
            if not isinstance(self.sequence, SearchSequence) or self.sequence.target.key != self.target.key:
                raise ValueError("Present membership requires its exact native sequence")
        elif self.sequence is not None:
            raise ValueError("Unknown/absent membership cannot assign a native sequence")

    @classmethod
    def from_dict(cls, data):
        values = dict(data)
        target = dict(values["target"])
        target["reference"] = ReferenceContext(**target["reference"])
        values["target"] = ScopedTarget(**target)
        if values.get("sequence") is not None:
            values["sequence"] = SearchSequence.from_dict(values["sequence"])
        return cls(**values)


@dataclass(frozen=True)
class SearchReference:
    """An opened, verified source FASTA with immutable entries and exact lookups."""

    path: Path
    spec: SearchReferenceSpec
    _entries: Mapping = field(init=False, repr=False, compare=False)

    def __init__(self, path, spec):
        if not isinstance(spec, SearchReferenceSpec):
            raise TypeError("SearchReference requires a declared SearchReferenceSpec")
        object.__setattr__(self, "path", Path(path).resolve())
        object.__setattr__(self, "spec", spec)
        digest = hashlib.sha256()
        size = 0
        with self.path.open("rb") as handle:
            for block in iter(lambda: handle.read(1048576), b""):
                size += len(block)
                digest.update(block)
        if size != spec.source.size or digest.hexdigest() != spec.source.sha256:
            raise ValueError("Search FASTA size/SHA-256 disagrees with pinned source")
        object.__setattr__(self, "_entries", MappingProxyType(self._read()))

    @property
    def reference(self):
        return self.spec.reference

    @property
    def source(self):
        return self.spec.source

    def __len__(self):
        return len(self._entries)

    def entries(self):
        """Iterate every original sequence, including unreported representatives."""
        return iter(self._entries.values())

    def get(self, identifier):
        """Look up an exact accession within this source's declared namespace."""
        _identifier(identifier, "search sequence identifier")
        return self._entries.get(identifier)

    def membership(self, target):
        """Check a native protein ID; unresolved groups and other namespaces stay unknown."""
        if not isinstance(target, ScopedTarget):
            raise TypeError("Search membership requires an explicitly scoped target")
        if target.reference.key != self.reference.key:
            raise ValueError("Search membership target belongs to a different reference")
        if target.kind != "protein":
            return SearchMembership(target, "unknown", "Protein/group/ORF attribution requires an explicit native-protein mapping")
        if target.namespace != self.spec.namespace:
            return SearchMembership(target, "unknown", "Identifier namespace requires an explicit mapping")
        sequence = self.get(target.identifier)
        if sequence is None:
            return SearchMembership(target, "not_in_search", "Exact accession absent from this pinned FASTA")
        return SearchMembership(target, "present", "Exact accession present in this pinned FASTA; detection is not established", sequence)

    def _read(self):
        entries, header, pieces, source_line = {}, None, [], None

        def save():
            if header is None:
                return
            fasta_id = header.split()[0]
            fields = fasta_id.split("|")
            taxon = re.search(r"(?:^| )OX=(\d+)(?: |$)", header)
            if len(fields) != 3 or fields[0] not in ("sp", "tr") or not all(fields[1:]) or taxon is None:
                raise ValueError("Invalid declared UniProt search FASTA header")
            if int(taxon.group(1)) != self.reference.taxon_id:
                raise ValueError("Search entry taxon differs from the declared native reference")
            identifier = fields[1]
            sequence = "".join(pieces)
            if identifier in entries or not sequence or not re.fullmatch("[A-Z]+", sequence):
                raise ValueError("Duplicate/empty/invalid native search sequence")
            gene = re.search(r"(?:^| )GN=(\S+)", header)
            entries[identifier] = SearchSequence(
                ScopedTarget(self.reference, "protein", identifier, self.spec.namespace), header, fasta_id,
                sequence, hashlib.sha256(sequence.encode("ascii")).hexdigest(), int(taxon.group(1)),
                gene.group(1) if gene else None, fields[0], source_line,
            )

        with self.path.open(encoding="utf-8") as handle:
            for line, text in enumerate(handle, 1):
                text = text.rstrip("\r\n")
                if not text.strip():
                    continue
                if text.startswith(">"):
                    save()
                    header, pieces, source_line = text[1:], [], line
                    if not header.strip():
                        raise ValueError("Empty search FASTA header")
                elif header is None:
                    raise ValueError("Search FASTA sequence precedes its header")
                else:
                    pieces.append(text.strip())
        save()
        if not entries or (self.spec.expected_sequences is not None and len(entries) != self.spec.expected_sequences):
            raise ValueError("Search FASTA sequence count disagrees with declared coverage")
        return entries


def reference_registry(spec, cache_root=None):
    """Use a separate immutable DataCache store without altering abundance caches."""
    from .datasets import _datacache, _digest, get_data_dir

    asset = spec.source
    definitions = {spec.dataset_id: dict(default_version=spec.dataset_version, versions={
        spec.dataset_version: {asset.asset_id: dict(url=asset.url, sha256=asset.sha256,
                                                   size=asset.size, decompress=False)},
    })}
    destination = (get_data_dir(cache_root) / "references" / spec.dataset_id /
                   spec.dataset_version / _digest(spec.to_dict()))
    return _datacache().VersionedDatasetRegistry(definitions, store_path=lambda name, version: destination)
