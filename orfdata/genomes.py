"""Gene annotations that PyEnsembl's ``EnsemblRelease`` doesn't cover, pinned so orfdata can install them.

One example is Ensembl's rapid-release annotation of T2T-CHM13v2.0. As with
chain files, :func:`install_genome` (or ``orfdata genomes install``) downloads
and checks the pinned files once; :func:`load_genome` then opens the PyEnsembl
genome offline, with its species set so ``Proteome.build`` and
``extra_orfs`` work on it as on any Ensembl release.
"""

import json
from dataclasses import dataclass
from importlib.resources import files
from typing import Mapping, Optional

_ROLES = ("gtf", "transcripts", "proteins")


@dataclass(frozen=True)
class PinnedGenome:
    """An annotation orfdata can install: its assembly name, annotation and files pinned by URL, size and SHA-256."""

    name: str
    species: str
    reference_name: str
    annotation_name: str
    annotation_version: str
    title: str
    assets: Mapping
    provider: str
    terms: str
    terms_url: str
    assembly_accession: Optional[str] = None


def _pins():
    return json.loads(files("orfdata").joinpath("data/genomes.json").read_text())


def pinned_genomes(*, species=None):
    """The pinned annotations, without downloading any."""
    from .model import _species

    found = tuple(PinnedGenome(name, **entry) for name, entry in sorted(_pins().items()))
    return found if species is None else tuple(g for g in found if g.species == _species(species))


def _pinned(name):
    found = {g.name: g for g in pinned_genomes()}
    if name not in found:
        raise ValueError(f"No pinned genome {name!r}; pinned: {', '.join(found) or 'none'}")
    return found[name]


def _files(genome):
    """Each role's file name, as published: PyEnsembl reads the format from it."""
    return {role: genome.assets[role]["url"].rsplit("/", 1)[-1] for role in _ROLES}


def _registry(genome, cache_root):
    from .datasets import _datacache, _digest, get_data_dir

    assets = {_files(genome)[role]: dict(genome.assets[role]) for role in _ROLES}
    return _datacache().VersionedDatasetRegistry(
        {genome.name: dict(default_version="1", versions={"1": assets})},
        store_path=lambda name, version: get_data_dir(cache_root) / "genomes" / name / _digest(assets))


def install_genome(name, *, force=False, source_paths=None, cache_root=None):
    """Download a pinned annotation, check each file's size and SHA-256, and index it for PyEnsembl.

    ``source_paths`` maps ``"gtf"``, ``"transcripts"`` or ``"proteins"`` to
    local copies to install instead of downloading, checked the same way.
    Returns the opened genome.
    """
    genome = _pinned(name)
    paths = None
    if source_paths:
        unknown = set(source_paths) - set(_ROLES)
        if unknown:
            raise ValueError(f"Unknown genome files {sorted(unknown)}; use {', '.join(_ROLES)}")
        paths = {_files(genome)[role]: str(path) for role, path in source_paths.items()}
    _registry(genome, cache_root).download(genome.name, "1", force=force, source_paths=paths)
    return _open(genome, cache_root, rebuild=force)


def load_genome(name, *, cache_root=None):
    """The installed pinned annotation as a PyEnsembl genome, offline; indexed on first use."""
    genome = _pinned(name)
    status = _registry(genome, cache_root).inspect(genome.name, "1", verify_files=False).status
    if status == "missing":
        raise ValueError(f"The {name} annotation is not installed; install it once with install_genome({name!r}) "
                         f"or 'orfdata genomes install {name}'")
    if status != "available":
        raise ValueError(f"The installed {name} annotation is {status}; install it again with "
                         f"install_genome({name!r}, force=True) or 'orfdata genomes install {name} --force'")
    return _open(genome, cache_root)


def _open(genome, cache_root, *, rebuild=False):
    import sqlite3

    from ._pinned_genome import PinnedEnsemblGenome
    from .datasets import _digest, get_data_dir

    registry = _registry(genome, cache_root)
    try:
        paths = {role: str(registry.local_path(genome.name, "1", asset=f)) for role, f in _files(genome).items()}
    except ValueError as error:
        raise ValueError(f"The installed {genome.name} annotation is damaged ({error}); install it again with "
                         f"'orfdata genomes install {genome.name} --force'") from error
    # The index belongs to these exact files, so a new pin never reuses an old one.
    index = get_data_dir(cache_root) / "genomes" / genome.name / "pyensembl" / _digest(dict(genome.assets))
    opened = PinnedEnsemblGenome(genome.reference_name, genome.annotation_name, genome.annotation_version,
                                 gtf_path_or_url=paths["gtf"], transcript_fasta_paths_or_urls=[paths["transcripts"]],
                                 protein_fasta_paths_or_urls=[paths["proteins"]], cache_directory_path=str(index),
                                 species=genome.species, assembly_accession=genome.assembly_accession)
    if rebuild:  # A forced install rebuilds the index from scratch, even a damaged one.
        import shutil

        shutil.rmtree(index, ignore_errors=True)
    try:
        opened.index()  # Local files only; done once.
    except sqlite3.Error as error:
        raise ValueError(f"The {genome.name} index is damaged ({error}); rebuild it with "
                         f"'orfdata genomes install {genome.name} --force'") from error
    return opened
