"""Gene annotations that PyEnsembl's ``EnsemblRelease`` doesn't cover, pinned so orfdata can install them.

One example is Ensembl's rapid-release annotation of T2T-CHM13v2.0. As with
chain files, :func:`install_genome` (or ``orfdata genomes install``) downloads
and checks the pinned files once; :func:`load_genome` then opens the PyEnsembl
genome offline, with its species set so ``Proteome.build`` and
``extra_orfs`` work on it as on any Ensembl release.
"""

import json
from dataclasses import dataclass, field
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
    assets: Mapping = field(hash=False)
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


def install_genome(name, *, force=False, source_paths=None, cache_root=None, download_options=None,
                   genome_fasta_path_or_url=None):
    """Download a pinned annotation, check each file's size and SHA-256, and index it for PyEnsembl.

    ``source_paths`` maps ``"gtf"``, ``"transcripts"`` or ``"proteins"`` to
    local copies to install instead of downloading, checked the same way.
    ``genome_fasta_path_or_url`` optionally installs chromosome DNA through
    PyEnsembl. This caller-selected source is separate from the checksum-pinned
    annotation assets. Omitting it preserves the annotation-only installation.
    Returns the opened genome.
    """
    genome = _pinned(name)
    paths = None
    if source_paths:
        unknown = set(source_paths) - set(_ROLES)
        if unknown:
            raise ValueError(f"Unknown genome files {sorted(unknown)}; use {', '.join(_ROLES)}")
        paths = {_files(genome)[role]: str(path) for role, path in source_paths.items()}
    from .datasets import _resumable

    _registry(genome, cache_root).download(genome.name, "1", force=force, source_paths=paths,
                                         **_resumable(download_options))
    opened = _with_genome_fasta(_open(genome, cache_root, rebuild=force), genome_fasta_path_or_url)
    if genome_fasta_path_or_url is not None:
        opened.download_genome_fasta(overwrite=force)
        opened.index_genome_fasta()
    return opened


def _install_command(name, cache_root, force=False):
    root = "" if cache_root is None else f" --cache-root {cache_root}"
    return f"'orfdata genomes{root} install {name}{' --force' if force else ''}'"


def load_genome(name, *, cache_root=None, genome_fasta_path_or_url=None):
    """Open an installed annotation offline, optionally with native chromosome DNA.

    Pass the same ``genome_fasta_path_or_url`` used for installation. Loading
    never downloads or opens DNA; a missing installation fails on DNA access
    with PyEnsembl's explicit installation instructions.
    """
    genome = _pinned(name)
    status = _registry(genome, cache_root).inspect(genome.name, "1", verify_files=False).status
    root = "" if cache_root is None else f", cache_root={str(cache_root)!r}"
    if status == "missing":
        raise ValueError(f"The {name} annotation is not installed; install it once with install_genome({name!r}"
                         f"{root}) or {_install_command(name, cache_root)}")
    if status != "available":
        raise ValueError(f"The installed {name} annotation is {status}; install it again with "
                         f"install_genome({name!r}, force=True{root}) or {_install_command(name, cache_root, True)}")
    return _with_genome_fasta(_open(genome, cache_root), genome_fasta_path_or_url)


def _with_genome_fasta(genome, source):
    if source is None:
        return genome
    return type(genome)(**dict(genome.to_dict(), genome_fasta_path_or_url=source))


def _open(genome, cache_root, *, rebuild=False):
    import pickle
    import sqlite3
    from pathlib import Path

    from ._pinned_genome import PinnedEnsemblGenome
    from .datasets import _digest, get_data_dir

    registry = _registry(genome, cache_root)
    # Absolute paths, so the genome works from any directory and when unpickled elsewhere.
    try:
        paths = {role: str(Path(registry.local_path(genome.name, "1", asset=f)).absolute())
                 for role, f in _files(genome).items()}
    except ValueError as error:
        raise ValueError(f"The installed {genome.name} annotation is damaged ({error}); install it again with "
                         f"{_install_command(genome.name, cache_root, True)}") from error
    # The index belongs to these exact files, so a new pin never reuses an old one.
    index = (get_data_dir(cache_root) / "genomes" / genome.name / "pyensembl" / _digest(dict(genome.assets))).absolute()
    opened = PinnedEnsemblGenome(genome.reference_name, genome.annotation_name, genome.annotation_version,
                                 gtf_path_or_url=paths["gtf"], transcript_fasta_paths_or_urls=[paths["transcripts"]],
                                 protein_fasta_paths_or_urls=[paths["proteins"]], cache_directory_path=str(index),
                                 species=genome.species, assembly_accession=genome.assembly_accession,
                                 pinned_name=genome.name)
    if rebuild:  # A forced install rebuilds the index from scratch, even a damaged one.
        import shutil

        shutil.rmtree(index, ignore_errors=True)
    try:
        opened.index()  # Local files only; done once.
    except (sqlite3.Error, EOFError, pickle.UnpicklingError) as error:
        raise ValueError(f"The {genome.name} index is damaged ({error or type(error).__name__}); rebuild it with "
                         f"{_install_command(genome.name, cache_root, True)}") from error
    return opened
