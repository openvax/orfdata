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

from .cache import content_digest, get_data_dir, pinned_registry
from .errors import NotInstalledError

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
    #: Other names for the assembly's sequences, e.g. ``hs1``: each gives the
    #: accession of the file's mitochondrion and, under ``without``, the
    #: contigs its files lack or carry as a different molecule.
    names: Mapping = field(default_factory=dict, hash=False)
    #: Names of other versions of the assembly, each with why it isn't this one.
    other_versions: Mapping = field(default_factory=dict, hash=False)
    #: The assembly's sequences: the annotation's name for each, its UCSC,
    #: GenBank and RefSeq names, and its length.
    contigs: tuple = field(default=(), hash=False)
    #: Where :attr:`contigs` came from, and its SHA-256.
    contigs_source: Mapping = field(default_factory=dict, hash=False)
    #: Mitochondria some of the assembly's files carry instead of its own, by
    #: accession, each with its annotation packaged in ORFdata. The genome
    #: includes their genes, on a contig named by the accession so that they
    #: can't be mistaken for the assembly's own mitochondrion.
    mitochondria: Mapping = field(default_factory=dict, hash=False)


def _pins():
    return json.loads(files("orfdata").joinpath("data/genomes.json").read_text())


def pinned_genomes(*, species=None):
    """The pinned annotations, without downloading any."""
    from .model import _species

    found = tuple(PinnedGenome(name, **dict(entry, contigs=tuple(entry.get("contigs", ()))))
                  for name, entry in sorted(_pins().items()))
    return found if species is None else tuple(g for g in found if g.species == _species(species))


def reference_names():
    """Every name the pinned genomes answer to, and which genome each one opens.

    A name is matched exactly, ignoring case: finding one inside a file path
    is the caller's job. A pinned genome answers to its own name, its
    ``reference_name``, its assembly accession and the names in its ``names``
    (for T2T-CHM13v2.0: ``hs1``, ``GCF_009914755.1``, the T2T Consortium's
    ``chm13v2.0*`` analysis sets and T2T-CHM13v1.1). Names of other versions
    of an assembly, such as ``T2T-CHM13v1.0``, are listed with no genome.

    Returns
    -------
    pandas.DataFrame
        One row per name, with columns:

        - ``name``: the name, as its source spells it.
        - ``genome``: the pinned genome it opens, e.g. ``"t2t-chm13"``, or
          None for another version of the assembly.
        - ``mitochondrion``: the accession of the mitochondrion in files with
          this name, or None when they have none. For T2T-CHM13v2.0 it is
          CHM13's own (``CP068254.1``), except in the ``_rCRS`` analysis set,
          whose chrM is the rCRS (``NC_012920.1``). Both are 16,569 bp, but
          CP068254.1 begins at what the rCRS numbers 577.
        - ``note``: what the name is or, with no genome, why it isn't one.
    """
    import pandas as pd

    rows = []
    for genome in pinned_genomes():
        mitochondrion = next((c["genbank"] for c in genome.contigs if c["name"] == "MT"), None)
        own = {genome.name: "orfdata's name", genome.reference_name: "The annotation's reference_name"}
        if genome.assembly_accession:
            own[genome.assembly_accession] = "The assembly's accession"
        for name, note in own.items():
            if name not in genome.names:
                rows.append(dict(name=name, genome=genome.name, mitochondrion=mitochondrion, note=note))
        rows += [dict(name=name, genome=genome.name, mitochondrion=entry["mitochondrion"], note=entry["note"])
                 for name, entry in genome.names.items()]
        rows += [dict(name=name, genome=None, mitochondrion=None, note=reason)
                 for name, reason in genome.other_versions.items()]
    return pd.DataFrame(rows, columns=["name", "genome", "mitochondrion", "note"], dtype=object)


def _find(name):
    """``(genome, its names entry, None)`` for a name a pinned genome answers to, ``(None, None, why not)``
    for another version of an assembly, and None for a name orfdata doesn't know."""
    if not isinstance(name, str):
        raise TypeError(f"A genome name is a string, not {type(name).__name__}")
    key = name.strip().lower()
    for genome in pinned_genomes():
        for other, reason in genome.other_versions.items():
            if other.lower() == key:
                return None, None, reason
        for alias, entry in genome.names.items():
            if alias.lower() == key:
                return genome, entry, None
        if key in {genome.name.lower(), genome.reference_name.lower(), (genome.assembly_accession or "").lower()}:
            return genome, {}, None
    return None


def _lookup(name):
    """The pinned genome ``name`` opens, and its ``names`` entry if it has one."""
    found = _find(name)
    if found is None:
        pinned = ", ".join(g.name for g in pinned_genomes()) or "none"
        raise ValueError(f"No pinned genome {name!r}; pinned: {pinned}. 'orfdata genomes names' lists every "
                         "name each one answers to")
    genome, entry, reason = found
    if genome is None:
        raise ValueError(f"{name} has no pinned genome. {reason}")
    return genome, entry


def _pinned(name):
    return _lookup(name)[0]


def contig_aliases(name):
    """The sequences in files named ``name``, under each of their names, with their lengths.

    Parameters
    ----------
    name : str
        Any name :func:`reference_names` gives a genome, e.g. ``"t2t-chm13"``,
        ``"hs1"`` or ``"chm13v2.0_maskedY_rCRS"``.

    Returns
    -------
    pandas.DataFrame
        One row per sequence, with columns ``name`` (the annotation's name,
        e.g. ``"1"``), ``ucsc`` (``"chr1"``), ``genbank`` (``"CP068277.2"``),
        ``refseq`` (``"NC_060925.1"``, or None) and ``length``. Sequences the
        name's files lack, or carry as a different molecule, are left out: no
        ``Y`` for ``chm13v2.0_noY``, and no ``MT`` for the ``_rCRS`` analysis
        set, whose ``chrM`` is the rCRS rather than CHM13's mitochondrion, or
        for ``GCF_009914755.1``, which has none. A mitochondrion in
        :attr:`PinnedGenome.mitochondria` is listed instead, under the
        annotation's name for it: the ``_rCRS`` set's ``chrM`` is
        ``NC_012920.1``. ``frame.attrs`` holds the genome's name and the
        source of the table.

    Raises
    ------
    ValueError
        For a name of another version of the assembly, such as
        ``T2T-CHM13v1.0``, with why its coordinates differ.
    """
    import pandas as pd

    genome, entry = _lookup(name)
    without = set(entry.get("without", ()))
    rows = [c for c in genome.contigs if c["name"] not in without]
    if entry.get("mitochondrion") in genome.mitochondria:
        rows.append(genome.mitochondria[entry["mitochondrion"]]["contig"])
    frame = pd.DataFrame(rows, columns=["name", "ucsc", "genbank", "refseq", "length"], dtype=object)
    frame = frame.astype({"length": "int64"})
    frame.attrs.update(genome=genome.name, source=dict(genome.contigs_source))
    return frame


def _files(genome):
    """Each role's file name, as published: PyEnsembl reads the format from it."""
    return {role: genome.assets[role]["url"].rsplit("/", 1)[-1] for role in _ROLES}


def genome_registry(genome, cache_root=None):
    """Return the DataCache registry that installs a pinned annotation's files.

    The release is named ``genome.name``, version ``"1"``, and its store is
    ``genomes/<name>/<digest of the file pins>`` in the cache, so a new pin
    never reuses old files. Nothing is downloaded until the registry is asked to.

    Parameters
    ----------
    genome : PinnedGenome
        The annotation, as :func:`pinned_genomes` lists it.
    cache_root : str or os.PathLike, optional
        A cache root instead of ``ORFDATA_CACHE``, ``OPENVAX_DATA_CACHE`` or the
        default.

    Returns
    -------
    datacache.VersionedDatasetRegistry
        Holds the GTF, transcript FASTA and protein FASTA under their published
        file names.
    """
    assets = {_files(genome)[role]: dict(genome.assets[role]) for role in _ROLES}
    return pinned_registry(genome.name, "1", assets,
                           get_data_dir(cache_root) / "genomes" / genome.name / content_digest(assets))


def install_genome(name, *, force=False, source_paths=None, cache_root=None, download_options=None,
                   genome_fasta_path_or_url=None):
    """Download a pinned annotation, check each file's size and SHA-256, and index it for PyEnsembl.

    ``name`` is any name :func:`reference_names` gives a genome, such as
    ``"t2t-chm13"`` or ``"hs1"``. ``source_paths`` maps ``"gtf"``,
    ``"transcripts"`` or ``"proteins"`` to local copies to install instead of
    downloading, checked the same way.
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
    genome_registry(genome, cache_root).download(genome.name, "1", force=force, source_paths=paths,
                                                 **(download_options or {}))
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

    ``name`` is any name :func:`reference_names` gives a genome: ``"hs1"`` and
    ``"chm13v2.0_maskedY_rCRS"`` both open ``"t2t-chm13"``. Pass the same
    ``genome_fasta_path_or_url`` used for installation. Loading never
    downloads or opens DNA; a missing installation fails on DNA access with
    PyEnsembl's explicit installation instructions.

    Raises
    ------
    NotInstalledError
        When the annotation isn't installed, naming the command that installs it.
    ValueError
        For an unknown name, a name of another version of the assembly, or a
        damaged installation.
    """
    genome = _pinned(name)
    pinned = genome.name
    status = genome_registry(genome, cache_root).inspect(pinned, "1", verify_files=False).status
    root = "" if cache_root is None else f", cache_root={str(cache_root)!r}"
    named = "" if name == pinned else f" (for {name})"
    if status == "missing":
        raise NotInstalledError(f"The {pinned} annotation{named} is not installed; install it once with "
                                f"install_genome({pinned!r}{root}) or {_install_command(pinned, cache_root)}")
    if status != "available":
        raise ValueError(f"The installed {pinned} annotation is {status}; install it again with "
                         f"install_genome({pinned!r}, force=True{root}) or "
                         f"{_install_command(pinned, cache_root, True)}")
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

    registry = genome_registry(genome, cache_root)
    # Absolute paths, so the genome works from any directory and when unpickled elsewhere.
    try:
        paths = {role: str(Path(registry.local_path(genome.name, "1", asset=f)).absolute())
                 for role, f in _files(genome).items()}
    except ValueError as error:
        raise ValueError(f"The installed {genome.name} annotation is damaged ({error}); install it again with "
                         f"{_install_command(genome.name, cache_root, True)}") from error
    # The index belongs to these exact files, so a new pin never reuses an old one.
    pinned = dict(genome.assets)
    if genome.mitochondria:
        pinned["mitochondria"] = {accession: dict(m["assets"]) for accession, m in genome.mitochondria.items()}
    index = (get_data_dir(cache_root) / "genomes" / genome.name / "pyensembl" / content_digest(pinned)).absolute()
    if rebuild:  # A forced install rebuilds the index from scratch, even a damaged one.
        import shutil

        shutil.rmtree(index, ignore_errors=True)
    gtf, transcripts, proteins = paths["gtf"], [paths["transcripts"]], [paths["proteins"]]
    metadata = None
    if genome.mitochondria:
        gtf, extra = _with_mitochondria(genome, paths["gtf"], index)
        transcripts += extra["transcripts"]
        proteins += extra["proteins"]
        # PyEnsembl recognizes MT and chrM by name; an accession needs saying.
        metadata = dict(contig_compartments={m["contig"]["name"]: "mitochondrial"
                                             for m in genome.mitochondria.values()})
    opened = PinnedEnsemblGenome(genome.reference_name, genome.annotation_name, genome.annotation_version,
                                 gtf_path_or_url=gtf, transcript_fasta_paths_or_urls=transcripts,
                                 protein_fasta_paths_or_urls=proteins, cache_directory_path=str(index),
                                 species=genome.species, assembly_accession=genome.assembly_accession,
                                 pinned_name=genome.name, reference_metadata=metadata)
    try:
        opened.index()  # Local files only; done once.
    except (sqlite3.Error, EOFError, pickle.UnpicklingError) as error:
        raise ValueError(f"The {genome.name} index is damaged ({error or type(error).__name__}); rebuild it with "
                         f"{_install_command(genome.name, cache_root, True)}") from error
    return opened


def _with_mitochondria(genome, gtf, index):
    """Write the genome's GTF with its packaged mitochondria appended, and their FASTAs, into ``index``.

    The files are copied out of the package so the genome still opens, and
    unpickles, after ORFdata is upgraded. Each is read only when its bytes are
    the pinned ones. The GTF is the Ensembl file's gzip stream followed by each
    mitochondrion's, which gzip readers read as one file.
    """
    import os

    from .reference_table import packaged_bytes

    index.mkdir(parents=True, exist_ok=True)
    extra = {"transcripts": [], "proteins": []}
    layers = []
    for m in genome.mitochondria.values():
        for role, asset in m["assets"].items():
            raw = packaged_bytes(asset["resource"], asset["sha256"])
            if role == "gtf":
                layers.append(raw)
                continue
            path = index / asset["resource"].rsplit("/", 1)[-1]
            if not path.exists() or path.stat().st_size != len(raw):
                path.write_bytes(raw)
            extra[role].append(str(path))
    name = gtf.rsplit("/", 1)[-1].removesuffix(".gtf.gz")
    merged = index / f"{name}+{'+'.join(genome.mitochondria)}.gtf.gz"
    size = os.path.getsize(gtf) + sum(len(raw) for raw in layers)
    if not merged.exists() or merged.stat().st_size != size:
        partial = merged.with_name(merged.name + ".partial")
        with open(gtf, "rb") as source, open(partial, "wb") as target:
            while chunk := source.read(1 << 20):
                target.write(chunk)
            for raw in layers:
                target.write(raw)
        os.replace(partial, merged)
    return str(merged), extra
