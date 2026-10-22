"""Discover ORFdata resources and their local installation state."""

from dataclasses import asdict, dataclass
from typing import Optional

from .reference import canonical_species, species_taxon


@dataclass(frozen=True)
class DownloadStatus:
    """Identity, availability and acquisition route for one ORFdata resource."""

    kind: str
    name: str
    version: Optional[str]
    title: str
    species: tuple[str, ...]
    taxon_ids: tuple[int, ...]
    contexts: tuple[str, ...]
    status: str
    path: Optional[str]
    origin: Optional[str]
    acquire: Optional[str]
    command: Optional[str]
    restricted: bool = False

    @property
    def installed(self):
        """Whether the resource is usable from local files."""
        return self.status == "available"

    def to_dict(self):
        """Return a JSON-ready representation, as other OpenVax status objects do."""
        return dict(asdict(self), installed=self.installed)


# The earlier development name remains importable while DownloadStatus is the
# canonical public name.
Downloadable = DownloadStatus


def _viral(cache_root, verify_files):
    from .viral_coverage import PANEL
    from .viral_panel import (
        BUNDLE_VERSION,
        VIRAL_COLLECTION,
        VIRAL_CONTEXTS,
        VIRAL_TARGET_SPECIES,
        VIRAL_TARGET_TAXON_ID,
        _bundle_path,
        _existing_bundle,
        viral_reference_bundle,
    )

    path = _bundle_path(cache_root, VIRAL_TARGET_SPECIES, VIRAL_COLLECTION)
    if not path.is_dir():
        status = "missing"
    else:
        try:
            if verify_files:
                viral_reference_bundle(
                    species=VIRAL_TARGET_SPECIES, collection=VIRAL_COLLECTION,
                    cache_root=cache_root)
            else:
                _existing_bundle(
                    path, species=VIRAL_TARGET_SPECIES, collection=VIRAL_COLLECTION)
        except FileNotFoundError:
            status = "missing"
        except OSError:
            status = "inaccessible"
        except ValueError:
            status = "corrupt"
        else:
            status = "available"
    return (DownloadStatus(
        kind="viral",
        name=VIRAL_COLLECTION,
        version=BUNDLE_VERSION,
        title=f"GDC {PANEL} viral genomes and paired coding/protein sequences",
        species=(VIRAL_TARGET_SPECIES,),
        taxon_ids=(VIRAL_TARGET_TAXON_ID,),
        contexts=VIRAL_CONTEXTS,
        status=status,
        path=str(path),
        origin="download" if status == "available" else None,
        acquire="download",
        command="orfdata download viral --species human",
    ),)


def _catalogs(cache_root, verify_files):
    from .catalog import catalog_status
    from .catalogs import CatalogRegistry

    registry = CatalogRegistry(cache_root=cache_root)
    rows = []
    for item in catalog_status(registry=registry):
        inspection = None
        if item.acquire is None:
            status = "available"
        else:
            inspection = registry.inspect(item.name, item.version, verify_files=verify_files)
            status = inspection.status
        command = None
        if item.acquire is not None:
            command = f"orfdata download catalog {item.name} --version {item.version}"
        rows.append(DownloadStatus(
            kind="catalog",
            name=item.name,
            version=item.version,
            title=item.title,
            species=(item.species,),
            taxon_ids=(species_taxon(item.species),),
            contexts=(),
            status=status,
            path=str(inspection.path) if inspection is not None else None,
            origin=item.installed if status == "available" else None,
            acquire=item.acquire,
            command=command,
        ))
    return tuple(rows)


def _datasets(cache_root, verify_files):
    from .datasets import DatasetRegistry

    registry = DatasetRegistry(cache_root=cache_root)
    rows = []
    for manifest in registry.list():
        name, version = manifest.dataset_id, manifest.version
        published = registry.published(name, version)
        inspection = registry.inspect(name, version, verify_files=verify_files)
        species = tuple(sorted({canonical_species(reference.species)
                                for reference in manifest.references}))
        taxa = tuple(sorted({reference.taxon_id for reference in manifest.references}))
        acquire = "download" if published else "build"
        rows.append(DownloadStatus(
            kind="dataset",
            name=name,
            version=version,
            title=manifest.title,
            species=species,
            taxon_ids=taxa,
            contexts=(),
            status=inspection.status,
            path=str(inspection.path),
            origin=("download" if published else "build") if inspection.status == "available" else None,
            acquire=acquire,
            command=f"orfdata download dataset {name} --version {version}",
        ))
    return tuple(rows)


def _references(cache_root, verify_files):
    from .datasets import DatasetRegistry

    registry = DatasetRegistry(cache_root=cache_root)
    rows = []
    for manifest in registry.list():
        for spec in registry.search_references(manifest.dataset_id, manifest.version):
            inspection = registry.inspect_reference(
                manifest.dataset_id, manifest.version, reference=spec.reference.key,
                asset_id=spec.source.asset_id, verify_files=verify_files)
            rows.append(DownloadStatus(
                kind="reference",
                name=f"{manifest.dataset_id}:{spec.source.asset_id}",
                version=manifest.version,
                title=f"{manifest.title}: native {spec.namespace} protein search reference",
                species=(canonical_species(spec.reference.species),),
                taxon_ids=(spec.reference.taxon_id,),
                contexts=(),
                status=inspection.status,
                path=str(inspection.path),
                origin="download" if inspection.status == "available" else None,
                acquire="download",
                command=(f"orfdata download reference {manifest.dataset_id} --version {manifest.version} "
                         f"--reference {spec.reference.key} --asset {spec.source.asset_id}"),
            ))
    return tuple(rows)


def _genomes(cache_root, verify_files):
    from .genomes import _registry, pinned_genomes

    rows = []
    for genome in pinned_genomes():
        inspection = _registry(genome, cache_root).inspect(
            genome.name, "1", verify_files=verify_files)
        rows.append(DownloadStatus(
            kind="genome",
            name=genome.name,
            version=genome.annotation_version,
            title=genome.title,
            species=(genome.species,),
            taxon_ids=(species_taxon(genome.species),),
            contexts=(),
            status=inspection.status,
            path=str(inspection.path),
            origin="download" if inspection.status == "available" else None,
            acquire="download",
            command=f"orfdata download genome {genome.name}",
        ))
    return tuple(rows)


def _chains(cache_root, verify_files):
    from .liftover import _chain_registry, pinned_chains

    rows = []
    for chain in pinned_chains():
        inspection = _chain_registry(chain, cache_root).inspect(
            chain.name, "1", verify_files=verify_files)
        accept = " --accept-terms" if chain.restricted else ""
        rows.append(DownloadStatus(
            kind="chain",
            name=chain.name,
            version="1",
            title=f"{chain.source_build} to {chain.target_build} chain ({chain.provider})",
            species=(chain.species,),
            taxon_ids=(species_taxon(chain.species),),
            contexts=(),
            status=inspection.status,
            path=str(inspection.path),
            origin="download" if inspection.status == "available" else None,
            acquire="download",
            command=(f"orfdata download chains {chain.source_build} {chain.target_build} "
                     f"--species {chain.species}{accept}"),
            restricted=chain.restricted,
        ))
    return tuple(rows)


def list_downloads(*, kind=None, species=None, contexts=None, context_match="any",
                   installed=None, cache_root=None, verify_files=False):
    """List known data resources and local state without acquiring anything.

    ``kind`` may be ``viral``, ``catalog``, ``dataset``, ``reference``,
    ``genome`` or ``chain``. ``species`` accepts ORFdata's species aliases.
    ``installed`` can be true or false to select local resources or
    unavailable ones. ``verify_files=True`` hashes pinned files and
    validates the viral bundle.
    """
    kinds = ("viral", "catalog", "dataset", "reference", "genome", "chain")
    if kind == "all":
        kind = None
    if kind is not None and kind not in kinds:
        raise ValueError(f"Unknown data kind {kind!r}; choose from {', '.join(kinds)}")
    if installed is not None and type(installed) is not bool:
        raise ValueError("installed must be True, False or None")
    if contexts is None:
        selected_contexts = ()
    else:
        from .viral_coverage import normalize_viral_contexts

        selected_contexts = normalize_viral_contexts(contexts, allow_unknown=True)
    if context_match not in ("any", "all"):
        raise ValueError("context_match must be 'any' or 'all'")
    selected_species = None
    if species is not None:
        if species_taxon(species) is None:
            raise ValueError(f"Unknown species alias {species!r}")
        selected_species = canonical_species(species)

    finders = {
        "viral": _viral,
        "catalog": _catalogs,
        "dataset": _datasets,
        "reference": _references,
        "genome": _genomes,
        "chain": _chains,
    }
    rows = []
    for current in kinds if kind is None else (kind,):
        rows.extend(finders[current](cache_root, verify_files))
    if selected_species is not None:
        rows = [row for row in rows if selected_species in row.species]
    if selected_contexts:
        predicate = any if context_match == "any" else all
        rows = [row for row in rows
                if predicate(context in row.contexts for context in selected_contexts)]
    if installed is not None:
        rows = [row for row in rows if row.installed is installed]
    return tuple(sorted(rows, key=lambda row: (row.kind, row.name, row.version or "")))


def downloadable_data(**kwargs):
    """Compatibility alias for :func:`list_downloads`."""
    return list_downloads(**kwargs)
