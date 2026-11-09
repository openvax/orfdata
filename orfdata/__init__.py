"""Reference ORFs and evidence, independent of disease and patient variants.

Every public name is importable from here, as ``from orfdata import Proteome``
or ``orfdata.Proteome``, but its submodule is imported the first time the name
is used (PEP 562). ``import orfdata`` and the ``orfdata`` command's tab
completion therefore load neither pandas, NumPy nor PyArrow until something
needs them. Submodules are reachable as attributes too, as ``orfdata.liftover``.
"""

from importlib import import_module
from typing import TYPE_CHECKING

from .version import __version__

# ``viral_subsets`` names both a function and its module. Importing the module here, before anything else
# can, lets the function keep the name, as it did when every name was imported eagerly; the module is cheap.
from .viral_subsets import normalize_viral_subsets, viral_subset_abbreviations, viral_subsets

if TYPE_CHECKING:  # What type checkers and IDEs read; at run time, __getattr__ imports these on first use.
    from .catalog import (
        EVIDENCE_METHODS,
        CatalogStatus,
        ORFCatalog,
        ORFEvidence,
        ORFRecord,
        PeptideEvidence,
        UORFCatalog,
        UORFPeptide,
        UORFRecord,
        catalog_status,
        load_catalog,
        load_uorf_catalog,
        merge_catalogs,
    )
    from .catalog_manifest import CatalogManifest
    from .catalogs import CatalogRegistry, build_catalog
    from .data_resources import data_info, data_path, fetch_data
    from .dataset import QuantitativeDataset
    from .datasets import DatasetRegistry, build_dataset
    from .downloads import DataStatus, list_data
    from .errors import NotInstalledError
    from .evidence import MS, EvidenceRule, RiboSeq
    from .evidence_summary import EvidenceSummary
    from .genomes import (
        PinnedGenome,
        contig_aliases,
        install_genome,
        load_genome,
        pinned_genomes,
        reference_names,
    )
    from .gtf import orf_gtf_lines
    from .initiation_context import (
        KOZAK_PREFERENCES,
        InitiationContext,
        annotate_initiation_context,
    )
    from .install import install_species
    from .liftover import (
        ChainFile,
        Liftover,
        PinnedChain,
        chain_routes,
        download_chain,
        download_chains,
        lift_catalog,
        lift_to,
        load_chain,
        pinned_chains,
    )
    from .measurements import (
        Assay,
        DatasetManifest,
        Measurement,
        MeasurementReference,
        Sample,
        SourceAsset,
        TargetMapping,
    )
    from .model import ORFModel, ReferenceORF, TranscriptORF
    from .pairing import (
        AssayPair,
        MatchedMeasurements,
        OrthologyComparison,
        OrthologyMapping,
        TargetRelation,
        compare_orthologs,
    )
    from .proteome import (
        BindingFailure,
        DerivedProtein,
        ExtraORFs,
        Occurrence,
        ProducerInfo,
        ProducerLink,
        ProducerRecord,
        ProteinDerivation,
        ProteinSequence,
        ProteinSet,
        Proteome,
        ReferenceTranscriptome,
        annotated_orfs,
        extra_orfs,
    )
    from .reference import ReferenceContext, ScopedTarget
    from .reference_mapping import MappedMeasurement, ReferenceMapping
    from .rna import IRES, ProteinProduct, RNAFeature, SourceRNA, products_from_uniprot, rfam_ires
    from .search_reference import (
        SearchMembership,
        SearchReference,
        SearchReferenceSpec,
        SearchSequence,
    )
    from .sequence import GeneticCode, genetic_code
    from .viral import (
        ViralProduct,
        ViralReference,
        viral_orfs,
        viral_reference,
        viral_reference_sequence,
        viral_references,
    )
    from .viral_coverage import (
        download_viral_inventory,
        viral_contexts,
        viral_panel_coverage,
    )
    from .viral_panel import (
        ViralReferenceFiles,
        download_viral_references,
        open_viral_references,
        prune_viral_references,
        select_viral_references,
        viral_annotation_survey,
        viral_genome,
        viral_reference_context,
        viral_reference_summary,
    )
    from .viral_registry import ViralCollection, viral_collection, viral_collections

# The submodule that defines each public name not imported above, as {name: submodule}.
EXPORTED_FROM = {name: module for module, names in {
    "catalog": (
        "EVIDENCE_METHODS", "CatalogStatus", "ORFCatalog", "ORFEvidence", "ORFRecord", "PeptideEvidence",
        "UORFCatalog", "UORFPeptide", "UORFRecord", "catalog_status", "load_catalog", "load_uorf_catalog",
        "merge_catalogs",
    ),
    "catalog_manifest": ("CatalogManifest",),
    "catalogs": ("CatalogRegistry", "build_catalog"),
    "data_resources": ("data_info", "data_path", "fetch_data"),
    "dataset": ("QuantitativeDataset",),
    "datasets": ("DatasetRegistry", "build_dataset"),
    "downloads": ("DataStatus", "list_data"),
    "errors": ("NotInstalledError",),
    "evidence": ("MS", "EvidenceRule", "RiboSeq"),
    "evidence_summary": ("EvidenceSummary",),
    "genomes": ("PinnedGenome", "install_genome", "load_genome", "pinned_genomes", "reference_names",
                "contig_aliases"),
    "gtf": ("orf_gtf_lines",),
    "initiation_context": ("KOZAK_PREFERENCES", "InitiationContext", "annotate_initiation_context"),
    "install": ("install_species",),
    "liftover": (
        "ChainFile", "Liftover", "PinnedChain", "chain_routes", "download_chain", "download_chains",
        "lift_catalog", "lift_to", "load_chain", "pinned_chains",
    ),
    "measurements": (
        "Assay", "DatasetManifest", "Measurement", "MeasurementReference", "Sample", "SourceAsset",
        "TargetMapping",
    ),
    "model": ("ORFModel", "ReferenceORF", "TranscriptORF"),
    "pairing": (
        "AssayPair", "MatchedMeasurements", "OrthologyComparison", "OrthologyMapping", "TargetRelation",
        "compare_orthologs",
    ),
    "proteome": (
        "BindingFailure", "DerivedProtein", "ExtraORFs", "Occurrence", "ProducerInfo", "ProducerLink",
        "ProducerRecord", "ProteinDerivation", "ProteinSequence", "ProteinSet", "Proteome",
        "ReferenceTranscriptome", "annotated_orfs", "extra_orfs",
    ),
    "reference": ("ReferenceContext", "ScopedTarget"),
    "reference_mapping": ("MappedMeasurement", "ReferenceMapping"),
    "rna": ("IRES", "ProteinProduct", "RNAFeature", "SourceRNA", "products_from_uniprot", "rfam_ires"),
    "search_reference": ("SearchMembership", "SearchReference", "SearchReferenceSpec", "SearchSequence"),
    "sequence": ("GeneticCode", "genetic_code"),
    "viral": (
        "ViralProduct", "ViralReference", "viral_orfs", "viral_reference", "viral_reference_sequence",
        "viral_references",
    ),
    "viral_coverage": ("download_viral_inventory", "viral_contexts", "viral_panel_coverage"),
    "viral_panel": (
        "ViralReferenceFiles", "download_viral_references", "open_viral_references", "prune_viral_references",
        "select_viral_references", "viral_annotation_survey", "viral_genome", "viral_reference_context",
        "viral_reference_summary",
    ),
    "viral_registry": ("ViralCollection", "viral_collection", "viral_collections"),
}.items() for name in names}


def __getattr__(name):
    """Import a public name's submodule, or a submodule itself, the first time it is used (PEP 562).

    Python calls this only for names not yet in the package's namespace; the
    name found is stored there, so later uses don't come back here.

    Parameters
    ----------
    name : str
        A name in ``__all__``, or the name of a submodule such as
        ``"liftover"``.

    Returns
    -------
    object
        The public object, or the submodule.

    Raises
    ------
    AttributeError
        If ``name`` is neither a public name nor a submodule.
    """
    module = EXPORTED_FROM.get(name)
    if module is not None:
        value = getattr(import_module(f".{module}", __name__), name)
        globals()[name] = value
        return value
    if not name.startswith("__"):
        try:
            return import_module(f".{name}", __name__)  # The import system sets it as an attribute too.
        except ModuleNotFoundError as error:
            if error.name != f"{__name__}.{name}":  # The submodule exists but needs something missing.
                raise
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """List the package's namespace with every public name, imported yet or not.

    Returns
    -------
    list of str
        Sorted names: those in ``__all__`` and those already defined.
    """
    return sorted(set(globals()) | set(__all__))


__all__ = [
    "ViralProduct",
    "ViralReference",
    "viral_orfs",
    "viral_reference",
    "viral_reference_sequence",
    "viral_references",
    "viral_reference_context",
    "download_viral_inventory",
    "viral_contexts",
    "viral_panel_coverage",
    "viral_subsets",
    "normalize_viral_subsets",
    "viral_subset_abbreviations",
    "ViralReferenceFiles",
    "download_viral_references",
    "open_viral_references",
    "prune_viral_references",
    "select_viral_references",
    "viral_genome",
    "viral_annotation_survey",
    "viral_reference_summary",
    "__version__",
    "Proteome",
    "ProteinSet",
    "ProteinSequence",
    "ReferenceTranscriptome",
    "Occurrence",
    "extra_orfs",
    "annotated_orfs",
    "orf_gtf_lines",
    "ViralCollection",
    "viral_collection",
    "viral_collections",
    "ExtraORFs",
    "BindingFailure",
    "DerivedProtein",
    "ProteinDerivation",
    "ProducerInfo",
    "ProducerLink",
    "ProducerRecord",
    "GeneticCode",
    "genetic_code",
    "SourceRNA",
    "RNAFeature",
    "IRES",
    "ProteinProduct",
    "products_from_uniprot",
    "rfam_ires",
    "EVIDENCE_METHODS",
    "RiboSeq",
    "MS",
    "EvidenceRule",
    "ORFCatalog",
    "ORFRecord",
    "ORFEvidence",
    "CatalogManifest",
    "CatalogRegistry",
    "build_catalog",
    "EvidenceSummary",
    "PeptideEvidence",
    "ORFModel",
    "TranscriptORF",
    "ReferenceORF",
    "ReferenceContext",
    "ScopedTarget",
    "QuantitativeDataset",
    "DatasetRegistry",
    "build_dataset",
    "DataStatus",
    "list_data",
    "data_info",
    "data_path",
    "fetch_data",
    "NotInstalledError",
    "DatasetManifest",
    "SourceAsset",
    "Sample",
    "Assay",
    "Measurement",
    "MeasurementReference",
    "TargetMapping",
    "AssayPair",
    "TargetRelation",
    "MatchedMeasurements",
    "OrthologyMapping",
    "OrthologyComparison",
    "compare_orthologs",
    "ReferenceMapping",
    "MappedMeasurement",
    "SearchReferenceSpec",
    "SearchReference",
    "SearchSequence",
    "SearchMembership",
    "load_catalog",
    "catalog_status",
    "CatalogStatus",
    "merge_catalogs",
    "install_species",
    "PinnedGenome",
    "pinned_genomes",
    "install_genome",
    "load_genome",
    "reference_names",
    "contig_aliases",
    "ChainFile",
    "Liftover",
    "lift_catalog",
    "PinnedChain",
    "pinned_chains",
    "download_chain",
    "download_chains",
    "chain_routes",
    "lift_to",
    "load_chain",
    "UORFCatalog",
    "UORFRecord",
    "UORFPeptide",
    "load_uorf_catalog",
    "KOZAK_PREFERENCES",
    "InitiationContext",
    "annotate_initiation_context",
]
