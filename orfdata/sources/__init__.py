"""Source-specific importers; source measurements retain their native semantics."""


def importer(name):
    """Resolve a packaged importer, never execute an import path from metadata."""
    if name == "pfizer-protein-v1":
        from .pfizer import import_protein

        return import_protein
    if name == "copy-yeast-v1":
        from .copy import import_copy

        return import_copy
    if name == "pfizer-rna-v1":
        from .pfizer_rna import import_rna

        return import_rna
    raise ValueError(f"No packaged importer for {name!r}")


def search_references(manifest):
    """Return source-declared native search sets; calibration FASTAs are excluded."""
    if manifest.processing.get("importer") == "pfizer-protein-v1":
        from .pfizer import search_references as pfizer_references

        return pfizer_references(manifest)
    return ()
