"""Source-specific importers; source measurements retain their native semantics."""


def importer(name):
    """Resolve a packaged importer, never execute an import path from metadata."""
    if name == "pfizer-protein-v1":
        from .pfizer import import_protein

        return import_protein
    if name == "pfizer-canine-evidence-v1":
        from .pfizer_native import import_native

        return import_native
    if name == "copy-yeast-v1":
        from .copy import import_copy

        return import_copy
    if name == "pfizer-rna-v1":
        from .pfizer_rna import import_rna

        return import_rna
    if name == "wang-human-tissue-v1":
        from .wang import import_tissue

        return import_tissue
    if name == "ingolia-yeast-ribo-v1":
        from .ingolia import import_ribo

        return import_ribo
    if name in ("kaabinejadian-dla-v1", "kaabinejadian-dla-v2"):
        from .kaabinejadian import import_ligands

        return import_ligands
    if name == "peptideatlas-build-v1":
        from .peptideatlas import import_build

        return import_build
    raise ValueError(f"No packaged importer for {name!r}")


def search_references(manifest):
    """Return source-declared native search sets; calibration FASTAs are excluded."""
    if manifest.processing.get("importer") == "pfizer-protein-v1":
        from .pfizer import search_references as pfizer_references

        return pfizer_references(manifest)
    if manifest.processing.get("importer") == "pfizer-canine-evidence-v1":
        from .pfizer_native import search_references as native_references

        return native_references(manifest)
    if manifest.processing.get("importer") == "wang-human-tissue-v1":
        from .wang import search_references as wang_references

        return wang_references(manifest)
    return ()
