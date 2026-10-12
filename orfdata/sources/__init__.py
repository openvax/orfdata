"""Source-specific importers; source measurements retain their native semantics."""


def importer(name):
    """Resolve a packaged importer, never execute an import path from metadata."""
    if name == "pfizer-protein-v1":
        from .pfizer import import_protein

        return import_protein
    raise ValueError(f"No packaged importer for {name!r}")
