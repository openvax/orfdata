"""Backwards-compatible Genome subclass for legacy serialized pinned annotations."""

from pyensembl import Genome


class PinnedEnsemblGenome(Genome):
    """Legacy subclass kept so previously serialized pinned genomes still load.

    ``pyensembl.Genome`` (>=2.35.0) natively accepts and serializes ``species``
    and ``assembly_accession``; :func:`orfdata.load_genome` returns plain
    ``pyensembl.Genome`` instances.
    """

    def __init__(self, *args, pinned_name=None, **kwargs):
        self.pinned_name = pinned_name
        super().__init__(*args, **kwargs)

    def to_dict(self):
        return dict(super().to_dict(), pinned_name=self.pinned_name)
