"""A PyEnsembl genome class for pinned annotations; imported only where PyEnsembl is installed."""

from pyensembl import Genome
from pyensembl.species import find_species_by_name


class PinnedEnsemblGenome(Genome):
    """A PyEnsembl genome that keeps its species and assembly accession when pickled or copied."""

    def __init__(self, *args, species=None, assembly_accession=None, pinned_name=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.species = find_species_by_name(species) if isinstance(species, str) else species
        self.assembly_accession = assembly_accession
        self.pinned_name = pinned_name

    def to_dict(self):
        return dict(super().to_dict(), species=self.species.latin_name, assembly_accession=self.assembly_accession,
                    pinned_name=self.pinned_name)
