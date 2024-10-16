# Class for configuration objects in the python code based configuration file
# your config.py should import this and construct a instance of
# PrimersJuJuConfig in a variable named 'config'
#
from primersjuju.uniqueness_query import UniquenessQuery, IsPcrServerSpec
from primersjuju.genome_data import GenomeData

class GenomeConfig:
    "Configuration for a particular assembly"

    def __init__(self,
                 genome_data: GenomeData,
                 genome_ispcr_spec: IsPcrServerSpec,
                 transcriptome_ispcr_spec: IsPcrServerSpec):
        self.genome_data = genome_data
        self.genome_ispcr_spec = genome_ispcr_spec
        self.transcriptome_ispcr_spec = transcriptome_ispcr_spec
        self.__uniqueness_query = None   # lazy

    @property
    def genome_name(self):
        return self.genome_data.genome_name

    @property
    def uniqueness_query(self) -> UniquenessQuery:
        "lazy get/create UniquenessQuery object, or None if not configured"
        if ((self.__uniqueness_query is None) and
            ((self.genome_ispcr_spec is not None) or (self.transcriptome_ispcr_spec is not None))):
            self.__uniqueness_query = UniquenessQuery(self.genome_data,
                                                      self.genome_ispcr_spec, self.transcriptome_ispcr_spec)
        return self.__uniqueness_query

class Primer3Config:
    """Options to pass to primer3, see primer3 manual for a description.
    """
    def __init__(self):
        self.num_5_prime_strong_match = 0
        self.min_5_prime_overlap_of_junction = 8
        self.min_3_prime_overlap_of_junction = 8
        self.max_poly_x = None
        self.misprime_lib = None
        self.mishyb_lib = None

_default_primer3_config = Primer3Config()

class PrimersJuJuConfig:
    """Main configuration object. An instance of this object must be create
    and stored in a variable Config"""
    def __init__(self, *, primer3_config=_default_primer3_config):
        self.genomes = {}
        self.primer3 = primer3_config

        # derived, not set in config file
        self.genome = None  # selected genome, set at load from genomes

    def add_genome(self, genome_config: GenomeConfig):
        self.genomes[genome_config.genome_name] = genome_config
