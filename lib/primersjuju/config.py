# Class for configuration objects in the python code based configuration file
# your config.py should import this and construct a instance of
# PrimersJuJuConfig in a variable named 'config'
#
from primersjuju.uniqueness_query import IsPcrServerSpec

class AssemblyConfig:
    "Configuration for a particular assembly"

    def __init__(self, asm_name: str,
                 genome_ispcr_spec: IsPcrServerSpec,
                 transcriptome_ispcr_spec: IsPcrServerSpec):
        self.asm_name = asm_name
        self.genome_ispcr_spec = genome_ispcr_spec
        self.transcriptome_ispcr_spec = transcriptome_ispcr_spec

class PrimersJuJuConfig:
    """Main configuration object. An instance of this object must be create
    and stored in a variable Config"""
    def __init__(self):
        self.assemblies = {}

    def add_assembly(self, asm_config: AssemblyConfig):
        self.assemblies[asm_config.asm_name] = asm_config
