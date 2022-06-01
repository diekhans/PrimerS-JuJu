"""
Query for uniqueness in genome and transcriptome
"""


class UniquenessQuery:
    """interfaces to UCSC isPCR to query for uniqueness."""
    def __init__(self, ucscAssembly, genomeHost=None, genomePort=None,
                 transcriptome=None, transcriptomePort=None):
        self.ucscAssembly = ucscAssembly
        self.genomeHost = genomeHost
        self.genomePort = genomePort
        self.transcriptome = transcriptome
        self.transcriptomePort = transcriptomePort
