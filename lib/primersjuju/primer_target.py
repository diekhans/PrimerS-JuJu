"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuDataError
from .target_transcripts import TargetTranscripts

@dataclass
class PrimerTarget:
    """one primer-pair target"""
    target_transcripts: TargetTranscripts
    region_seq_5p: str
    region_seq_3p: str

    @property
    def region_5p(self) -> Coords:
        return self.region_5p

    @property
    def region_3p(self) -> Coords:
        return self.region_3p

def primer_target_build(genome_data, target_transcripts):
    region_seq_5p = genome_data.get_genome_seq(target_transcripts.region_5p)
    region_seq_3p = genome_data.get_genome_seq(target_transcripts.region_3p)

    return PrimerTarget(target_transcripts, region_seq_5p, region_seq_3p)
