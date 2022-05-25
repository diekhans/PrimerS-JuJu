"""
Primer selection for a target.
"""
from dataclasses import dataclass
from pycbio.sys.objDict import ObjDict
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

def primer3_query(primer_target):
    seq_args = ObjDict()
    seq_args.SEQUENCE_ID = primer_target.target_id

    #    SEQUENCE_TEMPLATE (nucleotide sequence; default empty)

# The sequence from which to choose primers. The sequence must be presented 5'
# -> 3' (i.e, in the normal way). In general, the bases may be upper or lower
# case, but lower case letters are treated specially if
# PRIMER_LOWERCASE_MASKING is set.  The entire sequence MUST be all on a
# single line. (In other words, the sequence cannot span several lines.)
# SEQUENCE_INCLUDED_REGION (interval list; default empty)
