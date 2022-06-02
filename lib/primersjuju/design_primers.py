"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuDataError
from .primer3_interface import Primer3Result, primer3_design
from primersjuju.target_transcripts import TargetTranscripts, TargetTranscript, ExonFeature, transcript_range_to_genome

@dataclass
class PrimerDesign:
    """information collected on one primer pair"""
    design_id: str
    primer3_result: Primer3Result
    features_5p: Sequence[ExonFeature]
    features_3p: Sequence[ExonFeature]

@dataclass
class PrimerDesigns:
    """information collected on all primers for a target"""
    target_id: str
    target_transcripts: TargetTranscripts
    target_transcript: TargetTranscript
    designs: Sequence[PrimerDesign]


def _get_exon_features(target_transcript, primer3_coords):
    some_feat = target_transcript.features_5p.region
    start = primer3_coords[0] - 1  # one-based
    end = start + primer3_coords[1]  # has length
    trans_coords = Coords(some_feat.trans.name, start, end,
                          strand=some_feat.trans.strand, size=some_feat.trans.size)
    return transcript_range_to_genome(target_transcript, trans_coords)


def _build_primer_design(target_transcript, design_id, primer3_result):
    return PrimerDesign(design_id, primer3_result,
                        _get_exon_features(target_transcript, primer3_result.PRIMER_LEFT),
                        _get_exon_features(target_transcript, primer3_result.PRIMER_RIGHT))

def _build_primer_designs(target_transcripts, target_transcript, primer3_results):
    return PrimerDesigns(target_transcripts.target_id,
                         target_transcripts, target_transcript,
                         [_build_primer_design(target_transcript, "pp{}".format(i + 1), r)
                          for i, r in enumerate(primer3_results)])

def design_primers(genome_data, target_transcripts, uniqueness_query=None):
    """design transcripts """
    target_transcript = target_transcripts.transcripts[0]
    primer3_results = primer3_design(target_transcript)

    return _build_primer_designs(target_transcripts, target_transcript,
                                 primer3_results)

_primer_bed_columns = (
    'PRIMER_PAIR_PRODUCT_SIZE',    # 1133
    'PRIMER_PAIR_COMPL_ANY_TH',    # 12.498152433966595
    'PRIMER_PAIR_COMPL_END_TH',    # 5.033689592343535
    'PRIMER_PAIR_PENALTY',         # 0.42111006752037383
    'PRIMER_LEFT',                 # (27, 20)
    'PRIMER_LEFT_SEQUENCE',        # 'GGTTCTTCTGCGCTACTGCT'
    'PRIMER_LEFT_END_STABILITY',   # 4.24
    'PRIMER_LEFT_GC_PERCENT',      # 55.0
    'PRIMER_LEFT_HAIRPIN_TH',      # 36.97289132676252
    'PRIMER_LEFT_PENALTY',         # 0.3904716385882807
    'PRIMER_LEFT_SELF_ANY_TH',     # 9.732052967213463
    'PRIMER_LEFT_SELF_END_TH',     # 0.0
    'PRIMER_LEFT_TM',              # 60.39047163858828
    'PRIMER_RIGHT',                # (1159, 20)
    'PRIMER_RIGHT_SEQUENCE',       # 'CAAAAACCCACGCAGACAGG'
    'PRIMER_RIGHT_END_STABILITY',  # 4.0
    'PRIMER_RIGHT_GC_PERCENT',     # 55.0
    'PRIMER_RIGHT_HAIRPIN_TH',     # 0.0
    'PRIMER_RIGHT_PENALTY',        # 0.03063842893209312
    'PRIMER_RIGHT_SELF_ANY_TH',    # 0.0
    'PRIMER_RIGHT_SELF_END_TH',    # 0.0
    'PRIMER_RIGHT_TM',             # 59.96936157106791
)

def _make_bed(primer_design, color):
    "{:.2f}"

def write_primers_bed(primer_designs, color, bed_fh):
    pass
