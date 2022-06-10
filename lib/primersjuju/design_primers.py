"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuError
from .primer3_interface import Primer3Result, primer3_design
from .target_transcripts import TargetTranscripts, TargetTranscript, ExonFeature
from .transcript_features import transcript_range_to_features

@dataclass
class PrimerDesign:
    """information collected on one primer pair"""
    ppair_id: str
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
    return transcript_range_to_features(target_transcript.features, trans_coords)

def _validate_primer_features(features_5p, features_3p):
    for feature_5p in features_5p:
        for feature_3p in features_3p:
            if feature_5p.trans.overlaps(feature_3p.trans):
                raise PrimersJuJuError(f"primer3 pairs overlap in transcript space {features_5p} and {features_3p}")
            if feature_5p.genome.overlaps(feature_3p.genome):
                raise PrimersJuJuError(f"primer3 pairs overlap in genome space {features_5p} and {features_3p}")


def _build_primer_design(target_transcript, target_id, result_num, primer3_result):
    ppair_id = "{}+pp{}".format(target_id, result_num)
    features_5p = _get_exon_features(target_transcript, primer3_result.PRIMER_LEFT)
    features_3p = _get_exon_features(target_transcript, primer3_result.PRIMER_RIGHT)
    _validate_primer_features(features_5p, features_3p)
    return PrimerDesign(ppair_id, primer3_result, features_5p, features_3p)

def _build_primer_designs(target_transcripts, target_transcript, primer3_results):
    return PrimerDesigns(target_transcripts.target_id,
                         target_transcripts, target_transcript,
                         [_build_primer_design(target_transcript, target_transcripts.target_id, i + 1, r)
                          for i, r in enumerate(primer3_results)])

def design_primers(genome_data, target_transcripts, *, uniqueness_query=None, dump_fh=None):
    """design transcripts """
    target_transcript = target_transcripts.transcripts[0]
    if dump_fh is not None:
        print(64*"=", file=dump_fh)
        target_transcripts.dump(dump_fh)
    primer3_results = primer3_design(target_transcript, dump_fh=dump_fh)

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

def _get_extra_cols(primer_design):
    """get extra BED columns from primer design"""
    extra_cols = []
    for col_name in _primer_bed_columns:
        col = getattr(primer_design.primer3_result, col_name)
        if isinstance(col, float):
            col = "{:.2f}".format(col)
        extra_cols.append(col)
    return extra_cols

def _genomic_order_regions(primer_design):
    "get pair regions genomic order as need for bed"
    if primer_design.features_5p[-1].genome.end < primer_design.features_3p[0].genome.start:
        return primer_design.features_5p, primer_design.features_3p
    else:
        return primer_design.features_3p, primer_design.features_5p

def _build_bed(primer_design, color):
    features_5p, features_3p = _genomic_order_regions(primer_design)
    first = features_5p[0].genome
    last = features_3p[-1].genome

    bed = Bed(first.name, first.start, last.end, primer_design.ppair_id,
              strand=first.strand, thickStart=first.start, thickEnd=last.end,
              itemRgb=color.toRgb8Str(), numStdCols=12, blocks=[],
              extraCols=_get_extra_cols(primer_design))
    for feat in features_5p + features_3p:
        if isinstance(feat, ExonFeature):
            bed.addBlock(feat.genome.start, feat.genome.end)
    return bed

def build_primer_beds(primer_designs, color):
    return sorted([_build_bed(pd, color) for pd in primer_designs.designs],
                  key=Bed.genome_sort_key)
