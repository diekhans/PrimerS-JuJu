"""
Analysis of target transcript.  Includes validation, and trimming of primer
regions to match exons.
"""

class TranscriptRegions:
    def __init__(self, track_name, trans, region_5p_features, region_3p_features):
        self.track_name = track_name
        self.trans = trans
        self.region_5p_features = region_5p_features
        self.region_3p_features = region_3p_features


def _build_transcript_regions(genome_data, primer_target_spec, trans_spec):
    pass  # get_transcript_crange_features();
