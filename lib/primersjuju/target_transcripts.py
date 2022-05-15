"""
Target transcripts analysis.  Includes validation, and trimming of primer
regions to match exons.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuDataError
from pycbio.hgdata.bed import Bed


class Feature(Coords):
    "annotation feature"

    def intersect(self, other):
        if isinstance(self, ExonRegion):
            return ExonRegion(super().intersect(other))
        elif isinstance(self, IntronRegion):
            return ExonRegion(super().intersect(other))
        else:
            raise ValueError(f"bad object type: {type(other)}")

class ExonRegion(Feature):
    "exon in a model"
    pass

class IntronRegion(Feature):
    "intron in a model"
    pass

@dataclass
class TargetTranscriptFeatures:
    """Coords and features of a transcript for 5' or 3' region.  """
    region: Coords
    features: Sequence[Feature]

@dataclass
class TargetTranscript:
    "target transcripts for RT-PCR with features within primer regions"
    track_name: str
    trans: Bed
    features_5p: TargetTranscriptFeatures
    features_3p: TargetTranscriptFeatures

    def __str__(self):
        return f"({self.track_name}, {self.trans.name})"

def _get_features_bounds(features):
    f0 = features[0]
    return Coords(f0.name, f0.start, features[-1].end, f0.strand, f0.size)

def _block_features(trans, region, csize, prev_blk, blk, features):
    "get intron and exon feature intersection with region"

    def _mk_feature(feat_cls, start, end):
        "create Feature for intersecting the region"
        return feat_cls(trans.chrom,
                        max(start, region.start), min(end, region.end),
                        strand=trans.strand, size=csize)

    if prev_blk is not None:
        if (prev_blk.end < region.end) and (blk.start > region.start):
            features.append(_mk_feature(IntronRegion, prev_blk.end, blk.start))
    if (blk.start < region.end) and (blk.end > region.start):
        features.append(_mk_feature(ExonRegion, blk.start, blk.end))

def get_transcript_region_features(genome_data, trans, region):
    """Given a chromosome region in a transcript, generate of a list feature
    coords in that region.
    """
    csize = genome_data.get_chrom_size(trans.chrom)

    features = []
    prev_blk = None
    for blk in trans.blocks:
        _block_features(trans, region, csize, prev_blk, blk, features)
        prev_blk = blk
    return features

def _get_regions_strand(trans, region_5p, region_3p):
    # swap regions if needed so they match strand
    pos_order = region_5p.end < region_3p.start
    if ((trans.strand == '+' and not pos_order) or
        (trans.strand == '-' and pos_order)):
        region_5p, region_3p = region_3p, region_5p
    return region_5p, region_3p

def _build_transcript_features(genome_data, target_transcript, region):
    features = get_transcript_region_features(genome_data, target_transcript, region)
    exon_cnt = filter(lambda f: isinstance(f, ExonRegion), features)
    if exon_cnt == 0:
        raise PrimersJuJuDataError(f"transcript {target_transcript} has no exons in region {region}")
    return TargetTranscriptFeatures(_get_features_bounds(features), features)

def _build_target_transcript(genome_data, primer_target_spec, trans_spec):
    "build transcript with initial regions trimmed to exons"
    trans = genome_data.get_track(trans_spec.track_name).read_by_name(trans_spec.trans_id)
    region_5p, region_3p = _get_regions_strand(trans, primer_target_spec.region_5p,
                                               primer_target_spec.region_3p)
    return TargetTranscript(primer_target_spec.track_name, trans,
                            _build_transcript_features(genome_data, trans, region_5p),
                            _build_transcript_features(genome_data, trans, region_3p))

def _target_transcripts_build(genome_data, primer_target_spec):
    target_transcripts = []
    for trans_track in primer_target_spec.tracks.values():
        for trans_spec in primer_target_spec.tracks[trans_track]:
            target_transcripts.append(_build_target_transcript(genome_data, primer_target_spec, trans_spec))
    return target_transcripts

def _validate_strand(target_transcripts):
    # must do before other checks, as primer region swap could confuse error messages
    ttrans0 = target_transcripts[0]
    for ttrans in target_transcripts[1:]:
        if ttrans.region_5p.strand != ttrans0.trans.strand:
            raise PrimersJuJuDataError(f"transcript on different strands: {ttrans} vs {ttrans0}")

def _find_transcripts_common_region(target_transcripts, feats_func):
    region = feats_func(target_transcripts[0])
    for ttrans in target_transcripts[1:]:
        region = region.intersect(feats_func(ttrans).region)
    return region

def _adjust_transcript_region_features(target_transcript, common_region, feats_func):
    orig_features = feats_func(target_transcript).features
    adj_features = []
    for ofeat in orig_features:
        if ofeat.overlaps(common_region):
            adj_features.append(ofeat.intersect(common_region))

    # must still end in exons
    if not (isinstance(adj_features[0], ExonRegion) and
            isinstance(adj_features[-1], ExonRegion)):
        raise PrimersJuJuDataError(f"transcript {target_transcript} does not end in exons once trimmed to the common region between all transcripts in {common_region}")

    # this updates transcript features
    feats_func(target_transcript,
               TargetTranscriptFeatures(common_region, adj_features))

def _adjust_transcripts_region_features(target_transcripts, feats_func):
    common_region = _find_transcripts_common_region(target_transcripts, feats_func)
    for target_transcript in target_transcripts:
        _adjust_transcript_region_features(target_transcript, common_region, feats_func)
    return common_region

def _get_introns(target_transcript_features):
    return [feat for feat in target_transcript_features.features
            if isinstance(feat, IntronRegion)]

def _validate_intron_count(target_transcripts, region, transcripts_introns, feats_func):
    intron_count = None
    for trans, introns in target_transcripts, transcripts_introns:
        if len(introns) > 1:
            raise PrimersJuJuDataError(f"target transcript may span no more than two exons, {trans} spans ")
        if intron_count is None:
            intron_count = len(introns)
        if len(introns) != intron_count:
            raise PrimersJuJuDataError(f"target transcript may span no more than two exons, {trans} spans ")
    return intron_count

def _validate_intron_bounds(target_transcripts, region, transcripts_introns, feats_func):
    trans0 = target_transcripts[0]
    introns0 = transcripts_introns[0][0]
    for trans, introns in target_transcripts[1:], transcripts_introns[1:]:
        if introns[0] != introns0:
            raise PrimersJuJuDataError(f"intron-spanning transcripts must all have the same introns {trans} has {introns[0]} which {trans0} has {introns0}")

def _validate_intron_spanning(target_transcripts, region, feats_func):
    # all transcripts must span the same single intron or all must
    # be contained in exon

    transcripts_introns = [_get_introns(feats_func(tt))
                           for tt in target_transcripts]
    intron_count = _validate_intron_count(target_transcripts, region, transcripts_introns, feats_func)
    if intron_count > 0:
        _validate_intron_bounds(target_transcripts, region, transcripts_introns, feats_func)

def _adjust_transcripts_features(target_transcripts):
    """Adjust region to be the same for all transcript, which maybe have been
    adjusted to exon bounds.  Validate that they are still ending in exons"""

    # these get or set region
    def features_5p_access(t, new_value=None):
        if new_value is not None:
            t.features_5p = new_value
        return t.features_5p

    def features_3p_access(t, new_value=None):
        if new_value is not None:
            t.features_3p = new_value
        return t.features_3p

    region_5p = _adjust_transcripts_region_features(target_transcripts, features_5p_access)
    region_3p = _adjust_transcripts_region_features(target_transcripts, features_3p_access)

    _validate_intron_spanning(target_transcripts, region_5p, features_5p_access)
    _validate_intron_spanning(target_transcripts, region_3p, features_3p_access)


@dataclass
class TargetTranscripts:
    """Target tracksuit for primer regions with, with features in regions for
    each transcript.  Used when validating consistency of transcripts.  Target
    transcript for primer experiment.
    """
    target_id: str
    region_5p: Coords
    region_3p: Coords
    transcripts: [TargetTranscript] = None

def _do_target_transcripts_build(genome_data, primer_target_spec):
    target_transcripts = _target_transcripts_build(genome_data, primer_target_spec)

    _validate_strand(target_transcripts)
    _adjust_transcripts_features(target_transcripts)

    return TargetTranscripts(primer_target_spec.target_id, target_transcripts)

def target_transcripts_build(genome_data, primer_target_spec):
    """build TargetTranscripts object to a give primer and validity and consistency of
    the transcripts.
    """
    try:
        return _do_target_transcripts_build(genome_data, primer_target_spec)
    except PrimersJuJuDataError as ex:
        raise PrimersJuJuDataError(f"target {primer_target_spec.target_id} failed") from ex
