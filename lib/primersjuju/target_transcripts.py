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
            return ExonRegion(*super().intersect(other))
        elif isinstance(self, IntronRegion):
            return IntronRegion(*super().intersect(other))
        else:
            raise ValueError(f"bad object type: {type(other)}")

class ExonRegion(Feature):
    "exon in a model"
    pass

class IntronRegion(Feature):
    "intron in a model"
    pass

def _features_get_bounds(features):
    f0 = features[0]
    return Coords(f0.name, f0.start, features[-1].end, f0.strand, f0.size)

def _features_to_str(features):
    return '[' + ("\n ".join([repr(f) for f in features])) + ']'

def _features_sort(features):
    "sort into transcription order"
    return sorted(features, key=lambda f: f.start if f.strand == '+' else -f.end)

def _features_count(features, ftype):
    # True == 1, False == 0
    return sum([isinstance(f, ftype) for f in features])

def _primer_region_check_features(desc, track_name, trans_id, features):
    "check that features are sane, desc is used in error messages"
    exon_cnt = _features_count(features, ExonRegion)
    intron_cnt = _features_count(features, IntronRegion)
    if not (((exon_cnt == 1) and (intron_cnt == 0)) or ((exon_cnt == 2) and (intron_cnt == 1))):
        raise PrimersJuJuDataError(f"{desc} for transcript ({track_name}, {trans_id}) must contain either one exon, or two exons and an intron: {_features_to_str(features)}")
    if not (isinstance(features[0], ExonRegion) and
            isinstance(features[-1], ExonRegion)):
        raise PrimersJuJuDataError(f"{desc} transcript ({track_name}, {trans_id}) primer region does not end in exons: {_features_to_str(features)}")


@dataclass
class PrimerRegionFeatures:
    """Coords and features of a transcript for 5' or 3' region.  """
    region: Coords
    features: Sequence[Feature]

@dataclass
class TargetTranscript:
    "a target transcript for RT-PCR with features within primer regions"
    track_name: str
    bed: Bed
    features_5p: PrimerRegionFeatures
    features_3p: PrimerRegionFeatures
    amplicon: str = None  # includes target regions

    @property
    def trans_id(self):
        return self.bed.name

    @property
    def strand(self):
        return self.bed.strand

    @property
    def region_5p(self):
        return self.features_5p.region

    @property
    def region_3p(self):
        return self.features_3p.region

    def __str__(self):
        return f"({self.track_name}, {self.trans_id})"

    @property
    def amplicon_range(self):
        if self.strand == '+':
            return Coords(self.region_5p.name,
                          self.region_5p.start, self.region_3p.end,
                          self.region_5p.strand, self.region_5p.size)
        else:
            return Coords(self.region_5p.name,
                          self.region_3p.start, self.region_5p.end,
                          self.region_5p.strand, self.region_5p.size)


@dataclass
class TargetTranscripts:
    """Target tracksuit for primer regions with, with features in regions for
    each transcript.  Used when validating consistency of transcripts.  Target
    transcript for primer experiment.
    """
    target_id: str
    region_5p: Coords
    region_3p: Coords
    sequence_5p: str
    sequence_3p: str
    transcripts: [TargetTranscript]

    def get_transcript(self, track_name, trans_id):
        for t in self.transcripts:
            if (t.track_name == track_name) and (t.bed.name == trans_id):
                return t
        raise PrimersJuJuDataError(f"({track_name}, {trans_id}) not found in {self.target_id}")

def _block_features(trans_bed, region, csize, prev_blk, blk, features):
    "get intron and exon feature intersection with a genomic range"

    def _mk_feature(feat_cls, start, end):
        "create Feature for intersecting the region"
        return feat_cls(trans_bed.chrom,
                        max(start, region.start), min(end, region.end),
                        strand=trans_bed.strand, size=csize)

    if prev_blk is not None:
        if (prev_blk.end < region.end) and (blk.start > region.start):
            features.append(_mk_feature(IntronRegion, prev_blk.end, blk.start))
    if (blk.start < region.end) and (blk.end > region.start):
        features.append(_mk_feature(ExonRegion, blk.start, blk.end))

def get_transcript_region_features(genome_data, trans_bed, region):
    """Given a chromosome region in a transcript, generate of a list feature
    coords in that region.
    """
    csize = genome_data.get_chrom_size(trans_bed.chrom)

    features = []
    prev_blk = None
    for blk in trans_bed.blocks:
        _block_features(trans_bed, region, csize, prev_blk, blk, features)
        prev_blk = blk
    return features

def _get_regions_strand(trans, region_5p, region_3p):
    # swap regions if needed so they match strand
    pos_order = region_5p.end < region_3p.start
    if ((trans.strand == '+' and not pos_order) or
        (trans.strand == '-' and pos_order)):
        region_5p, region_3p = region_3p, region_5p
    return region_5p, region_3p

def _build_transcript_features(genome_data, track_name, trans_bed, region):
    features = get_transcript_region_features(genome_data, trans_bed, region)
    _primer_region_check_features("initially specified primer region", track_name, trans_bed.name, features)
    return PrimerRegionFeatures(_features_get_bounds(features), features)

def _build_target_transcript(genome_data, primer_target_spec, trans_spec):
    "build transcript with initial regions trimmed to exons"
    trans_bed = genome_data.get_track(trans_spec.trans_track).read_by_name(trans_spec.trans_id)
    region_5p, region_3p = _get_regions_strand(trans_bed, primer_target_spec.region_5p,
                                               primer_target_spec.region_3p)
    return TargetTranscript(trans_spec.trans_track, trans_bed,
                            _build_transcript_features(genome_data, trans_spec.trans_track, trans_bed, region_5p),
                            _build_transcript_features(genome_data, trans_spec.trans_track, trans_bed, region_3p))

def _target_transcripts_build(genome_data, primer_target_spec):
    target_transcripts = []
    for track in primer_target_spec.tracks.values():
        for trans_spec in track.values():
            target_transcripts.append(_build_target_transcript(genome_data, primer_target_spec, trans_spec))
    return target_transcripts

def _validate_strand(target_transcripts):
    # must do before other checks, as primer region swap could confuse error messages
    ttrans0 = target_transcripts[0]
    for ttrans in target_transcripts[1:]:
        if ttrans.region_5p.strand != ttrans0.strand:
            raise PrimersJuJuDataError(f"transcript on different strands: {ttrans} vs {ttrans0}")

def _find_transcripts_common_region(target_transcripts, feats_func):
    region = feats_func(target_transcripts[0]).region
    for ttrans in target_transcripts[1:]:
        region = region.intersect(feats_func(ttrans).region)
    return region

def _adjust_transcript_region_features(target_transcript, common_region, feats_func):
    orig_features = feats_func(target_transcript).features
    adj_features = []
    for ofeat in orig_features:
        if ofeat.overlaps(common_region):
            adj_features.append(ofeat.intersect(common_region))

    _primer_region_check_features("common adjusted primer region", target_transcript.track_name, target_transcript.trans_id, adj_features)

    # this updates transcript features
    feats_func(target_transcript,
               PrimerRegionFeatures(common_region, adj_features))

def _adjust_transcripts_region_features(target_transcripts, feats_func):
    common_region = _find_transcripts_common_region(target_transcripts, feats_func)
    for target_transcript in target_transcripts:
        _adjust_transcript_region_features(target_transcript, common_region, feats_func)
    return common_region

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
    return region_5p, region_3p

def _get_exon_intersect_seq(genome_data, bed, blk, crange):
    coords = Coords(bed.chrom, blk.start, blk.end, bed.strand, crange.size).intersect(crange)
    if coords is None:
        return None
    else:
        return genome_data.get_genome_seq(coords)

def _add_amplicaton_sequence(genome_data, target_transcript):
    amplicon_range = target_transcript.amplicon_range
    exon_seqs = []
    for blk in target_transcript.bed.blocks:
        seq = _get_exon_intersect_seq(genome_data, target_transcript.bed, blk, amplicon_range)
        if seq is not None:
            exon_seqs.append(seq)

    if target_transcript.strand == '-':
        exon_seqs = reversed(exon_seqs)
    target_transcript.amplicon = "".join(exon_seqs)

def _add_amplicaton_sequences(genome_data, target_transcripts):
    for target_transcript in target_transcripts:
        _add_amplicaton_sequence(genome_data, target_transcript)

def _do_target_transcripts_build(genome_data, primer_target_spec):
    target_transcripts = _target_transcripts_build(genome_data, primer_target_spec)

    _validate_strand(target_transcripts)
    region_5p, region_3p = _adjust_transcripts_features(target_transcripts)
    _add_amplicaton_sequences(genome_data, target_transcripts)

    return TargetTranscripts(primer_target_spec.target_id, region_5p, region_3p,
                             genome_data.get_genome_seq(region_5p),
                             genome_data.get_genome_seq(region_3p),
                             target_transcripts)

def target_transcripts_build(genome_data, primer_target_spec):
    """build TargetTranscripts object to a give primer and validity and consistency of
    the transcripts.
    """
    try:
        return _do_target_transcripts_build(genome_data, primer_target_spec)
    except PrimersJuJuDataError as ex:
        raise PrimersJuJuDataError(f"target {primer_target_spec.target_id} failed") from ex
