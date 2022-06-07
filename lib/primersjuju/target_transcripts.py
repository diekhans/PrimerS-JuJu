"""
Target transcripts analysis.  Includes validation, and trimming of primer
regions to match exons.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuError, PrimersJuJuDataError
from .transcript_features import Feature, IntronFeature, ExonFeature, Features, bed_to_features, features_intersect_genome, get_features_rna

@dataclass
class PrimerRegionFeatures:
    """Coords and features of a transcript for 5' or 3' region.  """
    region: Feature  # genomic and trans may not be the same length
    features: Sequence[Feature]

@dataclass
class TargetTranscript:
    "a target transcript for RT-PCR with features within primer regions"
    track_name: str
    bed: Bed
    features_5p: PrimerRegionFeatures
    features_3p: PrimerRegionFeatures
    # features of transcript
    features: Features
    rna: str

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
    def amplicon_genome_range(self) -> Coords:
        genome_5p = self.region_5p.genome
        genome_3p = self.region_3p.genome
        if self.strand == '+':
            return Coords(genome_5p.name,
                          genome_3p.start, genome_5p.end,
                          genome_5p.strand, genome_5p.size)
        else:
            return Coords(genome_5p.name,
                          genome_5p.start, genome_3p.end,
                          genome_5p.strand, genome_5p.size)

    @property
    def amplicon_trans_range(self) -> Coords:
        "coordinates on transcript (reverse of primers)"
        trans_5p = self.region_5p.trans
        trans_3p = self.region_3p.trans
        return Coords(trans_5p.name,
                      trans_3p.start, trans_5p.end,
                      trans_5p.strand, trans_5p.size)

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

def _primer_region_check_features(desc, track_name, trans_id, features):
    "check that features are sane, desc is used in error messages"
    exon_cnt = features.count_type(ExonFeature)
    intron_cnt = features.count_type(IntronFeature)
    if not (((exon_cnt == 1) and (intron_cnt == 0)) or ((exon_cnt == 2) and (intron_cnt == 1))):
        raise PrimersJuJuDataError(f"{desc} for transcript ({track_name}, {trans_id}) must contain either one exon, or two exons and an intron: {str(features)}")
    if not (isinstance(features[0], ExonFeature) and
            isinstance(features[-1], ExonFeature)):
        raise PrimersJuJuDataError(f"{desc} transcript ({track_name}, {trans_id}) primer region does not end in exons: {features}")


def _block_features(trans_bed, trans_off, trans_size, region, genome_size, prev_blk, blk, features):
    """get intron and exon feature intersection with a genomic range.  trans_off is in the genomic
    direction, not direction of transcription"""

    def _mk_genome(start, end):
        "create genomic range for intersecting the region"
        return Coords(trans_bed.chrom,
                      max(start, region.start), min(end, region.end),
                      strand=trans_bed.strand, size=genome_size)

    def _mk_trans(start, end):
        "create transcript range for intersecting the region"
        trans = Coords(trans_bed.name, start, end,
                       strand=trans_bed.strand, size=trans_size)
        if trans.strand == '-':
            trans = trans.reverse()
        return trans

    if (prev_blk is not None) and (prev_blk.end < region.end) and (blk.start > region.start):
        features.append(IntronFeature(_mk_genome(prev_blk.end, blk.start),
                                      _mk_trans(trans_off, trans_off)))
    if (blk.start < region.end) and (blk.end > region.start):
        genome = _mk_genome(blk.start, blk.end)
        trans_start = trans_off + (genome.start - blk.start)
        trans = _mk_trans(trans_start,
                          trans_start + len(genome))
        features.append(ExonFeature(genome, trans))

def _get_regions_primer_orient(trans, region_5p, region_3p):
    "swap regions if needed to be in primer orientation (reverse of transcript strand)"
    genome_orient = '+' if (region_5p.end < region_3p.start) else '-'
    if trans.strand == genome_orient:
        region_5p, region_3p = region_3p, region_5p
    return region_5p, region_3p

def _build_region_transcript_features(track_name, trans_name, features, region):
    region_features = features_intersect_genome(features, region)
    _primer_region_check_features("initially specified primer region", track_name, trans_name, region_features)
    return PrimerRegionFeatures(region_features.get_bounds(), region_features)

def _build_target_transcript(genome_data, primer_target_spec, trans_spec):
    "build transcript with initial regions trimmed to exons"
    trans_bed = genome_data.get_track(trans_spec.trans_track).read_by_name(trans_spec.trans_id)
    features = bed_to_features(genome_data, trans_bed)
    rna = get_features_rna(genome_data, features)
    region_5p, region_3p = _get_regions_primer_orient(trans_bed, primer_target_spec.region_5p,
                                                      primer_target_spec.region_3p)
    return TargetTranscript(trans_spec.trans_track, trans_bed,
                            _build_region_transcript_features(trans_spec.trans_track, trans_spec.trans_id, features, region_5p),
                            _build_region_transcript_features(trans_spec.trans_track, trans_spec.trans_id, features, region_3p),
                            features, rna)

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
    region = feats_func(target_transcripts[0]).region.genome
    for ttrans in target_transcripts[1:]:
        region = region.intersect(feats_func(ttrans).region.genome)
    return region

def _adjust_transcript_region_features(target_transcript, common_region, feats_func):
    orig_features = feats_func(target_transcript).features
    adj_features = Features()
    for ofeat in orig_features:
        if ofeat.genome.overlaps(common_region):
            adj_features.append(ofeat.intersect_genome(common_region))

    _primer_region_check_features("common adjusted primer region", target_transcript.track_name, target_transcript.trans_id, adj_features)

    # this updates transcript features
    feats_func(target_transcript,
               PrimerRegionFeatures(adj_features.get_bounds(),
                                    adj_features))

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

def _do_target_transcripts_build(genome_data, primer_target_spec):
    target_transcripts = _target_transcripts_build(genome_data, primer_target_spec)

    _validate_strand(target_transcripts)
    region_5p, region_3p = _adjust_transcripts_features(target_transcripts)

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
    except PrimersJuJuError as ex:
        raise PrimersJuJuError(f"target {primer_target_spec.target_id} failed") from ex


def _get_regions_genome_orient(target_transcripts):
    "swap regions if needed to be in genome orientation"
    if target_transcripts.region_5p.end < target_transcripts.region_3p.start:
        return target_transcripts.region_5p, target_transcripts.region_3p
    else:
        return target_transcripts.region_3p, target_transcripts.region_5p

def build_target_bed(target_transcripts, color):
    region_5p, region_3p = _get_regions_genome_orient(target_transcripts)
    bed = Bed(region_5p.name, region_5p.start, region_3p.end,
              target_transcripts.target_id,
              strand=region_5p.strand,
              thickStart=region_5p.start, thickEnd=region_3p.end,
              itemRgb=color.toRgb8Str(),
              blocks=[Bed.Block(region_5p.start, region_5p.end),
                      Bed.Block(region_3p.start, region_3p.end)])
    return bed
