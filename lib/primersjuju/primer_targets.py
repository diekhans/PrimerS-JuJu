"""
Target transcripts analysis.  Includes validation, and trimming of primer
regions to match exons.
"""
from dataclasses import dataclass
import pprint
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuError, PrimersJuJuDataError
from .transcript_features import Feature, IntronFeature, ExonFeature, Features, bed_to_features, features_intersect_genome, get_features_rna

@dataclass
class TargetTranscript:
    "a target transcript for RT-PCR with features within primer regions"
    track_name: str
    bed: Bed
    features_5p: Features  # 5' on transcript
    features_3p: Features
    # features of transcript
    features: Features
    rna: str

    @property
    def trans_id(self):
        return self.bed.name

    @property
    def bounds(self) -> Feature:
        "bounds of transcript "
        return self.features.bounds

    @property
    def strand(self):
        return self.bed.strand

    @property
    def region_5p(self):
        return self.features_5p.bounds

    @property
    def region_3p(self):
        return self.features_3p.bounds

    @property
    def trans_len(self):
        return self.features_5p[0].trans.size

    def get_genome_ordered_features(self):
        "pair of features with 5' genome first"
        if self.strand == '+':
            return self.features_5p, self.features_3p
        else:
            return self.features_3p, self.features_5p

    def __str__(self):
        return f"({self.track_name}, {self.trans_id})"

    def dump(self, dump_fh):
        pp = pprint.PrettyPrinter(stream=dump_fh, sort_dicts=False, indent=4)
        print("transcript:", self.track_name, self.trans_id, file=dump_fh)
        print("coords:", str(self.bounds), file=dump_fh)
        print("region_5p:", self.region_5p, file=dump_fh)
        print("region_3p:", self.region_3p, file=dump_fh)
        print("features_5p:", file=dump_fh)
        pp.pprint(self.features_5p)
        print("features_3p:", file=dump_fh)
        pp.pprint(self.features_3p)
        print("features:", file=dump_fh)
        pp.pprint(self.features)
        print("rna:", self.rna, file=dump_fh)


@dataclass
class PrimerTargets:
    """Primer target region, with transcript, features in regions for each
    transcript.
    """
    target_id: str
    genome_name: str
    # Adjusted original genomic specification based on transcript feature, in transcription order.
    # The genomic strand is positive.
    region_5p: Coords
    region_3p: Coords
    strand: str
    transcripts: [TargetTranscript]

    def get_transcript(self, track_name, trans_id):
        for t in self.transcripts:
            if (t.track_name == track_name) and (t.bed.name == trans_id):
                return t
        raise PrimersJuJuDataError(f"({track_name}, {trans_id}) not found in {self.target_id}")

    def dump(self, dump_fh):
        print("target_id:", self.target_id, file=dump_fh)
        print("region_5p:", self.region_5p, file=dump_fh)
        print("region_3p:", self.region_3p, file=dump_fh)
        print("strand:", self.strand, file=dump_fh)
        print("transcripts:", file=dump_fh)
        for t in self.transcripts:
            t.dump(dump_fh)

def _primer_region_check_features(desc, track_name, trans_id, region, features):
    "check that features are sane, desc is used in error messages"
    exon_cnt = features.count_type(ExonFeature)
    intron_cnt = features.count_type(IntronFeature)
    if (exon_cnt + intron_cnt) == 0:
        raise PrimersJuJuDataError(f"{desc} {region} does not overlap transcript ({track_name}, {trans_id})")
    elif not (((exon_cnt == 1) and (intron_cnt == 0)) or ((exon_cnt == 2) and (intron_cnt == 1))):
        raise PrimersJuJuDataError(f"{desc} {region} for transcript ({track_name}, {trans_id}) must contain either one exon, or two exons and an intron: {str(features)}")
    if not (isinstance(features[0], ExonFeature) and
            isinstance(features[-1], ExonFeature)):
        raise PrimersJuJuDataError(f"{desc} transcript ({track_name}, {trans_id}) primer region does not end in exons: {str(features)}")


def _get_regions_transcript_orient(trans, region_5p, region_3p):
    "swap regions if needed to be in transcript orientation"
    genome_orient = '+' if (region_5p.end < region_3p.start) else '-'
    if trans.strand != genome_orient:
        region_5p, region_3p = region_3p, region_5p
    return region_5p, region_3p

def _build_region_transcript_features(track_name, trans_name, features, region):
    region_features = features_intersect_genome(features, region)
    _primer_region_check_features("specified primer region", track_name, trans_name, region, region_features)
    return region_features

def _build_target_transcript(genome_data, primer_target_spec, trans_spec):
    "build transcript with initial regions trimmed to exons"
    trans_bed = genome_data.get_track(trans_spec.trans_track).read_by_name(trans_spec.trans_id)
    features = bed_to_features(genome_data, trans_bed)
    rna = get_features_rna(genome_data, features)
    region_5p, region_3p = _get_regions_transcript_orient(trans_bed, primer_target_spec.region_5p,
                                                          primer_target_spec.region_3p)
    return TargetTranscript(trans_spec.trans_track, trans_bed,
                            _build_region_transcript_features(trans_spec.trans_track, trans_spec.trans_id, features, region_5p),
                            _build_region_transcript_features(trans_spec.trans_track, trans_spec.trans_id, features, region_3p),
                            features, rna)

def transcripts_build(genome_data, primer_target_spec):
    transcripts = []
    for track in primer_target_spec.tracks.values():
        for trans_spec in track.values():
            transcripts.append(_build_target_transcript(genome_data, primer_target_spec, trans_spec))
    return transcripts

def _validate_strand(transcripts):
    # must do before other checks, as primer region swap could confuse error messages
    ttrans0 = transcripts[0]
    for ttrans in transcripts[1:]:
        if ttrans.region_5p.strand != ttrans0.strand:
            raise PrimersJuJuDataError(f"transcript on different strands: {ttrans} vs {ttrans0}")

def _find_transcripts_common_region(transcripts, feats_func):
    region = feats_func(transcripts[0]).bounds.genome
    for ttrans in transcripts[1:]:
        region = region.intersect(feats_func(ttrans).bounds.genome)
    return region

def _adjust_transcript_region_features(target_transcript, common_region, feats_func):
    orig_features = feats_func(target_transcript)
    adj_features = Features()
    for ofeat in orig_features:
        if ofeat.genome.overlaps(common_region):
            adj_features.append(ofeat.intersect_genome(common_region))

    _primer_region_check_features("common adjusted primer region", target_transcript.track_name, target_transcript.trans_id, common_region, adj_features)

    # this updates transcript features
    feats_func(target_transcript, adj_features)

def _adjust_transcripts_region_features(transcripts, feats_func):
    common_region = _find_transcripts_common_region(transcripts, feats_func)
    for target_transcript in transcripts:
        _adjust_transcript_region_features(target_transcript, common_region, feats_func)
    return common_region

def _adjust_transcripts_features(transcripts):
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

    region_5p = _adjust_transcripts_region_features(transcripts, features_5p_access)
    region_3p = _adjust_transcripts_region_features(transcripts, features_3p_access)
    return region_5p, region_3p

def _do_primer_targets_build(genome_data, primer_target_spec):
    transcripts = transcripts_build(genome_data, primer_target_spec)

    _validate_strand(transcripts)
    region_5p, region_3p = _adjust_transcripts_features(transcripts)

    return PrimerTargets(primer_target_spec.target_id, genome_data.genome_name, region_5p, region_3p,
                         transcripts[0].strand, transcripts)

def primer_targets_build(genome_data, primer_target_spec):
    """build PrimerTargets object to a give primer and validity and consistency of
    the transcripts.
    """
    try:
        return _do_primer_targets_build(genome_data, primer_target_spec)
    except PrimersJuJuDataError as ex:
        raise PrimersJuJuDataError(f"target {primer_target_spec.target_id} failed") from ex
    except PrimersJuJuError as ex:
        raise PrimersJuJuError(f"target {primer_target_spec.target_id} failed") from ex
