"""
Target transcripts analysis.  Includes validation, and trimming of primer
regions to match exons.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuDataError
from .genome_data import Transcript


class Feature(Coords):
    "annotation feature"
    pass

class ExonRegion(Feature):
    "exon in a model"
    pass

class IntronRegion(Feature):
    "intron in a model"
    pass


@dataclass
class TranscriptRegions:
    "primer regions and features for a given transcript"
    track_name: str
    trans: Transcript  # this is a BED record
    region_5p: Coords
    region_3p: Coords
    features_5p: Sequence[Feature]
    features_3p: Sequence[Feature]


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

def _build_region_features(genome_data, trans, region):
    features = get_transcript_region_features(genome_data, trans, region)
    exon_cnt = filter(lambda f: isinstance(f, ExonRegion), features)
    if exon_cnt == 0:
        raise PrimersJuJuDataError(f"transcript {trans.track_name} {trans.id} has no exons in region {region}")
    return features

def _build_transcript_regions(genome_data, primer_target_spec, trans_spec):
    "build transcript with initial regions trimmed to exons"
    trans = genome_data.get_track(trans_spec.track_name).read_by_name(trans_spec.trans_id)
    region_5p, region_3p = _get_regions_strand(trans, primer_target_spec.region_5p,
                                               primer_target_spec.region_3p)
    return TranscriptRegions(trans, region_5p, region_3p,
                             _build_region_features(genome_data, trans, region_5p),
                             _build_region_features(genome_data, trans, region_3p))
