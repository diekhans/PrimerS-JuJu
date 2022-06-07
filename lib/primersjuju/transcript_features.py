"""
Target transcripts features.  Loads BED into list of genomic and transcript feature coordinates
and handles conversions of sub0ranges.
"""
from collections import namedtuple
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuError

class Feature(namedtuple("Feature", ("genome", "trans"))):
    """annotation feature, both genome and transcript coordinates (for Exons)"""

    def intersect_genome(self, other):
        "intersect with genomic coordinates, None if no intersection"
        if not isinstance(other, Coords):
            raise ValueError(f"bad object type: {type(other)}, expected {Coords}")
        if other.name != self.genome.name:
            raise ValueError(f"mismatch genome sequence name '{other.name}', expected '{self.genome.name}'")
        genome_intr = self.genome.intersect(other)
        if len(genome_intr) == 0:
            return None
        elif isinstance(self, ExonFeature):
            assert other.contains(genome_intr), f"{other}.contains({genome_intr})"
            assert self.genome.contains(genome_intr), f"{self.genome}.contains({genome_intr})"
            genome_off = (genome_intr.start - self.genome.start)
            if self.trans.strand == self.genome.strand:
                trans_start = self.trans.start + genome_off
                trans_end = trans_start + len(genome_intr)
            else:
                trans_end = self.trans.end - genome_off
                trans_start = trans_end - len(genome_intr)

            assert trans_start < trans_end
            trans_intr = Coords(self.trans.name, trans_start, trans_end,
                                self.trans.strand, self.trans.size)
            assert len(trans_intr) == len(genome_intr)
            return ExonFeature(genome_intr, trans_intr)
        elif isinstance(self, IntronFeature):
            assert len(self.trans) == 0
            return IntronFeature(genome_intr, self.trans)
        else:
            raise PrimersJuJuError("intersect_genome not support on base Feature class")

    def intersect_transcript(self, other):
        "intersect with transcript coordinates, None if no intersection"
        if not isinstance(other, Coords):
            raise ValueError(f"bad object type: {type(other)}, expected {Coords}")
        if other.name != self.trans.name:
            raise ValueError(f"mismatch transcript name '{other.name}', expected '{self.trans.name}'")
        if isinstance(self, IntronFeature):
            return None
        trans_intr = self.trans.intersect(other)
        if len(trans_intr) == 0:
            return None
        if isinstance(self, ExonFeature):
            assert other.contains(trans_intr), f"{other}.contains({trans_intr})"
            assert self.trans.contains(trans_intr), f"{self.trans}.contains({trans_intr})"
            trans_off = (trans_intr.start - self.trans.start)
            if self.genome.strand == self.trans.strand:
                genome_start = self.genome.start + trans_off
                genome_end = genome_start + len(trans_intr)
            else:
                genome_end = self.genome.end - trans_off
                genome_start = genome_end - len(trans_intr)

            assert genome_start < genome_end
            genome_intr = Coords(self.genome.name, genome_start, genome_end,
                                self.genome.strand, self.genome.size)
            assert len(genome_intr) == len(trans_intr)
            return ExonFeature(genome_intr, trans_intr)
        else:
            raise PrimersJuJuError("intersect_transcript not support on base Feature class")

class ExonFeature(Feature):
    "exon in a model, with genome and trans coordinates "
    pass


class IntronFeature(Feature):
    "intron in a model with zero length transcript coordinates of the intron location"
    pass


class Features(list):
    "a list of transcript intron/exon features or a sub-range of those features"

    def __str__(self):
        return '[' + ("\n ".join([repr(f) for f in self])) + ']'

    def sort_genome_order(self):
        self.sort(key=lambda f: (f.genome.start, f.genome.end))

    def sort_transcript_order(self):
        self.sort(key=lambda f: (f.trans.start, f.trans.end))

    def count_type(self, ftype):
        # True == 1, False == 0
        return sum([isinstance(f, ftype) for f in self])

    def get_bounds(self):
        f0 = self[0]
        fN = self[-1]
        return Feature(f0.genome.adjrange(f0.genome.start, fN.genome.end),
                       f0.trans.adjrange(f0.trans.start, fN.trans.end))

def _bed_block_features(trans_bed, trans_off, trans_size, genome_size, prev_blk, blk, features):
    """get intron and exon features for a BED block.  trans_off is in the genomic
    direction, not direction of transcription
    """

    def _mk_genome(start, end):
        "create genomic range for intersecting the region"
        return Coords(trans_bed.chrom, start, end,
                      strand=trans_bed.strand, size=genome_size)

    def _mk_trans(start, end):
        "create transcript range for intersecting the region"
        trans = Coords(trans_bed.name, start, end,
                       strand=trans_bed.strand, size=trans_size)
        if trans.strand == '-':
            trans = trans.reverse()
        return trans

    if prev_blk is not None:
        features.append(IntronFeature(_mk_genome(prev_blk.end, blk.start),
                                      _mk_trans(trans_off, trans_off)))
    genome = _mk_genome(blk.start, blk.end)
    trans_start = trans_off + (genome.start - blk.start)
    trans = _mk_trans(trans_start,
                      trans_start + len(genome))
    features.append(ExonFeature(genome, trans))

def bed_to_features(genome_data, trans_bed) -> Features:
    """Generate features list from a transcript BED
    """
    genome_size = genome_data.get_chrom_size(trans_bed.chrom)
    trans_off = 0
    trans_size = trans_bed.coverage

    features = Features()
    prev_blk = None
    for blk in trans_bed.blocks:
        _bed_block_features(trans_bed, trans_off, trans_size, genome_size, prev_blk, blk, features)
        trans_off += len(blk)
        prev_blk = blk
    features_contig_assert(features)
    return features

def features_contig_assert(features):
    "are features contiguous in both genome and transcriptiome?"

    def _contig(prev, cur):
        "check for contiguous in either sort order"
        return (cur.start == prev.end) or (cur.end == prev.start)

    # this loop structures handles zero and one length lists
    prev_feat = None
    for feat in features:
        if prev_feat is not None:
            assert _contig(feat.trans, prev_feat.trans), f"trans coords not contiguous {prev_feat} => {feat}"
            assert _contig(feat.genome, prev_feat.genome), f"genome coords not contiguous {prev_feat} => {feat}"
        prev_feat = feat

def features_intersect_genome(features, region) -> Features:
    """get a sub-range of features based on genomic coordinates, sorted in genomic order"""
    subfeatures = Features([])
    for feat in features:
        subfeat = feat.intersect_genome(region)
        if subfeat is not None:
            subfeatures.append(subfeat)
    features_contig_assert(subfeatures)
    return subfeatures

def get_features_rna(genome_data, features):
    exon_seqs = []
    for feat in features:
        if isinstance(feat, ExonFeature):
            exon_seqs.append(genome_data.get_genome_seq(feat.genome))
    if features[0].genome.strand == '-':
        exon_seqs = reversed(exon_seqs)
    return "".join(exon_seqs)

def transcript_range_to_features(features, trange):
    """convert transcript coordinates to one or more Feature coordinates"""
    exon_regions = Features()
    for feat in features:
        if isinstance(feat, ExonFeature):
            exon_region = feat.intersect_transcript(trange)
            if exon_region is not None:
                assert isinstance(exon_region, ExonFeature)
                exon_regions.append(exon_region)
    return exon_regions

def genome_range_to_featyres(features, grange):
    """convert transcript coordinates to one or more Feature coordinates"""
    exon_regions = Features()
    for feat in features:
        if isinstance(feat, ExonFeature):
            exon_region = feat.intersect_genome(grange)
            if exon_region is not None:
                assert isinstance(exon_region, ExonFeature)
                exon_regions.append(exon_region)
    return exon_regions
