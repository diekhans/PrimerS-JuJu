"""
Target transcripts features.  Loads BED into list of genomic and transcript feature coordinates
and handles conversions of sub0ranges.
"""
from dataclasses import dataclass, KW_ONLY
import pprint
from collections import namedtuple
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuError

class Feature(namedtuple("Feature", ("genome", "trans"))):
    """annotation feature, both genome and transcript coordinates (for Exons)"""

    def __str__(self):
        return f"{self.__class__.__name__}(genome={str(self.genome)}/{self.genome.strand}, trans={str(self.trans)}/{self.trans.strand})"

    def intersect_genome(self, other):
        """intersect with genomic coordinates, None if no intersection,
        if other is on the opposite strand, switch other to match."""
        if not isinstance(other, Coords):
            raise ValueError(f"bad object type: {type(other)}, expected {Coords}")
        if other.name != self.genome.name:
            raise ValueError(f"mismatch genome sequence name '{other.name}', expected '{self.genome.name}'")
        if other.strand != self.genome.strand:
            other = other.reverse()
        genome_intr = self.genome.intersect(other)
        if len(genome_intr) == 0:
            return None
        elif isinstance(self, ExonFeature):
            assert other.contains(genome_intr), f"{other}.contains({genome_intr})"
            assert self.genome.contains(genome_intr), f"{self.genome}.contains({genome_intr})"
            genome_off = (genome_intr.start - self.genome.start)
            trans_start = self.trans.start + genome_off
            trans_end = trans_start + len(genome_intr)
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
        """intersect with transcript coordinates, None if no intersection,
        if other is on the opposite strand, switch other to match."""
        if not isinstance(other, Coords):
            raise ValueError(f"bad object type: {type(other)}, expected {Coords}")
        if other.name != self.trans.name:
            raise ValueError(f"mismatch transcript name '{other.name}', expected '{self.trans.name}'")
        if isinstance(self, IntronFeature):
            return None
        if other.strand != self.trans.strand:
            other = other.reverse()
        trans_intr = self.trans.intersect(other)
        if len(trans_intr) == 0:
            return None
        if isinstance(self, ExonFeature):
            assert other.contains(trans_intr), f"{other}.contains({trans_intr})"
            assert self.trans.contains(trans_intr), f"{self.trans}.contains({trans_intr})"
            trans_off = (trans_intr.start - self.trans.start)
            genome_start = self.genome.start + trans_off
            genome_end = genome_start + len(trans_intr)
            assert genome_start < genome_end
            genome_intr = Coords(self.genome.name, genome_start, genome_end,
                                 self.genome.strand, self.genome.size)
            assert len(genome_intr) == len(trans_intr)
            return ExonFeature(genome_intr, trans_intr)
        else:
            raise PrimersJuJuError("intersect_transcript not support on base Feature class")

    def reverse(self):
        """strand-reverse all features, creating instances of derived class"""
        return type(self)(self.genome.reverse(), self.trans.reverse())


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

    def count_type(self, ftype):
        # True == 1, False == 0
        return sum([isinstance(f, ftype) for f in self])

    def iter_type(self, ftype):
        "generator over a feature type"
        for f in self:
            if isinstance(f, ftype):
                yield f

    # FIXME: weird names
    def genome_coords_type(self, ftype):
        return [f.genome for f in self.iter_type(ftype)]

    @property
    def bounds(self):
        f0 = self[0]
        fN = self[-1]
        return Feature(f0.genome.adjrange(f0.genome.start, fN.genome.end),
                       f0.trans.adjrange(f0.trans.start, fN.trans.end))

    def reverse(self):
        """new Features object with strand-reverse coordinates and reversed order"""
        return Features([f.reverse() for f in reversed(self)])

class TranscriptId(namedtuple("TranscriptId", ("track", "name"))):
    "uniquely identifies a transcript by track and name"
    __slots__ = ()

    def __str__(self):
        return f"{self.track}/{self.name}"

@dataclass
class Transcript:
    "A transcript with features and optional RNA sequence"
    _: KW_ONLY
    trans_id: TranscriptId
    bed: Bed
    features: Features
    rna: str = None

    def __str__(self):
        return str(self.trans_id)

    @property
    def bounds(self) -> Feature:
        "bounds of transcript "
        return self.features.bounds

    @property
    def strand(self):
        return self.bed.strand

    @property
    def trans_len(self):
        return self.features[0].trans.size

    def dump(self, dump_fh):
        pp = pprint.PrettyPrinter(stream=dump_fh, sort_dicts=False, indent=4)
        print("transcript:", self.trans_id, file=dump_fh)
        print("coords:", str(self.bounds), file=dump_fh)
        print("features:", file=dump_fh)
        pp.pprint(self.features)
        print("rna:", self.rna, file=dump_fh)
        print("rna exons:", ','.join(self._features_to_trans_seqs(self.features)), file=dump_fh)


def _bed_block_features(trans_bed, trans_start, trans_size, genome_size, prev_blk, blk, features):
    """get intron and exon features for a BED block.  trans_start is in the
    genomic direction, not the direction of transcription.  """

    def _mk_genome(start, end):
        "create genomic range for a region"
        return Coords(trans_bed.chrom, start, end,
                      strand='+', size=genome_size)

    def _mk_trans(start, end):
        "create transcript range for a region"
        trans = Coords(trans_bed.name, start, end,
                       strand=trans_bed.strand, size=trans_size)
        return trans

    if prev_blk is not None:
        features.append(IntronFeature(_mk_genome(prev_blk.end, blk.start),
                                      _mk_trans(trans_start, trans_start)))
    genome = _mk_genome(blk.start, blk.end)
    trans = _mk_trans(trans_start, trans_start + len(genome))
    features.append(ExonFeature(genome, trans))

def bed_to_features(genome_data, trans_bed) -> Features:
    """Generate features list from a transcript BED
    """
    genome_size = genome_data.get_chrom_size(trans_bed.chrom)
    trans_size = trans_bed.coverage
    trans_start = 0

    features = Features()
    prev_blk = None
    for blk in trans_bed.blocks:
        _bed_block_features(trans_bed, trans_start, trans_size, genome_size, prev_blk, blk, features)
        trans_start += len(blk)
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

def features_sort_genome(features: Features):
    features.sort(key=lambda f: (f.genome.start, f.genome.end))

def features_sort_transcript(features: Features):
    features.sort(key=lambda f: (f.trans.start, f.trans.end))

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
            exon_seqs.append(genome_data.get_genome_seq(feat.genome, strand=feat.trans.strand))
    if features[0].trans.strand == '-':
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

def genome_range_to_features(features, grange):
    """convert transcript coordinates to one or more Feature coordinates"""
    exon_regions = Features()
    for feat in features:
        if isinstance(feat, ExonFeature):
            exon_region = feat.intersect_genome(grange)
            if exon_region is not None:
                assert isinstance(exon_region, ExonFeature)
                exon_regions.append(exon_region)
    return exon_regions

def features_to_genomic_coords(features, feature_filter=Feature):
    """create list of genomic coordinates, only keep features of time feature_filter"""
    return [f.genome for f in features if isinstance(f, feature_filter)]
