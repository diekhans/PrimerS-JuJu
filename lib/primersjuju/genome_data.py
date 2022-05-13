"""
Genome sequences and track data.
"""
import twobitreader
import pipettor
from pycbio.hgdata.bed import Bed, BedReader
from pycbio.hgdata import dnaOps
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuDataError

class Transcript(Bed):
    "transcript is a BED record with extra attributes"
    __slots__ = ("track_name")

class Track:
    """access annotations in a bigBed track by name.  Must be a bed12.  The
    bedBed maybe a file path or URL.  srcUrl is for doc and error message"""
    def __init__(self, track_name, big_bed, src_url):
        self.track_name = track_name
        self.big_bed = big_bed
        self.src_url = src_url

    def read_by_names(self, names):
        transcript_beds = {}
        with pipettor.Popen(['bigBedNamedItems', '-nameFile', self.big_bed, '/dev/stdin', '/dev/stdout'],
                            stdin=pipettor.DataWriter('\n'.join(names) + '\n')) as fh:
            for b in BedReader(fh, bedClass=Transcript):
                transcript_beds[b.name] = b
                b.track_name = self.track_name

        missing_names = set(names) - set(transcript_beds.keys())
        if len(missing_names) > 0:
            raise PrimersJuJuDataError(f"{len(missing_names)} record(s) not found in track {self.track_name}: " + ", ".join(sorted(missing_names)))
        return transcript_beds

class GenomeData:
    "genome sequence and annotations tracks"

    def __init__(self, genome_name, genome2bit, genome_url):
        self.genome_name = genome_name
        self.genome2bit = genome2bit
        self.genome_seqs = twobitreader.TwoBitFile(genome2bit)
        self.genome_url = genome_url
        self.tracks = {}

    def add_track(self, track_name, bigBed, srcUrl):
        self.tracks[track_name] = Track(track_name, bigBed, srcUrl)

    def get_genome_seq(self, coords):
        "reverse-complements if coords has '-' strand"
        bases = self.genome_seqs[coords.name][coords.start:coords.end]
        if coords.strand == '-':
            bases = dnaOps.reverseComplement(bases)
        return bases

    def get_chrom_size(self, chrom):
        return self.genome_seqs.sequence_sizes()[chrom]

    def get_track(self, track_name):
        try:
            return self.tracks[track_name]
        except KeyError:
            raise PrimersJuJuDataError(f"unknown annotation track: '{track_name}'")


class Feature(Coords):
    "annotation feature"
    pass

class ExonRegion(Feature):
    "exon in a model"
    pass

class IntronRegion(Feature):
    "intron in a model"
    pass


def _block_features(trans, crange, csize, prev_blk, blk, features):
    "get intron and exon feature intersection with range"

    def _mk_feature(feat_cls, start, end):
        "create Feature for intersecting the range"
        return feat_cls(trans.chrom,
                        max(start, crange.start), min(end, crange.end),
                        strand=trans.strand, size=csize)

    if prev_blk is not None:
        if (prev_blk.end < crange.end) and (blk.start > crange.start):
            features.append(_mk_feature(IntronRegion, prev_blk.end, blk.start))
    if (blk.start < crange.end) and (blk.end > crange.start):
        features.append(_mk_feature(ExonRegion, blk.start, blk.end))

def get_transcript_crange_features(genome_data, trans, crange):
    """Given a chromosome range in a transcript, generate of a list feature
    coords in that range.
    """
    csize = genome_data.get_chrom_size(trans.chrom)

    features = []
    prev_blk = None
    for blk in trans.blocks:
        _block_features(trans, crange, csize, prev_blk, blk, features)
        prev_blk = blk
    return features
