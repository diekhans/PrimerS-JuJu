"""
Genome sequences and track data.
"""
from dataclasses import dataclass
import twobitreader
import pipettor
from pycbio.hgdata.bed import Bed, BedReader
from pycbio.hgdata import dnaOps
from . import PrimersJuJuDataError

class Transcript(Bed):
    "transcript is a BED record with extra attributes"
    __slots__ = ("track_name")

@dataclass
class Track:
    """access annotations in a bigBed track by name.  Must be a bed12.  The
    bedBed maybe a file path or URL.  srcUrl is for doc and error message"""
    track_name: str
    big_bed: str
    src_url: str

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

    def read_by_name(self, name):
        transcript_beds = self.read_by_names([name])
        return transcript_beds[0]


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
