"""
Genome sequences and track data.
"""
import twobitreader
import pipettor
from pycbio.hgdata.bed import BedReader
from pycbio.hgdata import dnaOps
from . import PrimersJuJuDataError

class Track:
    """access annotations in a bigBed track by name.  Must be a bed12.  The
    bedBed maybe a file path or URL.  srcUrl is for doc and error message"""
    def __init__(self, track_name, bigBed, srcUrl):
        self.track_name = track_name
        self.bigBed = bigBed
        self.srcUrl = srcUrl

    def read_by_names(self, names):
        transcript_beds = {}
        with pipettor.Popen(['bigBedNamedItems', '-nameFile', self.bigBed, '/dev/stdin', '/dev/stdout'],
                            stdin=pipettor.DataWriter('\n'.join(names) + '\n')) as fh:
            for b in BedReader(fh):
                transcript_beds[b.name] = b

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

    def get_track(self, track_name):
        try:
            return self.tracks[track_name]
        except KeyError:
            raise PrimersJuJuDataError(f"unknown annotation track: '{track_name}'")
