"""
Genome sequences and track data.
"""
import os
from dataclasses import dataclass
import twobitreader
import pipettor
from pycbio.hgdata.bed import BedReader
from pycbio.hgdata import dnaOps
from pycbio.sys import fileOps
from . import PrimersJuJuDataError

@dataclass
class Track:
    """access annotations in a bigBed track by name.  Must be a bed12.  The
    bedBed maybe a file path or URL.  srcUrl is for doc and error message"""
    track_name: str
    bigbed: str
    src_url: str

    def read_by_name(self, name):
        try:
            return bigbed_read_by_name(self.bigbed, name)
        except Exception as ex:
            raise PrimersJuJuDataError(f"failed to read {name} from track {self.track_name}") from ex

    def read_by_names(self, names):
        try:
            return bigbed_read_by_names(self.bigbed, names)
        except Exception as ex:
            raise PrimersJuJuDataError(f"track {self.track_name}: {ex}") from ex

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

    def get_genome_seq(self, coords, *, strand='+'):
        """reverses coordinates if has '-' strand; with reverse-complement of sequence will
        be don't if strand is '-' """
        if coords.strand == '-':
            coords = coords.reverse()
        bases = self.genome_seqs[coords.name][coords.start:coords.end]
        if strand == '-':
            bases = dnaOps.reverseComplement(bases)
        return bases

    def get_chrom_size(self, chrom):
        return self.genome_seqs.sequence_sizes()[chrom]

    def get_track(self, track_name):
        try:
            return self.tracks[track_name]
        except KeyError:
            raise PrimersJuJuDataError(f"unknown annotation track: '{track_name}'")


def bigbed_read_by_names(bigbed, names):
    # FIXME:
    # stdin=pipettor.DataWriter('\n'.join(names) + '\n')
    # for some reason work on MacOS, but not Linux
    transcript_beds = {}

    tmpNamesFile = fileOps.tmpFileGet()
    fileOps.writeLines(tmpNamesFile, names)
    with pipettor.Popen(['bigBedNamedItems', '-nameFile', bigbed, tmpNamesFile, '/dev/stdout'],
                        stdin=tmpNamesFile) as fh:
        for b in BedReader(fh):
            transcript_beds[b.name] = b
    os.unlink(tmpNamesFile)

    missing_names = set(names) - set(transcript_beds.keys())
    if len(missing_names) > 0:
        plural = "s" if len(missing_names) > 1 else ""
        raise PrimersJuJuDataError(f"record{plural} not found in bigBed {bigbed}: " + ", ".join(sorted(missing_names)))
    return transcript_beds

def bigbed_read_by_name(bigbed, name):
    transcript_beds = bigbed_read_by_names(bigbed, [name])
    return transcript_beds[name]
