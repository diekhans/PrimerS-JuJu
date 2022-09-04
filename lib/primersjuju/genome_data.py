"""
Genome sequences and track data.
"""
import os
from dataclasses import dataclass
from twobitreader import TwoBitFile
import pipettor
from pycbio.hgdata.bed import BedReader
from pycbio.hgdata import dnaOps
from pycbio.sys import fileOps
from pycbio.ncbi.assembly import AssemblyReport
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
    """Genome sequence and annotations tracks.  Files are opened in a lazy manner,
    so this can be specified in a configuration for multiple genomes with overhead
    of opening all of the files"""

    def __init__(self, genome_name, genome2bit, *, assembly_report=None):
        self.genome_name = genome_name
        self.genome2bit = genome2bit
        self.__genome_seqs = None  # lazy
        self.assembly_report = assembly_report
        self.__assembly_info = None  # lazy
        self.tracks = {}

    @property
    def genome_seqs(self) -> TwoBitFile:
        if self.__genome_seqs is None:
            self.__genome_seqs = TwoBitFile(self.genome2bit)
        return self.__genome_seqs

    @property
    def assembly_info(self) -> AssemblyReport:
        if (self.__assembly_info is None) and (self.assembly_report is not None):
            self.__assembly_info = AssemblyReport(self.assembly_report)
        return self.__assembly_info

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

def _bigbed_read_with_names(bigbed, names):
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
    return transcript_beds

def _bigbed_read_check_names(bigbed, names, transcript_beds):
    missing_names = set(names) - set(transcript_beds.keys())
    if len(missing_names) > 0:
        plural = "s" if len(missing_names) > 1 else ""
        raise PrimersJuJuDataError(f"record{plural} not found in bigBed {bigbed}: " + ", ".join(sorted(missing_names)))

def bigbed_read_by_names(bigbed, names):
    transcript_beds = _bigbed_read_with_names(bigbed, names)
    _bigbed_read_check_names(bigbed, names, transcript_beds)
    return transcript_beds

def bigbed_read_by_name(bigbed, name):
    return bigbed_read_by_names(bigbed, [name])[name]

def bigbed_fetch_by_name(bigbed, name):
    "returns None if not found"
    return _bigbed_read_with_names(bigbed, [name]).get(name)

def bigbed_read_by_range(bigbed, coords):
    "read by genomic range"
    with pipettor.Popen(['bigBedToBed',
                         '-chrom=' + coords.name, '-start=' + str(coords.start), '-end=' + str(coords.end),
                         bigbed, '/dev/stdout']) as fh:
        return [b for b in BedReader(fh)]
