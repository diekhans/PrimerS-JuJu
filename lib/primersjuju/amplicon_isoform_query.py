"""
Find unspecified transcripts that will be amplified by an amplicon.
"""

from pycbio.hgdata.coords import Coords
from .genome_data import bigbed_read_by_range
from .transcript_features import TranscriptId, Transcript, bed_to_features

class AmpliconIsoformTracks:
    """
    Tracks to be check for isoforms to be amplified.
    For a given target, this will create an AmpliconIsoformQuery to
    avoid looking reloading the isoforms for each query.
    """
    def __init__(self, genome_data, track_names):
        self.genome_data = genome_data
        self.tracks = [genome_data.get_track(t) for t in track_names]

    def _bed_to_transcript(self, track, bed):
        return Transcript(TranscriptId(track.name, bed.name), bed,
                          bed_to_features(self.genome_data, bed))

    def _load_track_transcripts(self, gcoords, strand, track):
        return [self._bed_to_transcript(track, bed)
                for bed in bigbed_read_by_range(track.bigbed, gcoords)
                if bed.strand == strand]

    def get_amplicon_isoformat_query(self, gregion1, gregion2, strand):
        gcoords = Coords.adjrange(min(gregion1.start, gregion2.start),
                                  max(gregion1.end, gregion2.end))
        transcripts = []
        for track in self.tracks:
            transcripts.extend(self._load_track_transcripts(gcoords, strand, track))

class AmpliconIsoformQuery:
    """
    Query for amplicons for a given pair of target regiosn
    """
    def __init__(self, transcripts):
        self.transcripts = transcripts

    #def _is_amplified(self, primer_design, transcript):
