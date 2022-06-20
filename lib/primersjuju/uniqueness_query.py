"""
Query for uniqueness in genome and transcriptome
"""
from typing import Sequence
from dataclasses import dataclass, KW_ONLY

import pipettor
from pycbio.hgdata.psl import Psl, PslReader
from pycbio.hgdata.coords import Coords

from .genome_data import bigbed_read_by_names
from .transcript_features import Features, bed_to_features, transcript_range_to_features
from . import PrimersJuJuError

def _coords_range(coords_list):
    "convert a list of coordinates to the range that it spans"
    coords_list = sorted(coords_list)
    c0 = coords_list[0]
    return Coords(c0.name, coords_list[0].start, coords_list[-1].end, c0.strand, c0.size)

@dataclass
class IsPcrServerSpec:
    """Specification of an isPCR server, static or dynamic;
    for genome or transcriptome"""
    host: str
    port: str
    target_seq_dir: str  # contains 2bit matching one returned by server
    _: KW_ONLY
    dyn_name: str = None  # target or transcriptome name for dynamic server
    dyn_data_dir: str = None   # dynamic blat data dir
    trans_bigbed: str = None  # big bed file/URL for transcriptome

@dataclass
class GenomeHit:
    "one isPcr hit to the genome, in positive genomic coords and order"
    left_coords: Coords
    right_coords: Coords
    alignment: Psl

    def __str__(self):
        return f"{self.__class__.__name__}(left={str(self.left_coords)},right={str(self.right_coords)})"

    def get_genome_range(self):
        """get the range covering both alignments"""
        return _coords_range([self.left_coords, self.right_coords])

@dataclass
class TranscriptomeHit:
    "one isPcr hit to a transcript, with mappings back to genome, in positive genomic coordinates and order"
    trans_id: str
    gene_name: str
    left_features: Features
    right_features: Features
    alignment: Psl

    def __str__(self):
        def _lfmt(l):
            return '[' + ",\n\t".join([str(v) for v in l]) + '\n    ]'
        return f"{self.__class__.__name__}(left={_lfmt(self.left_features)},right={_lfmt(self.right_features)})"

    def get_genome_range(self):
        """get the range covering both  """
        return _coords_range([self.left_features.bounds.genome, self.right_features.bounds.genome])

def _gfPcr(spec, name, left_primer, right_primer, max_size):
    "returns PSL records"
    cmd = ["gfPcr", f"-maxSize={max_size}", "-out=psl", f"-name={name}"]
    if spec.dyn_name is not None:
        cmd.append(f"-genome={spec.dyn_name}")
    if spec.dyn_data_dir is not None:
        cmd.extend(f"-genomeDataDir={spec.dyn_data_dir}")
    cmd.extend([spec.host, spec.port, spec.target_seq_dir, left_primer, right_primer, "/dev/stdout"])
    with pipettor.Popen(cmd) as fh:
        return [p for p in PslReader(fh)]

def _check_psl(psl):
    if len(psl.blocks) != 2:
        raise PrimersJuJuError(f"expected a two-block result back from isPcr, got: {psl}")

def _genome_psl_to_hit(psl):
    """create hit records in positive genomic coordinates"""
    _check_psl(psl)
    coords = [Coords(psl.tName, psl.blocks[i].tStart, psl.blocks[i].tEnd, psl.tStrand, psl.tSize).abs()
              for i in range(2)]
    return GenomeHit(*coords, psl)

def _split_transcriptome_id(ispcr_trans_id):
    """split id in the form ENST00000244050.3__SNAI1, second part is optional"""
    return ispcr_trans_id.split('__')

def _trans_to_features(genome_data, transcriptome_spec, trans_pcr_psls):
    """alignments of transcripts to the genome as features"""
    trans_ids = set([_split_transcriptome_id(p.tName)[0] for p in trans_pcr_psls])
    return {b.name: bed_to_features(genome_data, b)
            for b in bigbed_read_by_names(transcriptome_spec.trans_bigbed, trans_ids).values()}

def _trans_range_to_features(trans_features, coords):
    """map a transcript range to features, with positive genome coordinates"""
    features = transcript_range_to_features(trans_features, coords)
    if features[0].genome.strand == '-':
        features = features.reverse()
    return features

def _trans_psl_to_hit(trans_to_features, psl):
    _check_psl(psl)
    trans_id, gene_name = _split_transcriptome_id(psl.tName)
    trans_features = trans_to_features[trans_id]
    # transcript coordinates
    coords_list = [Coords(trans_id, psl.blocks[i].tStart, psl.blocks[i].tEnd, psl.tStrand, psl.tSize)
                   for i in range(2)]
    # multiple features per primer if crosses splice sites
    features_list = [_trans_range_to_features(trans_features, coords)
                     for coords in coords_list]
    if features_list[0][0].genome.start > features_list[1][0].genome.start:
        features_list.reverse()   # put in genomic order
    return TranscriptomeHit(trans_id, gene_name, *features_list, psl)


class UniquenessQuery:
    """Interface to UCSC isPCR server to query for uniqueness."""
    def __init__(self, genome_data, genome_spec, transcriptome_spec):
        assert transcriptome_spec.trans_bigbed is not None
        self.genome_data = genome_data
        self.genome_spec = genome_spec
        self.transcriptome_spec = transcriptome_spec

    def query_genome(self, name, left_primer, right_primer, max_size) -> Sequence[GenomeHit]:
        """query for primer hits in genome"""
        genome_pcr_psls = _gfPcr(self.genome_spec, name, left_primer, right_primer, max_size)
        return [_genome_psl_to_hit(psl) for psl in genome_pcr_psls]

    def query_transcriptome(self, name, left_primer, right_primer, max_size) -> Sequence[TranscriptomeHit]:
        """query for primer hits in transcriptome"""
        trans_pcr_psls = _gfPcr(self.transcriptome_spec, name, left_primer, right_primer, max_size)
        trans_to_features = _trans_to_features(self.genome_data, self.transcriptome_spec, trans_pcr_psls)
        return [_trans_psl_to_hit(trans_to_features, psl) for psl in trans_pcr_psls]
