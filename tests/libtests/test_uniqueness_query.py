"""
tests cover
   primersjuju.uniqueness_query
"""
from pycbio.hgdata.coords import Coords

from primersjuju.transcript_features import ExonFeature
from primersjuju.uniqueness_query import GenomeHit, TranscriptomeHit

SNAI1_PRIMER_LEFT_SEQUENCE = 'GGTTCTTCTGCGCTACTGCT'
SNAI1_PRIMER_RIGHT_SEQUENCE = 'CAAAAACCCACGCAGACAGG'

ZBTB45_PRIMER_LEFT_SEQUENCE = 'AAGAGGAGTCGGAGATGGCG'
ZBTB45_PRIMER_RIGHT_SEQUENCE = 'GCGTAGAGAGAAGGATCGCC'

def test_uniqueness_query_genome(hg38_uniqueness_query):
    hits = hg38_uniqueness_query.query_genome("SNAI1+1+pp1", SNAI1_PRIMER_LEFT_SEQUENCE, SNAI1_PRIMER_RIGHT_SEQUENCE, 200000)
    assert hits == [GenomeHit(left_coords=Coords(name='chr20', start=49983006, end=49983026, strand='+', size=64444167),
                              right_coords=Coords(name='chr20', start=49988321, end=49988341, strand='+', size=64444167))]

def test_uniqueness_query_transcriptome(hg38_uniqueness_query):
    hits = hg38_uniqueness_query.query_transcriptome("SNAI1+1+pp1", SNAI1_PRIMER_LEFT_SEQUENCE, SNAI1_PRIMER_RIGHT_SEQUENCE, 200000)
    assert hits == [TranscriptomeHit(trans_id='ENST00000244050.3', gene_name='SNAI1',
                                     left_features=[ExonFeature(genome=Coords(name='chr20', start=49983006, end=49983026, strand='+', size=64444167),
                                                                trans=Coords(name='ENST00000244050.3', start=27, end=47, strand='+', size=1705))],
                                     right_features=[ExonFeature(genome=Coords(name='chr20', start=49988321, end=49988341, strand='+', size=64444167),
                                                                 trans=Coords(name='ENST00000244050.3', start=1140, end=1160, strand='+', size=1705))])]

def test_uniqueness_query_transcriptome_junc(hg38_uniqueness_query):
    hits = hg38_uniqueness_query.query_transcriptome("ZBTB45+1", ZBTB45_PRIMER_LEFT_SEQUENCE, ZBTB45_PRIMER_RIGHT_SEQUENCE, 200000)
    assert hits == [TranscriptomeHit(trans_id='ENST00000594051.6', gene_name='ZBTB45',
                                     left_features=[ExonFeature(genome=Coords(name='chr19', start=58514235, end=58514255, strand='+', size=58617616),
                                                                trans=Coords(name='ENST00000594051.6', start=55, end=75, strand='-', size=2129))],
                                     right_features=[ExonFeature(genome=Coords(name='chr19', start=58517056, end=58517076, strand='+', size=58617616),
                                                                 trans=Coords(name='ENST00000594051.6', start=1378, end=1398, strand='-', size=2129))])]
    # this failed to find one spanning the intron
