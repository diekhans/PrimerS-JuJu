"""
tests cover
   primersjuju.uniqueness_query
"""
import pytest
import re
from pycbio.hgdata.coords import Coords
from pycbio.sys import fileOps
from primersjuju.transcript_features import ExonFeature, IntronFeature, bed_to_features, features_intersect_genome
from primersjuju.target_transcripts import target_transcripts_build


SNAI1_PRIMER_LEFT_SEQUENCE = 'GGTTCTTCTGCGCTACTGCT'
SNAI1_PRIMER_RIGHT_SEQUENCE = 'CAAAAACCCACGCAGACAGG'


def test_uniqueness_query_genome(hg38_uniqueness_query):
    ps = hg38_uniqueness_query.query_genome("SNAI1+1+pp1", SNAI1_PRIMER_LEFT_SEQUENCE, SNAI1_PRIMER_RIGHT_SEQUENCE, 200000)
    for p in ps:
        print(p)

def test_uniqueness_query_transcriptome(hg38_uniqueness_query):
    ps = hg38_uniqueness_query.query_transcriptome("SNAI1+1+pp1", SNAI1_PRIMER_LEFT_SEQUENCE, SNAI1_PRIMER_RIGHT_SEQUENCE, 200000)
    for p in ps:
        print(p)
