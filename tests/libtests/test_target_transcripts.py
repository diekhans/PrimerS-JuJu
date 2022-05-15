"""
tests cover
   primersjuju.target_transcripts
"""

import pytest
from pycbio.hgdata.coords import Coords
from primersjuju.target_transcripts import (
    ExonRegion, IntronRegion, get_transcript_region_features,
    target_transcripts_build
)

@pytest.fixture(scope="session")
def wtc11_fsm(gdata):
    wtc11 = gdata.get_track("WTC11_consolidated")
    wrecs = wtc11.read_by_names(["FSM-45093", "NNC-318304"])
    return wrecs["FSM-45093"]

def test_get_transcript_region_exon(gdata, wtc11_fsm):
    feats = get_transcript_region_features(gdata, wtc11_fsm, Coords.parse("chr17:7709209-7710259"))
    assert len(feats) == 1
    assert feats[0] == ExonRegion(name='chr17', start=7709209, end=7710259, strand='+', size=83257441)

def test_get_transcript_region_intron_exon(gdata, wtc11_fsm):
    feats = get_transcript_region_features(gdata, wtc11_fsm, Coords.parse("chr17:7,708,917-7,709,357"))
    assert len(feats) == 2
    assert feats == [IntronRegion(name='chr17', start=7708917, end=7709166, strand='+', size=83257441),
                     ExonRegion(name='chr17', start=7709166, end=7709357, strand='+', size=83257441)]

def test_get_transcript_region_exons(gdata, wtc11_fsm):
    feats = get_transcript_region_features(gdata, wtc11_fsm, Coords.parse(" chr17:7,708,471-7,709,256"))
    assert len(feats) == 5
    assert feats == [ExonRegion(name='chr17', start=7708471, end=7708527, strand='+', size=83257441),
                     IntronRegion(name='chr17', start=7708527, end=7708634, strand='+', size=83257441),
                     ExonRegion(name='chr17', start=7708634, end=7708739, strand='+', size=83257441),
                     IntronRegion(name='chr17', start=7708739, end=7709166, strand='+', size=83257441),
                     ExonRegion(name='chr17', start=7709166, end=7709256, strand='+', size=83257441)]

def test_target_build(gdata, example_wtc11_targets_specs_set1):
    target_spec = example_wtc11_targets_specs_set1.get_target("C4orf48+1")
    target_transcripts = target_transcripts_build(gdata, target_spec)
