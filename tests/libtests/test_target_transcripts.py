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
def wtc11(genome_data):
    return genome_data.get_track("WTC11_consolidated")

def transcript_region_check(genome_data, wtc11, trans_id, crange, expected_feats):
    fsm = wtc11.read_by_name(trans_id)
    feats = get_transcript_region_features(genome_data, fsm, crange)
    assert len(feats) == len(expected_feats)
    assert feats == expected_feats

def test_get_transcript_region_exon(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM-45093",
                            Coords("chr17", 7709209, 7710259),
                            [ExonRegion('chr17', 7709209, 7710259, strand='+', size=83257441)])

def test_get_transcript_region_intron_exon(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM-45093",
                            Coords("chr17", 7708917, 7709357),
                            [IntronRegion('chr17', 7708917, 7709166, strand='+', size=83257441),
                             ExonRegion('chr17', 7709166, 7709357, strand='+', size=83257441)])

def test_get_transcript_region_exons(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM-45093",
                            Coords("chr17", 7708471, 7709256),
                            [ExonRegion('chr17', 7708471, 7708527, strand='+', size=83257441),
                             IntronRegion('chr17', 7708527, 7708634, strand='+', size=83257441),
                             ExonRegion('chr17', 7708634, 7708739, strand='+', size=83257441),
                             IntronRegion('chr17', 7708739, 7709166, strand='+', size=83257441),
                             ExonRegion('chr17', 7709166, 7709256, strand='+', size=83257441)])

def test_get_transcript_region_FSM_45580(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM-45580",
                            Coords('chr19', 47228327, 47231039),
                            [ExonRegion('chr19', 47228327, 47228446, strand='-', size=58617616),
                             IntronRegion('chr19', 47228446, 47230928, strand='-', size=58617616),
                             ExonRegion('chr19', 47230928, 47231039, strand='-', size=58617616)])

def _target_build(genome_data, targets_specs, target_id, expect_num_trans, expect_region_5p, expect_region_3p):
    target_spec = targets_specs.get_target(target_id)
    target_transcripts = target_transcripts_build(genome_data, target_spec)
    assert target_transcripts.target_id == target_id
    assert len(target_transcripts.transcripts) == expect_num_trans
    assert target_transcripts.region_5p == expect_region_5p
    assert target_transcripts.region_3p == expect_region_3p
    return target_transcripts

def _check_target_transcript(target_transcripts, trans_track, trans_id, expect_features_5p, expect_features_3p):
    target_transcript = target_transcripts.get_transcript(trans_track, trans_id)
    assert target_transcript.features_5p.features == expect_features_5p
    assert target_transcript.features_3p.features == expect_features_3p

def test_target_build(genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _target_build(genome_data, example_wtc11_targets_specs_set1, "C4orf48+1", 1,
                                       Coords("chr4", 2042041, 2042366, strand='+', size=190214555),
                                       Coords("chr4", 2043862, 2043922, strand='+', size=190214555))
    _check_target_transcript(target_transcripts, "WTC11_consolidated", "FSM-48428",
                             [ExonRegion(name='chr4', start=2042041, end=2042068, strand='+', size=190214555),
                              IntronRegion(name='chr4', start=2042068, end=2042326, strand='+', size=190214555),
                              ExonRegion(name='chr4', start=2042326, end=2042366, strand='+', size=190214555)],
                             [ExonRegion(name='chr4', start=2043862, end=2043922, strand='+', size=190214555)])

def test_target_build_junc(genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _target_build(genome_data, example_wtc11_targets_specs_set1, "BBC3+1", 1,
                                       Coords('chr19', 47228327, 47231039, strand='-', size=58617616),
                                       Coords('chr19', 47221197, 47221914, strand='-', size=58617616))
    _check_target_transcript(target_transcripts, "WTC11_consolidated", "FSM-45580",
                             [ExonRegion('chr19', 47228327, 47228446, strand='-', size=58617616),
                              IntronRegion('chr19', 47228446, 47230928, strand='-', size=58617616),
                              ExonRegion('chr19', 47230928, 47231039, strand='-', size=58617616)],
                             [ExonRegion('chr19', 47221197, 47221914, strand='-', size=58617616)])
