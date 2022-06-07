"""
tests cover
   primersjuju.target_transcripts
   primersjuju.transcript_features
"""

import pytest
import re
from pycbio.hgdata.coords import Coords
from pycbio.sys import fileOps
from primersjuju.transcript_features import ExonFeature, IntronFeature, bed_to_features, features_intersect_genome
from primersjuju.target_transcripts import target_transcripts_build


@pytest.fixture(scope="session")
def wtc11(genome_data):
    return genome_data.get_track("WTC11_consolidated")

def _dump_trans_region(target_id, trans_track, trans_id, target_transcript):
    """use to write test info to files for testing in primer3"""
    with open(trans_id + ".txt", 'w') as fh:
        tt = target_transcript
        ok_region = "{},{},{},{}".format(tt.region_5p.trans.start + 1,
                                         len(tt.region_5p.trans),
                                         tt.region_3p.trans.start + 1,
                                         len(tt.region_3p.trans))
        fileOps.prRowv(fh, "target_id", target_id)
        fileOps.prRowv(fh, "trans_id", trans_track, trans_id)
        fileOps.prRowv(fh, "region_5p", tt.region_5p.genome, len(tt.region_5p.genome), "features:", len(tt.features_5p.features))
        fileOps.prRowv(fh, "region_3p", tt.region_3p.genome, len(tt.region_3p.genome), "features:", len(tt.features_3p.features))
        fileOps.prRowv(fh, "ampl_region_5p", tt.region_5p.trans.format(oneBased=True))
        fileOps.prRowv(fh, "ampl_region_3p", tt.region_3p.trans.format(oneBased=True))
        fileOps.prRowv(fh, "SEQUENCE_PRIMER_PAIR_OK_REGION", ok_region)
        fileOps.prRowv(fh, "amplicon", len(target_transcript.amplicon_trans_range), target_transcript.amplicon_trans_range)
        fileOps.prRowv(fh, "rna", len(tt.rna), tt.rna)

def transcript_region_check(genome_data, wtc11, trans_id, region, expected_feats):
    trans_bed = wtc11.read_by_name(trans_id)
    feats = bed_to_features(genome_data, trans_bed)
    subfeats = features_intersect_genome(feats, region)
    assert len(subfeats) == len(expected_feats)
    assert subfeats == expected_feats

def test_get_transcript_region_exon(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM_45093",
                            Coords("chr17", 7709209, 7710259),
                            [ExonFeature(genome=Coords(name='chr17', start=7709209, end=7710259, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=1051, end=2101, strand='+', size=3209))
                             ])

def test_get_transcript_region_intron_exon(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM_45093",
                            Coords("chr17", 7708917, 7709357),
                            [IntronFeature(genome=Coords(name='chr17', start=7708917, end=7709166, strand='+', size=83257441),
                                           trans=Coords(name='FSM_45093', start=1008, end=1008, strand='+', size=3209)),
                             ExonFeature(genome=Coords(name='chr17', start=7709166, end=7709357, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=1008, end=1199, strand='+', size=3209))
                             ])

def test_get_transcript_region_exons(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM_45093",
                            Coords("chr17", 7708471, 7709256),
                            [ExonFeature(genome=Coords(name='chr17', start=7708471, end=7708527, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=847, end=903, strand='+', size=3209)),
                             IntronFeature(genome=Coords(name='chr17', start=7708527, end=7708634, strand='+', size=83257441),
                                           trans=Coords(name='FSM_45093', start=903, end=903, strand='+', size=3209)),
                             ExonFeature(genome=Coords(name='chr17', start=7708634, end=7708739, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=903, end=1008, strand='+', size=3209)),
                             IntronFeature(genome=Coords(name='chr17', start=7708739, end=7709166, strand='+', size=83257441),
                                           trans=Coords(name='FSM_45093', start=1008, end=1008, strand='+', size=3209)),
                             ExonFeature(genome=Coords(name='chr17', start=7709166, end=7709256, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=1008, end=1098, strand='+', size=3209))
                             ])

def test_get_transcript_region_FSM_45580(genome_data, wtc11):
    transcript_region_check(genome_data, wtc11, "FSM_45580",
                            Coords('chr19', 47228327, 47231039),
                            [ExonFeature(genome=Coords(name='chr19', start=47228327, end=47228446, strand='-', size=58617616),
                                         trans=Coords(name='FSM_45580', start=266, end=385, strand='+', size=1841)),
                             IntronFeature(genome=Coords(name='chr19', start=47228446, end=47230928, strand='-', size=58617616),
                                           trans=Coords(name='FSM_45580', start=266, end=266, strand='+', size=1841)),
                             ExonFeature(genome=Coords(name='chr19', start=47230928, end=47231039, strand='-', size=58617616),
                                         trans=Coords(name='FSM_45580', start=155, end=266, strand='+', size=1841))
                             ])

def _target_build(genome_data, targets_specs, target_id, num_trans, region_5p, region_3p):
    target_spec = targets_specs.get_target(target_id)
    target_transcripts = target_transcripts_build(genome_data, target_spec)
    assert target_transcripts.target_id == target_id
    assert len(target_transcripts.transcripts) == num_trans
    return target_transcripts

def _check_target_transcripts(target_transcripts, region_5p, region_3p, seq_5p_re, seq_3p_re):
    assert target_transcripts.region_5p == region_5p
    assert target_transcripts.region_3p == region_3p
    assert len(target_transcripts.sequence_5p) == len(region_5p)
    assert len(target_transcripts.sequence_3p) == len(region_3p)
    assert re.match(seq_5p_re, target_transcripts.sequence_5p)
    assert re.match(seq_3p_re, target_transcripts.sequence_3p)

def _check_target_transcript(target_transcripts, trans_track, trans_id, rna_len,
                             features_5p, features_3p):
    target_transcript = target_transcripts.get_transcript(trans_track, trans_id)
    assert target_transcript.features_5p.features == features_5p
    assert target_transcript.features_3p.features == features_3p
    assert len(target_transcript.rna) == rna_len

def test_target_build_exonic(genome_data, example_wtc11_targets_specs_set1):
    region_5p = Coords("chr20", 49983001, 49983051, strand='+', size=64444167)
    region_3p = Coords("chr20", 49987886, 49988699, strand='+', size=64444167)
    target_transcripts = _target_build(genome_data, example_wtc11_targets_specs_set1, "SNAI1+1", 1,
                                       region_5p, region_3p)
    # coords were swapped
    _check_target_transcripts(target_transcripts, region_3p, region_5p,
                              "^CCTGTCCCCACTGCAGCCGTGCCTTCGCTG.+GGCTGTCACTTGTCGGGGGCCCAAGTGGGGTGCTCTG$",
                              "^CGAGTGGTTCTTCTGCGCTACTGCTGCGCGAATCGGCGACCCCAGTGCCT$")

    _check_target_transcript(target_transcripts, "WTC11_consolidated", "FSM_23673", 1703,
                             [ExonFeature(genome=Coords(name='chr20', start=49987886, end=49988699, strand='+', size=64444167),
                                          trans=Coords(name='FSM_23673', start=705, end=1518, strand='+', size=1703))],
                             [ExonFeature(genome=Coords(name='chr20', start=49983001, end=49983051, strand='+', size=64444167),
                                          trans=Coords(name='FSM_23673', start=22, end=72, strand='+', size=1703))])

def test_target_build_junc(genome_data, example_wtc11_targets_specs_set1):
    region_5p = Coords("chr4", 2043862, 2043922, strand='+', size=190214555)
    region_3p = Coords("chr4", 2042041, 2042366, strand='+', size=190214555)
    target_transcripts = _target_build(genome_data, example_wtc11_targets_specs_set1, "C4orf48+1", 1,
                                       region_5p, region_3p)

    # coords were not swapped
    _check_target_transcripts(target_transcripts, region_5p, region_3p,
                              "^GGAGACCCTACTGCTGCAGGCAGAGCGCCGTGCCCTGTGTGCCTGCTGGCCAGCGGGGCA$",
                              "^TGCTGGTCCCGGGGTCCCTGAACCGCGGTAAGGGCGGTGGTGCGGGCGTCCGAATGGGCGTTTTCTAGATACGGGGCGCGGACTAGAGGCTCGCTGGGCC.+TCTGCTGGCCATGGCCCCCCCGCcc$")

    _check_target_transcript(target_transcripts, "WTC11_consolidated", "FSM_48428", 422,
                             [ExonFeature(genome=Coords(name='chr4', start=2043862, end=2043922, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=321, end=381, strand='+', size=422))],
                             [ExonFeature(genome=Coords(name='chr4', start=2042041, end=2042068, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=45, end=72, strand='+', size=422)),
                              IntronFeature(genome=Coords(name='chr4', start=2042068, end=2042326, strand='+', size=190214555),
                                            trans=Coords(name='FSM_48428', start=72, end=72, strand='+', size=422)),
                              ExonFeature(genome=Coords(name='chr4', start=2042326, end=2042366, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=72, end=112, strand='+', size=422))])

def test_target_build_junc2(genome_data, example_wtc11_targets_specs_set1):
    region_5p = Coords('chr19', 47221197, 47221914, strand='-', size=58617616)
    region_3p = Coords('chr19', 47228327, 47231039, strand='-', size=58617616)

    target_transcripts = _target_build(genome_data, example_wtc11_targets_specs_set1, "BBC3+1", 1,
                                       region_5p, region_3p)
    # coords were not swapped
    _check_target_transcripts(target_transcripts, region_5p, region_3p,
                              "^AAGAGGAGCAGCAGCGGCACCGCCCCTCACCCTGGAGGGTCCTGTACAATC.+GTTCCAGCTGCAGGGGTGACACTGGGAGGGGGGGGCTCTCCTCTCGGTGCT.+AGCCAGCCGGCGGGTGGTGGGCATGCCTGCCTCACCTTCATCAGGGGGT$",
                              "^gacagccacagcagcagccgccgcggagagcggcGCTCGGCGGGCGCGCCC.+GGGAGGCTCTCCAGGCCAGCCAGGACCCGG.+CGGCCCGCGCCCCTTCCCGCTCGGCCGCCTGGTGCCCTCGGCAGT$")

    _check_target_transcript(target_transcripts, "WTC11_consolidated", "FSM_45580", 1841,
                             [ExonFeature(genome=Coords(name='chr19', start=47221197, end=47221914, strand='-', size=58617616),
                                          trans=Coords(name='FSM_45580', start=750, end=1467, strand='+', size=1841))],
                             [ExonFeature(genome=Coords(name='chr19', start=47228327, end=47228446, strand='-', size=58617616),
                                          trans=Coords(name='FSM_45580', start=266, end=385, strand='+', size=1841)),
                              IntronFeature(genome=Coords(name='chr19', start=47228446, end=47230928, strand='-', size=58617616),
                                            trans=Coords(name='FSM_45580', start=266, end=266, strand='+', size=1841)),
                              ExonFeature(genome=Coords(name='chr19', start=47230928, end=47231039, strand='-', size=58617616),
                                          trans=Coords(name='FSM_45580', start=155, end=266, strand='+', size=1841))])
