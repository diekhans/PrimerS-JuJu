"""
tests cover
   primersjuju.transcript_features
"""
from pycbio.hgdata.coords import Coords
from primersjuju.transcript_features import ExonFeature, IntronFeature, bed_to_features, features_intersect_genome, transcript_range_to_features

def transcript_region_check(genome_data_hg38, wtc11_track, trans_id, region, expected_feats):
    trans_bed = wtc11_track.read_by_name(trans_id)
    feats = bed_to_features(genome_data_hg38, trans_bed)
    subfeats = features_intersect_genome(feats, region)
    assert len(subfeats) == len(expected_feats)
    assert subfeats == expected_feats

def test_get_transcript_region_exon(genome_data_hg38, wtc11_track):
    transcript_region_check(genome_data_hg38, wtc11_track, "FSM_45093",
                            Coords("chr17", 7709209, 7710259, strand='+', size=83257441),
                            [ExonFeature(genome=Coords(name='chr17', start=7709209, end=7710259, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=1051, end=2101, strand='+', size=3209))
                             ])

def test_get_transcript_region_intron_exon(genome_data_hg38, wtc11_track):
    transcript_region_check(genome_data_hg38, wtc11_track, "FSM_45093",
                            Coords("chr17", 7708917, 7709357, strand='+', size=83257441),
                            [IntronFeature(genome=Coords(name='chr17', start=7708917, end=7709166, strand='+', size=83257441),
                                           trans=Coords(name='FSM_45093', start=1008, end=1008, strand='+', size=3209)),
                             ExonFeature(genome=Coords(name='chr17', start=7709166, end=7709357, strand='+', size=83257441),
                                         trans=Coords(name='FSM_45093', start=1008, end=1199, strand='+', size=3209))
                             ])

def test_get_transcript_region_exons(genome_data_hg38, wtc11_track):
    transcript_region_check(genome_data_hg38, wtc11_track, "FSM_45093",
                            Coords("chr17", 7708471, 7709256, strand='+', size=83257441),
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

def test_get_transcript_region_FSM_45580(genome_data_hg38, wtc11_track):
    transcript_region_check(genome_data_hg38, wtc11_track, "FSM_45580",
                            Coords('chr19', 47228327, 47231039, strand='+', size=83257441),
                            [ExonFeature(genome=Coords(name='chr19', start=47228327, end=47228446, strand='+', size=58617616),
                                         trans=Coords(name='FSM_45580', start=1456, end=1575, strand='-', size=1841)),
                             IntronFeature(genome=Coords(name='chr19', start=47228446, end=47230928, strand='+', size=58617616),
                                           trans=Coords(name='FSM_45580', start=1575, end=1575, strand='-', size=1841)),
                             ExonFeature(genome=Coords(name='chr19', start=47230928, end=47231039, strand='+', size=58617616),
                                         trans=Coords(name='FSM_45580', start=1575, end=1686, strand='-', size=1841))])

def test_genome_inter_pos_pos():
    feat = ExonFeature(Coords("chr24", 6, 12, strand='+', size=24),
                       Coords("tr12", 0, 12, strand='+', size=12))
    grange = Coords("chr24", 3, 8, strand='+', size=24)

    feat_intr1 = feat.intersect_genome(grange)
    expect = ExonFeature(genome=Coords(name='chr24', start=6, end=8, strand='+', size=24),
                         trans=Coords(name='tr12', start=0, end=2, strand='+', size=12))
    assert feat_intr1 == expect
    feat_intr2 = feat.intersect_genome(grange.reverse())
    assert feat_intr2 == expect

def test_trans_inter_pos_neg():
    feat = ExonFeature(Coords("chr24", 6, 12, strand='+', size=24),
                       Coords("tr12", 0, 12, strand='-', size=12))
    trange = Coords("tr12", 3, 8, strand='+', size=12)

    feat_intr1 = feat.intersect_transcript(trange)
    expect = ExonFeature(genome=Coords(name='chr24', start=10, end=15, strand='+', size=24),
                         trans=Coords(name='tr12', start=4, end=9, strand='-', size=12))
    assert feat_intr1 == expect
    feat_intr2 = feat.intersect_transcript(trange.reverse())
    assert feat_intr2 == expect

def test_transcript_mapping():
    # ZBTB45+1 / FSM_45682
    transcript_features = [
        ExonFeature(genome=Coords(name='chr19', start=58513529, end=58514310, strand='+', size=58617616),
                    trans=Coords(name='FSM_45682', start=0, end=781, strand='-', size=2136)),
        IntronFeature(genome=Coords(name='chr19', start=58514310, end=58516394, strand='+', size=58617616),
                      trans=Coords(name='FSM_45682', start=781, end=781, strand='-', size=2136)),
        ExonFeature(genome=Coords(name='chr19', start=58516394, end=58517673, strand='+', size=58617616),
                    trans=Coords(name='FSM_45682', start=781, end=2060, strand='-', size=2136)),
        IntronFeature(genome=Coords(name='chr19', start=58517673, end=58519741, strand='+', size=58617616),
                      trans=Coords(name='FSM_45682', start=2060, end=2060, strand='-', size=2136)),
        ExonFeature(genome=Coords(name='chr19', start=58519741, end=58519817, strand='+', size=58617616),
                    trans=Coords(name='FSM_45682', start=2060, end=2136, strand='-', size=2136))]

    primer_trans_coord_5p = Coords(name='FSM_45682', start=61, end=81, strand='+', size=2136)
    primer_trans_coord_3p = Coords(name='FSM_45682', start=1403, end=1423, strand='+', size=2136)

    genome_features_5p = transcript_range_to_features(transcript_features, primer_trans_coord_5p)

    assert genome_features_5p == [
        ExonFeature(genome=Coords(name='chr19', start=58517668, end=58517673, strand='+', size=58617616), trans=Coords(name='FSM_45682', start=2055, end=2060, strand='-', size=2136)),
        ExonFeature(genome=Coords(name='chr19', start=58519741, end=58519756, strand='+', size=58617616), trans=Coords(name='FSM_45682', start=2060, end=2075, strand='-', size=2136))
    ]
    genome_features_3p = transcript_range_to_features(transcript_features, primer_trans_coord_3p)
    assert genome_features_3p == [
        ExonFeature(genome=Coords(name='chr19', start=58514242, end=58514262, strand='+', size=58617616), trans=Coords(name='FSM_45682', start=713, end=733, strand='-', size=2136))
    ]
