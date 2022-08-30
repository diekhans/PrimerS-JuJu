"""
tests cover
   primersjuju.primer_targets
"""
import pytest_check as check
from pycbio.hgdata.coords import Coords
from primersjuju.transcript_features import ExonFeature, IntronFeature
from primersjuju.primer_targets import primer_targets_build


def _target_build(genome_data_hg38, targets_specs, target_id, num_trans, region_5p, region_3p):
    target_spec = targets_specs.get_target(target_id)
    primer_targets = primer_targets_build(genome_data_hg38, target_spec)
    check.equal(primer_targets.target_id, target_id)
    check.equal(len(primer_targets.transcripts), num_trans)
    return primer_targets

def _check_primer_targets(primer_targets, region_5p, region_3p):
    check.equal(primer_targets.region_5p, region_5p)
    check.equal(primer_targets.region_3p, region_3p)

def _check_target_transcript(primer_targets, trans_track, trans_id, rna_len,
                             features_5p, features_3p):
    target_transcript = primer_targets.get_transcript(trans_track, trans_id)
    check.equal(target_transcript.features_5p, features_5p)
    check.equal(target_transcript.features_3p, features_3p)
    check.equal(len(target_transcript.rna), rna_len)

def test_target_build_exonic(genome_data_hg38, wtc11_targets_specs_set1):
    region_5p = Coords("chr20", 49983001, 49983051, strand='+', size=64444167)
    region_3p = Coords("chr20", 49987886, 49988699, strand='+', size=64444167)
    primer_targets = _target_build(genome_data_hg38, wtc11_targets_specs_set1, "SNAI1+1", 1,
                                   region_5p, region_3p)
    # coords were not swapped
    _check_primer_targets(primer_targets, region_5p, region_3p)

    _check_target_transcript(primer_targets, "WTC11_consolidated", "FSM_23673", 1703,
                             [ExonFeature(genome=Coords(name='chr20', start=49983001, end=49983051, strand='+', size=64444167),
                                          trans=Coords(name='FSM_23673', start=22, end=72, strand='+', size=1703))],
                             [ExonFeature(genome=Coords(name='chr20', start=49987886, end=49988699, strand='+', size=64444167),
                                          trans=Coords(name='FSM_23673', start=705, end=1518, strand='+', size=1703))])

def test_target_build_junc(genome_data_hg38, wtc11_targets_specs_set1):
    region_5p = Coords("chr4", 2043862, 2043922, strand='+', size=190214555)
    region_3p = Coords("chr4", 2042041, 2042366, strand='+', size=190214555)
    primer_targets = _target_build(genome_data_hg38, wtc11_targets_specs_set1, "C4orf48+1", 1,
                                   region_5p, region_3p)

    # coords were swapped
    _check_primer_targets(primer_targets, region_3p, region_5p)
    _check_target_transcript(primer_targets, "WTC11_consolidated", "FSM_48428", 422,
                             [ExonFeature(genome=Coords(name='chr4', start=2042041, end=2042068, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=45, end=72, strand='+', size=422)),
                              IntronFeature(genome=Coords(name='chr4', start=2042068, end=2042326, strand='+', size=190214555),
                                            trans=Coords(name='FSM_48428', start=72, end=72, strand='+', size=422)),
                              ExonFeature(genome=Coords(name='chr4', start=2042326, end=2042366, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=72, end=112, strand='+', size=422))],
                             [ExonFeature(genome=Coords(name='chr4', start=2043862, end=2043922, strand='+', size=190214555),
                                          trans=Coords(name='FSM_48428', start=321, end=381, strand='+', size=422))])

def test_target_build_junc2(genome_data_hg38, wtc11_targets_specs_set1):
    region_5p = Coords('chr19', 47221197, 47221914, strand='+', size=58617616)
    region_3p = Coords('chr19', 47228327, 47231039, strand='+', size=58617616)

    primer_targets = _target_build(genome_data_hg38, wtc11_targets_specs_set1, "BBC3+1", 1,
                                   region_5p, region_3p)
    # coords were swapped
    _check_primer_targets(primer_targets, region_3p, region_5p)
    _check_target_transcript(primer_targets, "WTC11_consolidated", "FSM_45580", 1841,
                             [ExonFeature(genome=Coords(name='chr19', start=47228327, end=47228446, strand='+', size=58617616),
                                          trans=Coords(name='FSM_45580', start=1456, end=1575, strand='-', size=1841)),
                              IntronFeature(genome=Coords(name='chr19', start=47228446, end=47230928, strand='+', size=58617616),
                                            trans=Coords(name='FSM_45580', start=1575, end=1575, strand='-', size=1841)),
                              ExonFeature(genome=Coords(name='chr19', start=47230928, end=47231039, strand='+', size=58617616),
                                          trans=Coords(name='FSM_45580', start=1575, end=1686, strand='-', size=1841))],
                             [ExonFeature(genome=Coords(name='chr19', start=47221197, end=47221914, strand='+', size=58617616),
                                          trans=Coords(name='FSM_45580', start=374, end=1091, strand='-', size=1841))])

def test_target_build_trim_sj(genome_data_hg38, wtc11_targets_specs_set1):
    # splice junction with slop into introns being adjusted

    region_5p = Coords.parse("chr2:218217185-218225938", strand='+', size=242193529, oneBased=True)
    region_3p = Coords.parse("chr2:218238453-218240799", strand='+', size=242193529, oneBased=True)

    primer_targets = _target_build(genome_data_hg38, wtc11_targets_specs_set1, "ARPC2+1", 1,
                                   region_5p, region_3p)
    # coords were swapped
    _check_primer_targets(primer_targets,
                          Coords(name='chr2', start=218217184, end=218225938, strand='+', size=242193529),
                          Coords(name='chr2', start=218238663, end=218239484, strand='+', size=242193529))
    _check_target_transcript(primer_targets, "WTC11_consolidated", "NIC_485306", 1339,
                             [ExonFeature(genome=Coords(name='chr2', start=218217184, end=218217254, strand='+', size=242193529),
                                          trans=Coords(name='NIC_485306', start=6, end=76, strand='+', size=1339)),
                              IntronFeature(genome=Coords(name='chr2', start=218217254, end=218225919, strand='+', size=242193529),
                                            trans=Coords(name='NIC_485306', start=76, end=76, strand='+', size=1339)),
                              ExonFeature(genome=Coords(name='chr2', start=218225919, end=218225938, strand='+', size=242193529),
                                          trans=Coords(name='NIC_485306', start=76, end=95, strand='+', size=1339))],
                             [ExonFeature(genome=Coords(name='chr2', start=218238663, end=218238850, strand='+', size=242193529),
                                          trans=Coords(name='NIC_485306', start=270, end=457, strand='+', size=1339)),
                              IntronFeature(genome=Coords(name='chr2', start=218238850, end=218239390, strand='+', size=242193529),
                                            trans=Coords(name='NIC_485306', start=457, end=457, strand='+', size=1339)),
                              ExonFeature(genome=Coords(name='chr2', start=218239390, end=218239484, strand='+', size=242193529),
                                          trans=Coords(name='NIC_485306', start=457, end=551, strand='+', size=1339))])

def test_target_build_trim_exon1(genome_data_hg38, wtc11_targets_specs_set1):
    # slop on one side of an exon

    region_5p = Coords.parse("chr19:35944582-35945182", strand='+', size=242193529, oneBased=True)
    region_3p = Coords.parse("chr19:35936839-35937541", strand='+', size=242193529, oneBased=True)

    primer_targets = _target_build(genome_data_hg38, wtc11_targets_specs_set1, "LRFN3+2", 1,
                                   region_5p, region_3p)
    # coords were swapped
    _check_primer_targets(primer_targets,
                          Coords(name='chr19', start=35936838, end=35937541, strand='+', size=58617616),
                          Coords(name='chr19', start=35944581, end=35945182, strand='+', size=58617616))
    _check_target_transcript(primer_targets, "WTC11_consolidated", "FSM_45884", 3633,
                             [ExonFeature(genome=Coords(name='chr19', start=35936838, end=35937541, strand='+', size=58617616),
                                          trans=Coords(name='FSM_45884', start=460, end=1163, strand='+', size=3633))],
                             [ExonFeature(genome=Coords(name='chr19', start=35944581, end=35945182, strand='+', size=58617616),
                                          trans=Coords(name='FSM_45884', start=2642, end=3243, strand='+', size=3633))])
