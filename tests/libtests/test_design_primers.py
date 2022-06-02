"""
tests cover
   primersjuju.design_primers
"""

from primersjuju.design_primers import design_primers
from primersjuju.target_transcripts import target_transcripts_build

def _get_target_transcripts(genome_data, targets_specs, target_id):
    target_spec = targets_specs.get_target(target_id)
    return target_transcripts_build(genome_data, target_spec)

def test_exon_contained(genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _get_target_transcripts(genome_data, example_wtc11_targets_specs_set1, "SNAI1+1")

    primers = design_primers(genome_data, target_transcripts)
    assert primers.target_id == 'SNAI1+1'

    p3r = primers.designs[0].primer3_result
    assert p3r.PRIMER_LEFT == (27, 20)
    assert p3r.PRIMER_RIGHT == (1159, 20)
    assert p3r.PRIMER_LEFT_SEQUENCE == 'GGTTCTTCTGCGCTACTGCT'
    assert p3r.PRIMER_RIGHT_SEQUENCE == 'CAAAAACCCACGCAGACAGG'

    import pprint
    print('@@')
    pprint.pprint(vars(p3r))
