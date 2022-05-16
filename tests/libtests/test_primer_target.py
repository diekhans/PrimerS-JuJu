"""
tests cover
   primersjuju.primer_target
"""

import pytest
from pycbio.hgdata.coords import Coords
from primersjuju.primer_target import primer_target_build
from primersjuju.target_transcripts import target_transcripts_build


def _get_target_transcripts(genome_data, targets_specs, target_id):
    target_spec = targets_specs.get_target(target_id)
    return target_transcripts_build(genome_data, target_spec)


def test_exon_contained(genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _get_target_transcripts(genome_data, example_wtc11_targets_specs_set1, "C4orf48+1")
    primer_target = primer_target_build(genome_data, target_transcripts)

    assert len(primer_target.region_seq_5p) == 325
    assert primer_target.region_seq_5p == (
        "TGCTGGTCCCGGGGTCCCTGAACCGCGGTAAGGGCGGTGGTGCGGGCGTCCGAATGGGCGTTTTCTAGATACG"
        "GGGCGCGGACTAGAGGCTCGCTGGGCCCGGAGACCGGCGGACTGGAGTCGGGGAACCGGAGGTGGGGAGGGGG"
        "CTCCCGGGCCCGGGGTGGGTGGGTCCAGGGCTCCCAGGCCTGGGGCTTGGACAAGGGTCGTGGGGCCCgcggg"
        "aggggacgggggctcaccggcccggggcgggcggggcgggCGCCGCTGACCCCTCGCTGGCTTCAGGGCGGCCC"
        "CGCTCCCTCTGCTGGCCATGGCCCCCCCGCcc"
    )
    assert len(primer_target.region_seq_3p) == 60
    assert primer_target.region_seq_3p == (
        "GGAGACCCTACTGCTGCAGGCAGAGCGCCGTGCCCTGTGTGCCTGCTGGCCAGCGGGGCA"
    )
