"""
tests cover
   primersjuju.design_primers
"""
from .testfuncs import get_test_id, run_primer_design_test

def test_SNAI1(request, genome_data_hg38, wtc11_targets_specs_set1, hg38_uniqueness_query):
    # both regions in exons
    primer_designs = run_primer_design_test(get_test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "SNAI1+1", hg38_uniqueness_query)

    assert primer_designs.target_id == 'SNAI1+1'
    assert len(primer_designs.designs) == 5

    p3p = primer_designs.designs[0].primer3_pair
    assert p3p.PRIMER_LEFT == [22, 20]
    assert p3p.PRIMER_RIGHT == [1026, 20]
    assert p3p.PRIMER_LEFT_SEQUENCE == 'CGAGTGGTTCTTCTGCGCTA'
    assert p3p.PRIMER_RIGHT_SEQUENCE == 'TCATCAAAGTCCTGTGGGGC'


def test_BBC3(request, genome_data_hg38, wtc11_targets_specs_set1, hg38_uniqueness_query):
    # intron in one region, no primer3 results
    primer_designs = run_primer_design_test(get_test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "BBC3+1", hg38_uniqueness_query)
    assert primer_designs.target_id == 'BBC3+1'
    assert len(primer_designs.designs) == 0


def test_ZBTB45(request, genome_data_hg38, wtc11_targets_specs_set1, hg38_uniqueness_query):
    # intron in one region
    primer_designs = run_primer_design_test(get_test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "ZBTB45+1", hg38_uniqueness_query)
    assert primer_designs.target_id == 'ZBTB45+1'
    assert len(primer_designs.designs) == 5

def test_CERNA1(request, genome_data_hg38, wtc11_targets_specs_set1, hg38_uniqueness_query):
    # regression for incorrect uniqueness checks when results overlap target
    primer_designs = run_primer_design_test(get_test_id(request), genome_data_hg38, wtc11_targets_specs_set1, "CERNA1+1",
                                            hg38_uniqueness_query)
    assert primer_designs.target_id == 'CERNA1+1'
    assert len(primer_designs.designs) == 5
