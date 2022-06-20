"""
tests cover
   primersjuju.design_primers
these are faster tests that check various function that don't use isPcr queries.
"""
from .testfuncs import test_id, run_primer_design_test

def test_neg_strand(request, genome_data, wtc11_targets_specs_set1):
    # regression test for incorrect specifying primer regions.
    # negative strand gene with primer results were mapped to
    # the wrong place in the gene.
    primer_designs = run_primer_design_test(test_id(request), genome_data, wtc11_targets_specs_set1,
                                            "ZBTB45+1")
    assert len(primer_designs.designs) == 5

def test_FBXL16(request, genome_data, wtc11_targets_specs_set1):
    # Negative strand, create primer that went beyond 3' end of region.  This maybe something
    # primer3 does

    primer_designs = run_primer_design_test(test_id(request), genome_data, wtc11_targets_specs_set1,
                                            "FBXL16+1")
    assert len(primer_designs.designs) == 5
