"""
tests cover
   primersjuju.design_primers
these are faster tests that check various function that don't use isPcr queries.
"""
import re
import pytest
from primersjuju import PrimersJuJuDataError
from .testfuncs import test_id, run_primer_design_test

def test_pos_strand(request, genome_data_hg38, wtc11_targets_specs_set1):
    # simple positive strand case, both regions in exons
    primer_designs = run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "SNAI1+1")

    assert primer_designs.target_id == 'SNAI1+1'
    assert len(primer_designs.designs) == 5

    p3p = primer_designs.designs[0].primer3_pair
    assert p3p.PRIMER_LEFT == (22, 20)
    assert p3p.PRIMER_RIGHT == (1026, 20)
    assert p3p.PRIMER_LEFT_SEQUENCE == 'CGAGTGGTTCTTCTGCGCTA'
    assert p3p.PRIMER_RIGHT_SEQUENCE == 'TCATCAAAGTCCTGTGGGGC'

def test_neg_strand(request, genome_data_hg38, wtc11_targets_specs_set1):
    # simple negative strand case, both regions in exons
    primer_designs = run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "ENO1+1")

    assert primer_designs.target_id == 'ENO1+1'
    assert len(primer_designs.designs) == 5

    p3p = primer_designs.designs[0].primer3_pair
    assert p3p.PRIMER_LEFT == (244, 20)
    assert p3p.PRIMER_RIGHT == (1348, 20)
    assert p3p.PRIMER_LEFT_SEQUENCE == 'TATGAGGCCCTAGAGCTCCG'
    assert p3p.PRIMER_RIGHT_SEQUENCE == 'CAAGAGCACTGACTCAGGGG'

def test_neg_strand_sj(request, genome_data_hg38, wtc11_targets_specs_set1):
    # regression test for incorrect specifying primer regions.
    # negative strand gene with primer results were mapped to
    # the wrong place in the gene.
    primer_designs = run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "ZBTB45+1")
    assert len(primer_designs.designs) == 5

def test_FBXL16(request, genome_data_hg38, wtc11_targets_specs_set1):
    # Negative strand, create primer that went beyond 3' end of region.  This maybe something
    # primer3 does

    primer_designs = run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "FBXL16+1")
    assert len(primer_designs.designs) == 5

def test_SLCA1(request, genome_data_hg38, wtc11_targets_specs_set1):
    # both regions in exons, multiple transcript in target
    primer_designs = run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                                            "SLC46A1+1")

    assert primer_designs.target_id == 'SLC46A1+1'
    assert len(primer_designs.designs) == 5

    p3p = primer_designs.designs[0].primer3_pair
    assert p3p.PRIMER_LEFT == (815, 20)
    assert p3p.PRIMER_RIGHT == (1662, 20)
    assert p3p.PRIMER_LEFT_SEQUENCE == 'CTCTTCACGTTCCGTCACCA'
    assert p3p.PRIMER_RIGHT_SEQUENCE == 'TCCTAGACAGAGGCTGGGTC'

def test_too_short(request, genome_data_hg38, wtc11_targets_specs_set1):
    # primer3 would generate PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE,
    # but we catch this upfront
    with pytest.raises(PrimersJuJuDataError,
                       match=re.escape('transcript WTC11_consolidated/NNC_64139 region chr19:50658490-50658506 exon length 16 is less than PRIMER_MIN_SIZE 18')):
        run_primer_design_test(test_id(request), genome_data_hg38, wtc11_targets_specs_set1,
                               "C19orf81+1")
