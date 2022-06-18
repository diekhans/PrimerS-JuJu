"""
tests cover
   primersjuju.design_primers
these are faster tests that check various function that don't require isPcr queries.
"""
import os.path as osp
from pycbio.sys import fileOps
from primersjuju.primer_targets import primer_targets_build
from primersjuju.design_primers import design_primers, _genome_uniqueness_classify, _transcriptome_uniqueness_check
from primersjuju.output import output_target_beds
from .testfuncs import test_id, diff_expected

def _diff_target_beds(test_id, primer_targets):
    for suffix in (".primers.bed", ".target.bed",
                   ".genome-uniqueness.bed", ".transcriptome-uniqueness.bed"):
        diff_expected(osp.join(test_id, primer_targets.target_id + suffix))

def test_neg_strand(request, genome_data, example_wtc11_targets_specs_set1):
    # regression test for incorrect specifying primer regions.
    # negative strand gene with primer results were mapped to
    # the wrong place in the gene.

    target_spec = example_wtc11_targets_specs_set1.get_target("ZBTB45+1")
    primer_targets = primer_targets_build(genome_data, target_spec)
    # no uniqueness to make this faster
    primer_designs = design_primers(genome_data, primer_targets)
    output_target_beds(osp.join("output", test_id(request)), primer_targets, primer_designs)
    _diff_target_beds(test_id(request), primer_targets)
