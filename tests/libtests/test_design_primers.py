"""
tests cover
   primersjuju.design_primers
"""
import os.path as osp
from primersjuju.primer_targets import primer_targets_build
from primersjuju.design_primers import design_primers
from primersjuju.output import output_target_designs
from .testfuncs import test_id, diff_expected

hub_urls = ["http://test.org/JuJu/hub.txt"]

def _get_primer_targets(genome_data, targets_specs, target_id):
    target_spec = targets_specs.get_target(target_id)
    return primer_targets_build(genome_data, target_spec)

def _run_primer_design_test(test_id, genome_data, primer_targets, uniqueness_query):
    outdir = osp.join("output", test_id)
    primer_designs = design_primers(genome_data, primer_targets, uniqueness_query=uniqueness_query)

    output_target_designs(outdir, primer_targets, primer_designs, hub_urls=hub_urls)

    for suffix in (".debug.txt", ".designs.tsv", ".primers.bed", ".target.bed",
                   ".genome-uniqueness.bed", ".transcriptome-uniqueness.bed"):
        diff_expected(osp.join(test_id, primer_targets.target_id + suffix))
    return primer_designs

def test_SNAI1(request, genome_data, example_wtc11_targets_specs_set1, hg38_uniqueness_query):
    # both regions in exons
    primer_targets = _get_primer_targets(genome_data, example_wtc11_targets_specs_set1, "SNAI1+1")
    primer_designs = _run_primer_design_test(test_id(request), genome_data, primer_targets, hg38_uniqueness_query)

    assert primer_designs.target_id == 'SNAI1+1'
    assert len(primer_designs.designs) == 5

    p3p = primer_designs.designs[0].primer3_pair
    assert p3p.PRIMER_LEFT == (27, 20)
    assert p3p.PRIMER_RIGHT == (1159, 20)
    assert p3p.PRIMER_LEFT_SEQUENCE == 'GGTTCTTCTGCGCTACTGCT'
    assert p3p.PRIMER_RIGHT_SEQUENCE == 'CAAAAACCCACGCAGACAGG'


def test_BBC3(request, genome_data, example_wtc11_targets_specs_set1, hg38_uniqueness_query):
    # intron in one region, no primer3 results
    primer_targets = _get_primer_targets(genome_data, example_wtc11_targets_specs_set1, "BBC3+1")
    primer_designs = _run_primer_design_test(test_id(request), genome_data, primer_targets, hg38_uniqueness_query)
    assert primer_designs.target_id == 'BBC3+1'
    assert len(primer_designs.designs) == 0


def test_ZBTB45(request, genome_data, example_wtc11_targets_specs_set1, hg38_uniqueness_query):
    # intron in one region
    primer_targets = _get_primer_targets(genome_data, example_wtc11_targets_specs_set1, "ZBTB45+1")
    primer_designs = _run_primer_design_test(test_id(request), genome_data, primer_targets, hg38_uniqueness_query)
    assert primer_designs.target_id == 'ZBTB45+1'
    assert len(primer_designs.designs) == 5

def test_CERNA1(request, genome_data, example_wtc11_targets_specs_set1, hg38_uniqueness_query):
    # regression for incorrect uniqueness checks when results overlap target

    primer_targets = _get_primer_targets(genome_data, example_wtc11_targets_specs_set1, "CERNA1+1")
    primer_designs = _run_primer_design_test(test_id(request), genome_data, primer_targets, hg38_uniqueness_query)
    assert primer_designs.target_id == 'CERNA1+1'
    assert len(primer_designs.designs) == 5
