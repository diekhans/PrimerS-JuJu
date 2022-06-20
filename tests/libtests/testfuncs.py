"""
functions to support tests
"""
import os.path as osp
import pipettor
from primersjuju.primer_targets import primer_targets_build
from primersjuju.design_primers import design_primers
from primersjuju.output import output_target_designs

def test_id(request):
    return request.node.name

def diff_expected(rel_name):
    pipettor.run(["diff", "-u",
                  osp.join("expected", rel_name),
                  osp.join("output", rel_name)])

hub_urls = ["http://test.org/JuJu/hub.txt"]

def run_primer_design_test(test_id, genome_data, targets_specs, target_id, uniqueness_query=None):
    outdir = osp.join("output", test_id)
    target_spec = targets_specs.get_target(target_id)
    primer_targets = primer_targets_build(genome_data, target_spec)
    primer_designs = design_primers(genome_data, primer_targets, uniqueness_query=uniqueness_query)

    output_target_designs(outdir, primer_targets, primer_designs, hub_urls=hub_urls)

    output_suffixes = (".debug.txt", ".designs.tsv", ".primers.bed", ".target.bed")
    if uniqueness_query is not None:
        output_suffixes += (".genome-uniqueness.bed", ".transcriptome-uniqueness.bed")
    for suffix in output_suffixes:
        diff_expected(osp.join(test_id, primer_targets.target_id + suffix))
    return primer_designs
