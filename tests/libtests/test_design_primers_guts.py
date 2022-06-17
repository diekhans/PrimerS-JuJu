"""
tests cover
   primersjuju.design_primers
these are faster tests that check various function that don't require isPcr queries.
"""
import os.path as osp
import pickle
from pycbio.sys import fileOps
from primersjuju.primer_targets import primer_targets_build
from primersjuju.design_primers import design_primers, _genome_uniqueness_classify, _transcriptome_uniqueness_check
from primersjuju.output import output_target_beds
from .testfuncs import test_id, diff_expected

def _serialize_path(test_id):
    return osp.join("input", test_id + ".serial.pickle")

def _store_design(test_id, primer_targets, primer_designs):
    "store design for faster load"
    root = {"primer_targets": primer_targets,
            "primer_designs": primer_designs}
    with fileOps.AtomicFileCreate(_serialize_path(test_id)) as tmp_file:
        with open(tmp_file, 'wb') as fh:
            pickle.dump(root, fh)

def _load_design(test_id):
    with open(_serialize_path(test_id), 'rb') as fh:
        root = pickle.load(fh)
    return root["primer_targets"], root["primer_designs"]

def _update_hit_classification(genome_data, primer_targets, primer_designs):
    """Updates hits, which would be used in debugging classification problems using
    serialized data.  This is all done because it was so time consuming to do debug cycles
    with gfPcr """
    def _redo_classify(design):
        "merge hits and reclassify"
        target_transcript = primer_targets.transcripts[0]
        genome_hits = design.genome_on_targets + design.genome_off_targets + design.genome_non_targets
        transcriptome_hits = design.transcriptome_on_targets + design.transcriptome_off_targets + design.transcriptome_non_targets

        design.genome_on_targets, design.genome_off_targets, design.genome_non_targets = \
            _genome_uniqueness_classify(genome_data, target_transcript, genome_hits)
        design.transcriptome_on_targets, design.transcriptome_off_targets, design.transcriptome_non_targets = \
            _transcriptome_uniqueness_check(genome_data, target_transcript, transcriptome_hits)

    for design in primer_designs.designs:
        _redo_classify(design)


def test_CERNA1_uniqueness(request, genome_data, example_wtc11_targets_specs_set1, hg38_uniqueness_query):
    "test uniqueness checks"
    # this is a slow query, so query results are serialized and save in repo
    target_spec = example_wtc11_targets_specs_set1.get_target("CERNA1+1")

    if not osp.exists(_serialize_path(test_id(request))):
        primer_targets = primer_targets_build(genome_data, target_spec)
        primer_designs = design_primers(genome_data, primer_targets, uniqueness_query=hg38_uniqueness_query)
        _store_design(test_id(request), primer_targets, primer_designs)
    else:
        primer_targets, primer_designs = _load_design(test_id(request))

    _update_hit_classification(genome_data, primer_targets, primer_designs)
    output_target_beds(osp.join("output", test_id(request)), primer_targets, primer_designs)

    for suffix in (".primers.bed", ".target.bed",
                   ".genome-uniqueness.bed", ".transcriptome-uniqueness.bed"):
        diff_expected(osp.join(test_id(request), primer_targets.target_id + suffix))
