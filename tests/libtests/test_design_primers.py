"""
tests cover
   primersjuju.design_primers
"""
import os.path as osp
from pycbio.sys.svgcolors import SvgColors
from pycbio.sys import fileOps
from primersjuju.design_primers import design_primers, build_primer_beds
from primersjuju.primer3_interface import primer3_annotate_rna
from primersjuju.target_transcripts import target_transcripts_build, build_target_bed

def _get_target_transcripts(genome_data, targets_specs, target_id):
    target_spec = targets_specs.get_target(target_id)
    return target_transcripts_build(genome_data, target_spec)

def _write_beds(beds, bed_file):
    with open(bed_file, "w") as fh:
        for bed in beds:
            bed.write(fh)

def _run_primer_design_test(request, genome_data, target_transcripts):
    fileOps.ensureDir("output")
    dump_fh = open(osp.join("output", request.node.name + ".debug"), 'w')
    target_transcripts.dump(dump_fh)

    primers = design_primers(genome_data, target_transcripts, dump_fh=dump_fh)
    primer_beds = build_primer_beds(primers, SvgColors.darkmagenta)
    target_bed = build_target_bed(target_transcripts, SvgColors.blue)

    print("annotated_rna:", primer3_annotate_rna(target_transcripts.transcripts[0]), file=dump_fh)
    _write_beds(primer_beds, osp.join("output", request.node.name + ".primers.bed"))
    _write_beds([target_bed], osp.join("output", request.node.name + ".target.bed"))
    return primers, primer_beds, target_bed

def test_exon_contained(request, genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _get_target_transcripts(genome_data, example_wtc11_targets_specs_set1, "SNAI1+1")

    primers, primer_beds, target_bed = _run_primer_design_test(request, genome_data, target_transcripts)

    assert primers.target_id == 'SNAI1+1'
    assert len(primers.designs) == 5

    p3r = primers.designs[0].primer3_result
    assert p3r.PRIMER_LEFT == (27, 20)
    assert p3r.PRIMER_RIGHT == (1159, 20)
    assert p3r.PRIMER_LEFT_SEQUENCE == 'GGTTCTTCTGCGCTACTGCT'
    assert p3r.PRIMER_RIGHT_SEQUENCE == 'CAAAAACCCACGCAGACAGG'

    # BED testing
    assert len(primer_beds) == 5
    assert str(primer_beds[0]) == ("chr20	49983005	49988359	SNAI1+1+pp1	0	+	49983005	49988359	139,0,139	2	20,20,	0,5334,	"
                                   "1133	12.50	5.03	0.42	(27, 20)	GGTTCTTCTGCGCTACTGCT	4.24	55.00	36.97	0.39	9.73	0.00	60.39	(1159, 20)	CAAAAACCCACGCAGACAGG	4.00	55.00	0.00	0.03	0.00	0.00	59.97")

def test_intron_containing(request, genome_data, example_wtc11_targets_specs_set1):
    # this doesn't find anything
    target_transcripts = _get_target_transcripts(genome_data, example_wtc11_targets_specs_set1, "BBC3+1")

    primers, primer_beds, target_bed = _run_primer_design_test(request, genome_data, target_transcripts)
    assert primers.target_id == 'BBC3+1'
    assert len(primers.designs) == 0
    assert len(primer_beds) == 0
