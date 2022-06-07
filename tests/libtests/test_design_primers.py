"""
tests cover
   primersjuju.design_primers
"""
from pycbio.sys.svgcolors import SvgColors
from primersjuju.design_primers import design_primers, build_primer_beds
from primersjuju.target_transcripts import target_transcripts_build, build_target_bed

def _get_target_transcripts(genome_data, targets_specs, target_id):
    target_spec = targets_specs.get_target(target_id)
    return target_transcripts_build(genome_data, target_spec)

def _get_bed_lines(primers):
    beds = build_primer_beds(primers, SvgColors.darkmagenta)
    return [str(b) for b in beds]

def test_exon_contained(genome_data, example_wtc11_targets_specs_set1):
    target_transcripts = _get_target_transcripts(genome_data, example_wtc11_targets_specs_set1, "SNAI1+1")

    primers = design_primers(genome_data, target_transcripts)
    assert primers.target_id == 'SNAI1+1'
    assert len(primers.designs) == 5

    p3r = primers.designs[0].primer3_result
    assert p3r.PRIMER_LEFT == (27, 20)
    assert p3r.PRIMER_RIGHT == (1159, 20)
    assert p3r.PRIMER_LEFT_SEQUENCE == 'GGTTCTTCTGCGCTACTGCT'
    assert p3r.PRIMER_RIGHT_SEQUENCE == 'CAAAAACCCACGCAGACAGG'

    # BED testing
    bed_rows = _get_bed_lines(primers)
    assert len(bed_rows) == 5
    assert bed_rows[0] == ("chr20	49983005	49988359	SNAI1+1+pp1	0	+	49983005	49988359	139,0,139	2	20,20,	0,5334,	"
                           "1133	12.50	5.03	0.42	(27, 20)	GGTTCTTCTGCGCTACTGCT	4.24	55.00	36.97	0.39	9.73	0.00	60.39	(1159, 20)	CAAAAACCCACGCAGACAGG	4.00	55.00	0.00	0.03	0.00	0.00	59.97")

    from pycbio.sys import fileOps  #TMP
    fileOps.writeLines("primers.bed", bed_rows)

    tbed = build_target_bed(target_transcripts, SvgColors.blue)
    with open("targets.bed", "w") as fh:
        tbed.write(fh)
