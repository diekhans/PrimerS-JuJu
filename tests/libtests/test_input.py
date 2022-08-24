"""
tests cover
   primersjuju.primer_targets
   primersjuju.genome_data
"""


import pytest
import io
from . import assert_except_msg
from pycbio.hgdata.coords import Coords, CoordsError
from pycbio.hgdata.bed import Bed
from primersjuju import PrimersJuJuDataError
from primersjuju.primer_target_spec import primer_target_specs_read

PRIMER_TARGETS_1 = (
    "target_id\ttrans_track\ttrans_id\t\tregion_5p\tregion_3p",
    "EFNB3+1\tGENCODE_V39\tENST00000226091.3\t\tchr17:7709210-7710259\tchr17:7705244-7705707",
    "EFNB3+1\tWTC11_consolidated\tFSM_45093\t\t\t",
    "EFNB3+1\tWTC11_consolidated\tNNC_318304\t\t\t",
)

def _mk_tsv(rows):
    "lines of tab-separated strings"
    lines = "\n".join(rows) + '\n'
    return io.StringIO(lines)

def test_primer_targets_example(example_targets_specs):
    target_ids = sorted(example_targets_specs.targets.keys())
    assert target_ids == ["EFNB3+1"]

    target = example_targets_specs.get_target("EFNB3+1")
    trans_ids = sorted(target.get_tracks_trans())
    assert trans_ids == [('GENCODE_V39', 'ENST00000226091.3'),
                         ('WTC11_consolidated', 'FSM_45093'),
                         ('WTC11_consolidated', 'NNC_318304')]

    assert isinstance(target.region_5p, Coords)

    trans = target.get_transcript('GENCODE_V39', 'ENST00000226091.3')
    assert str(trans) == "(GENCODE_V39, ENST00000226091.3) [('trans_cat', ''), ('trans_notes', ''), ('trans_structural_cat', '')]"

    trans = target.get_transcript('WTC11_consolidated', 'FSM_45093')
    assert str(trans) == "(WTC11_consolidated, FSM_45093) [('trans_cat', 'FSM'), ('trans_notes', ''), ('trans_structural_cat', 'PB-specific')]"


def test_primer_targets_bad_header():
    tsv = list(PRIMER_TARGETS_1)
    tsv[0] = tsv[0].replace("region_5p", "fred_5p")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuDataError) as exinfo:
        primer_target_specs_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, PrimersJuJuDataError,
                      msg="required column is missing: 'region_5p'")

def test_primer_targets_bad_coords():
    tsv = list(PRIMER_TARGETS_1)
    tsv[1] = tsv[1].replace("chr17:7709210-7710259", "chr17:77XZ210-7710259")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuDataError) as exinfo:
        primer_target_specs_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, CoordsError,
                      msg="invalid coordinates: 'chr17:77XZ210-7710259'")

def test_primer_targets_overlapping_coords():
    tsv = list(PRIMER_TARGETS_1)
    # overlap with chr17:7709210-7710259
    tsv[1] = tsv[1].replace("chr17:7705244-7705707", "chr17:7709220-7710299")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuDataError) as exinfo:
        primer_target_specs_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, PrimersJuJuDataError,
                      msg="region_5p (chr17:7709209-7710259) overlaps region_3p (chr17:7709219-7710299)")

def test_primer_targets_diff_chrom_coords():
    tsv = list(PRIMER_TARGETS_1)
    tsv[1] = tsv[1].replace("chr17:7709210-7710259", "chr22:7709210-7710259")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuDataError) as exinfo:
        primer_target_specs_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, PrimersJuJuDataError,
                      msg="region_5p (chr22:7709209-7710259) is on a different chromosome than region_3p (chr17:7705243-7705707)")

def test_primer_targets_dup_target_id():
    tsv = list(PRIMER_TARGETS_1)
    tsv.append(PRIMER_TARGETS_1[1])
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuDataError) as exinfo:
        primer_target_specs_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, PrimersJuJuDataError,
                      msg="duplicate primer target_id 'EFNB3+1'")

def test_get_genome_seq(genome_data_hg38):
    coords = Coords.parse("chr17:7709209-7709259")
    assert genome_data_hg38.get_genome_seq(coords) == 'CCTGCCCCCTCCCAGCATGCCTGCAGTGGCTGGGGCAGCAGGGGGGCTGG'

def test_get_track_annot(genome_data_hg38):
    gencode = genome_data_hg38.get_track("gencodeV39")
    grecs = gencode.read_by_names(["ENST00000226091.3"])
    assert len(grecs) == 1

    wtc11 = genome_data_hg38.get_track("WTC11_consolidated")
    wrecs = wtc11.read_by_names(["FSM_45093", "NNC_318304"])
    assert len(wrecs) == 2
    assert isinstance(wrecs['FSM_45093'], Bed)

def test_get_track_missing(genome_data_hg38):
    wtc11 = genome_data_hg38.get_track("WTC11_consolidated")

    with pytest.raises(PrimersJuJuDataError,
                       match='track WTC11_consolidated: records not found in bigBed.*: Barney, Fred'):
        wtc11.read_by_names(["FSM_45093", "Fred", "Barney"])
