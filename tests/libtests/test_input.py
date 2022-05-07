import pytest
import os.path as osp
import io
from . import mydir, assert_except_msg
from pycbio.hgdata.coords import Coords, CoordsError
from primersjuju import PrimersJuJuUserError
from primersjuju.primer_targets import primer_targets_read

PRIMER_TARGETS_1 = (
    "target_id\ttrans_track\ttrans_id\t\tregion_5p\tregion_3p",
    "EFNB3+1\tGENCODE_V39\tENST00000226091.3\t\tchr17:7709210-7710259\tchr17:7705244-7705707",
    "EFNB3+1\tWTC11_consolidated\tFSM_45093\t\t\t",
    "EFNB3+1\tWTC11_consolidated\tNNC-318304\t\t\t",
)

def _mk_tsv(rows):
    "lines of tab-separated strings"
    lines = "\n".join(rows) + '\n'
    return io.StringIO(lines)

def test_primer_targets_example():
    targets = primer_targets_read(osp.join(mydir, "../../docs/primer-targets-example.tsv"))

    target_ids = sorted(targets.targets.keys())
    assert target_ids == ["EFNB3+1"]

    target = targets.access_target("EFNB3+1")
    trans_ids = sorted(target.get_tracks_trans())
    assert trans_ids == [('GENCODE_V39', 'ENST00000226091.3'),
                         ('WTC11_consolidated', 'FSM_45093'),
                         ('WTC11_consolidated', 'NNC-318304')]

    assert isinstance(target.region_5p, Coords)

    trans = target.access_transcript('GENCODE_V39', 'ENST00000226091.3')
    assert str(trans) == "(GENCODE_V39, ENST00000226091.3) [('trans_cat', ''), ('trans_notes', ''), ('trans_structural_cat', '')]"

    trans = target.access_transcript('WTC11_consolidated', 'FSM_45093')
    assert str(trans) == "(WTC11_consolidated, FSM_45093) [('trans_cat', 'FSM'), ('trans_notes', ''), ('trans_structural_cat', 'PB-specific')]"


def test_primer_targets_bad_header():
    tsv = list(PRIMER_TARGETS_1)
    tsv[0] = tsv[0].replace("region_5p", "fred_5p")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuUserError) as exinfo:
        primer_targets_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, PrimersJuJuUserError,
                      "required column is missing: 'region_5p'")

def test_primer_targets_bad_coords():
    tsv = list(PRIMER_TARGETS_1)
    tsv[1] = tsv[1].replace("chr17:7709210-7710259", "chr17:77XZ210-7710259")
    tsv_fh = _mk_tsv(tsv)
    with pytest.raises(PrimersJuJuUserError) as exinfo:
        primer_targets_read("<test>", in_fh=tsv_fh)
    assert_except_msg(exinfo, CoordsError, "invalid coordinates: 'chr17:77XZ210-7710259'")
