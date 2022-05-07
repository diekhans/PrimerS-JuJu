import os.path as osp
import io
from primersjuju.primer_targets import primer_targets_read

mydir = osp.join(osp.dirname(__file__))

def test_primer_targets_example():
    targets = primer_targets_read(osp.join(mydir, "../../docs/primer-targets-example.tsv"))

    target_ids = sorted(targets.targets.keys())
    assert target_ids == ["EFNB3+1"]

    target = targets.access_target("EFNB3+1")
    trans_ids = sorted(target.get_tracks_trans())
    assert trans_ids == [('GENCODE_V39', 'ENST00000226091.3'),
                         ('WTC11_consolidated', 'FSM_45093'),
                         ('WTC11_consolidated',  'NNC-318304')]
    trans = target.access_transcript('GENCODE_V39', 'ENST00000226091.3')
    assert str(trans) == "(GENCODE_V39, ENST00000226091.3) [('trans_cat', ''), ('trans_notes', ''), ('trans_structural_cat', '')]"

    trans = target.access_transcript('WTC11_consolidated', 'FSM_45093')
    assert str(trans) == "(WTC11_consolidated, FSM_45093) [('trans_cat', 'FSM'), ('trans_notes', ''), ('trans_structural_cat', 'PB-specific')]"
