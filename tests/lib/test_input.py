import os.path as osp
from primersjuju.primer_targets import primer_targets_read

root = osp.join(osp.dirname(__file__), "../..")



def test_primer_targets_example():
    targets = primer_targets_read(osp.join(root, "docs/primer-targets-example.tsv"))
