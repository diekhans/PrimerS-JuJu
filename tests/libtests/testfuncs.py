"""
functions to support tests
"""
import os.path as osp
import pipettor

def test_id(request):
    return request.node.name

def diff_expected(rel_name):
    pipettor.run(["diff", "-u",
                  osp.join("expected", rel_name),
                  osp.join("output", rel_name)])
