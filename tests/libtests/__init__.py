import re
import os.path as osp

mydir = osp.join(osp.dirname(__file__))

def assert_except_msg(exinfo, cls, msgre):
    "check causes"
    ex = exinfo.value
    while ex is not None:
        if isinstance(ex, cls) and re.match(msgre, str(ex)):
            return
        ex = ex.__cause__
    raise Exception(f"Caught exception can't match: {str(cls)} {msgre}") from exinfo.value
