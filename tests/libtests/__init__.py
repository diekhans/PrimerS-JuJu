import re
import os.path as osp

mydir = osp.join(osp.dirname(__file__))

def assert_except_msg(exinfo, cls, *, msgre=None, msg=None):
    "check causes"
    assert (msgre is not None) or (msg is not None)

    ex = exinfo.value
    while ex is not None:
        if (isinstance(ex, cls) and
            ((msgre is not None) and re.search(msgre, str(ex))) or
            ((msg is not None) and (str(ex).find(msg) >= 0))):
            return
        ex = ex.__cause__
    raise Exception(f"Caught exception can't match: {str(cls)} {msgre}") from exinfo.value
