"""
Interface to primer3 python module
"""
import re
from pycbio.sys.objDict import ObjDict
from . import PrimersJuJuDataError

_global_args_defaults = ObjDict()
_global_args_defaults.PRIMER_TASK = "generic"
_global_args_defaults.PRIMER_PICK_LEFT_PRIMER = 1
_global_args_defaults.PRIMER_PICK_INTERNAL_OLIGO = 0
_global_args_defaults.PRIMER_PICK_RIGHT_PRIMER = 1
_global_args_defaults.PRIMER_OPT_SIZE = 20
_global_args_defaults.PRIMER_MIN_SIZE = 18
_global_args_defaults.PRIMER_MAX_SIZE = 22
_global_args_defaults.PRIMER_PRODUCT_SIZE_RANGE = [[75, 150]]
_global_args_defaults.PRIMER_EXPLAIN_FLAG = 1


class Primer3Result:
    """One result set from primer3, files are primer three result names with
    number dropped.  So PRIMER_LEFT_0_GC_PERCENT becomes PRIMER_LEFT_GC_PERCENT.
    result_num has the result number"""
    def __init__(self, result_num):
        self.result_num = result_num

class Primer3Results(list):
    """results from primer3, an order list of Primer3Result objects"""
    def __init__(self, results):
        self.PRIMER_LEFT_EXPLAIN = results.PRIMER_LEFT_EXPLAIN
        self.PRIMER_PAIR_EXPLAIN = results.PRIMER_PAIR_EXPLAIN
        self.PRIMER_RIGHT_EXPLAIN = results.PRIMER_RIGHT_EXPLAIN

def _parse_result(results, result_num):
    result_name_re = re.compile(f"^([A-Z_]+)_{result_num}(_[A-Z_]+)?$")
    result = Primer3Result(result_num)
    for fld, val in results.items():
        m = result_name_re.match(fld)
        if m is not None:
            name = m.group(1) + (m.group(2) if m.group(2) is not None else '')
            setattr(result, name, val)
    return result

def primer3_parse_results(results):
    """convert primer3 results into a Primer3Results object"""
    if not isinstance(results, ObjDict):
        results = ObjDict(results)
    primer3_results = []
    if ((results.PRIMER_LEFT_NUM_RETURNED != results.PRIMER_PAIR_NUM_RETURNED) or
        (results.PRIMER_RIGHT_NUM_RETURNED != results.PRIMER_PAIR_NUM_RETURNED)):
        raise PrimersJuJuDataError(f"primer3 NUM_RETURN inconsistent, don't know how to deal with this: "
                                   f"PRIMER_LEFT_NUM_RETURNED ({results.PRIMER_LEFT_NUM_RETURNED}) != "
                                   f"PRIMER_PAIR_NUM_RETURNED ({results.PRIMER_PAIR_NUM_RETURNED}) != "
                                   f"PRIMER_RIGHT_NUM_RETURNED ({results.PRIMER_RIGHT_NUM_RETURNED})")
    for result_num in range(results.PRIMER_PAIR_NUM_RETURNED):
        primer3_results.append(_parse_result(results, result_num))

    return primer3_results

def primer3_design(primer_target):
    pass
