"""
Interface to primer3 python module
"""
import re
import copy
import pprint
import primer3
from pycbio.sys.objDict import ObjDict
from . import PrimersJuJuDataError

_global_args_defaults = ObjDict(PRIMER_TASK="generic",
                                PRIMER_PICK_LEFT_PRIMER=1,
                                PRIMER_PICK_INTERNAL_OLIGO=0,
                                PRIMER_PICK_RIGHT_PRIMER=1,
                                PRIMER_OPT_SIZE=20,
                                PRIMER_MIN_SIZE=18,
                                PRIMER_MAX_SIZE=22,
                                PRIMER_EXPLAIN_FLAG=1)

# these might me configured in the future
FIVE_PRIME_NUM_STRONG_MATCH = 2

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

def parse_results(results):
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

def make_ok_region(target_transcript):
    return [[target_transcript.region_3p.trans.start + 1,
             len(target_transcript.region_3p.trans),
             target_transcript.region_5p.trans.start + 1,
             len(target_transcript.region_5p.trans)]]

def _build_junction_overlap(features):
    # The junction associated with a given position is the space immediately
    # to the right of that position in the template sequence on the strand
    # given as input.
    if len(features) == 3:
        return (features[1].trans.start + 1,)
    else:
        return ()

def _build_seq_args(target_transcript):
    # likes oligos 5' with G or C will anneal with 3' of sequence
    # S = Strong (G or C)
    match_5p = FIVE_PRIME_NUM_STRONG_MATCH * 'S'

    ok_regions = make_ok_region(target_transcript)
    junction_overlaps = _build_junction_overlap(target_transcript.features_5p.features) + _build_junction_overlap(target_transcript.features_3p.features)
    seq_args = ObjDict(SEQUENCE_ID=target_transcript.trans_id,
                       SEQUENCE_TEMPLATE=target_transcript.rna,
                       SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=ok_regions,
                       PRIMER_MUST_MATCH_FIVE_PRIME=match_5p,
                       SEQUENCE_OVERLAP_JUNCTION_LIST=junction_overlaps)
    return seq_args

def _build_global_args(target_transcript, global_args):
    # set PRIMER_PRODUCT_SIZE_RANGE for specific sequence
    if "PRIMER_PRODUCT_SIZE_RANGE" not in global_args:
        global_args = copy.copy(global_args)
        # inside of primer regions and including primer regions
        global_args.PRIMER_PRODUCT_SIZE_RANGE = [[target_transcript.region_5p.trans.start - target_transcript.region_3p.trans.end,
                                                  target_transcript.region_5p.trans.end - target_transcript.region_3p.trans.start]]
    return global_args

def primer3_global_defaults():
    "get a copy to modify"
    return copy.copy(_global_args_defaults)

def _dump_primer3_info(target_transcript, global_args, seq_args, results, dump_fh):
    pp = pprint.PrettyPrinter(stream=dump_fh, sort_dicts=False, indent=4)
    print("global_args:", file=dump_fh)
    pp.pprint(global_args)
    print("seq_args:", file=dump_fh)
    pp.pprint(seq_args)
    print("results:", file=dump_fh)
    pp.pprint(results)


def primer3_design(target_transcript, *, global_args=_global_args_defaults, misprime_lib=None, mishyb_lib=None, debug=False,
                   dump_fh=None):
    "main entry to run primer3"

    global_args = _build_global_args(target_transcript, global_args)
    seq_args = _build_seq_args(target_transcript)

    results = parse_results(primer3.bindings.designPrimers(seq_args, global_args,
                                                           misprime_lib=misprime_lib, mishyb_lib=mishyb_lib,
                                                           debug=debug))
    if dump_fh is not None:
        _dump_primer3_info(target_transcript, global_args, seq_args, results, dump_fh)
    return results
