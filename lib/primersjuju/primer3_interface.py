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
        super().__init__(results)

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
    # positive transcript strand
    region_5p = target_transcript.region_5p.trans.abs()
    region_3p = target_transcript.region_3p.trans.abs()
    return [[region_5p.start + 1, len(region_5p),
             region_3p.start + 1, len(region_3p)]]

def _build_junction_overlap(features):
    # The junction associated with a given position is the space immediately
    # to the right of that position in the template sequence on the strand
    # given as input.
    assert len(features) in (1, 3)
    if len(features) == 3:
        return (features[1].trans.abs().start + 1,)
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
        # inside of primer regions and including primer regions, abs(), sorted() ignores strand.
        global_args.PRIMER_PRODUCT_SIZE_RANGE = [sorted([abs(target_transcript.region_5p.trans.start - target_transcript.region_3p.trans.end),
                                                         abs(target_transcript.region_5p.trans.end - target_transcript.region_3p.trans.start)])]
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
    for r in results:
        pp.pprint(vars(r))


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

def _make_point_char_inserts(point_list, mark_char):
    # primer3 point 1-based index position before junction
    return [(p - 1, mark_char) for p in point_list]

def _make_region_char_inserts(region_lists, start_char, end_char):
    # primer3 region 1-based index plus length, with multiple regions
    # in same list
    # [[1457, 230, 375, 717]]
    inserts = []
    for region_list in region_lists:
        for i in range(0, len(region_list), 2):
            start = region_list[i]
            end = start + region_list[i+1]
            inserts.append((start - 1, start_char))
            inserts.append((end, end_char))
    return inserts

def primer3_annotate_rna(target_transcript, *, global_args=_global_args_defaults):
    """Generate primer3-style annotated RNA, for debugging purposes """

    rna = target_transcript.rna
    seq_args = _build_seq_args(target_transcript)

    # build insert list [(before_index, char), ...]
    # include regions: {}, junctions: -
    inserts = sorted(_make_region_char_inserts(seq_args.SEQUENCE_PRIMER_PAIR_OK_REGION_LIST, '{', '}') +
                     _make_point_char_inserts(seq_args.SEQUENCE_OVERLAP_JUNCTION_LIST, '-'))
    # generate annotated sequence
    ann_rna_parts = []
    prev_end = 0
    for idx, char in inserts:
        ann_rna_parts.append(rna[prev_end:idx])
        ann_rna_parts.append(char)
        prev_end = idx
    ann_rna_parts.append(rna[prev_end:])
    return "".join(ann_rna_parts)
