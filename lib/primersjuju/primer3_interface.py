"""
Interface to primer3 python module
"""
import sys
import re
import copy
import pprint
import primer3
from pycbio.sys.objDict import ObjDict
from . import PrimersJuJuDataError, PrimersJuJuError

# Notes:
# forces the use of 0-based indexing

_global_args_defaults = ObjDict(PRIMER_TASK="generic",
                                PRIMER_FIRST_BASE_INDEX=0,
                                PRIMER_PICK_LEFT_PRIMER=1,
                                PRIMER_PICK_RIGHT_PRIMER=1,
                                PRIMER_EXPLAIN_FLAG=1)

class Primer3Pair(ObjDict):
    """One result pair set from primer3, files are primer three result names with
    number dropped.  So PRIMER_LEFT_0_GC_PERCENT becomes PRIMER_LEFT_GC_PERCENT.
    result_num has the result number.
    See https://www.primer3plus.com/primer3plusHelp.html#outputTags
    for tag description
    """
    def __init__(self, result_num):
        self.result_num = result_num

    def dump(self, fh=sys.stderr):
        pp = pprint.PrettyPrinter(stream=fh, sort_dicts=False, indent=4)
        print(f"Primer3Pair {self.result_num}:", file=fh)
        pp.pprint(self)


class Primer3Results(ObjDict):
    """results from primer3, contains an order Primer3Result objects,
    as well as all tags returned by primer3 as attributes"""

    def __init__(self, p3_output):
        # initialize with all result tags, include ones later split out into
        # Primer3Result objects
        for n in p3_output.keys():
            setattr(self, n, getattr(p3_output, n))
        self.pairs = []

    def dump(self, fh=sys.stderr):
        print(f">>> Primer3Results: {len(self.pairs)} results <<<", file=fh)
        for r in self.pairs:
            r.dump(fh)

def _parse_result(p3_output, result_num):
    "parse out the results for a give result number"
    pair_name_re = re.compile(f"^([A-Z_]+)_{result_num}(_[A-Z_]+)?$")
    pair = Primer3Pair(result_num)
    for fld, val in p3_output.items():
        m = pair_name_re.match(fld)
        if m is not None:
            name = m.group(1) + (m.group(2) if m.group(2) is not None else '')
            setattr(pair, name, val)
    return pair

def primer3_parse_output(p3_output):
    """convert primer3 results into a Primer3Results object"""
    if not isinstance(p3_output, ObjDict):
        p3_output = ObjDict(p3_output)
    results = Primer3Results(p3_output)
    if ((p3_output.PRIMER_LEFT_NUM_RETURNED != p3_output.PRIMER_PAIR_NUM_RETURNED) or
        (p3_output.PRIMER_RIGHT_NUM_RETURNED != p3_output.PRIMER_PAIR_NUM_RETURNED)):
        raise PrimersJuJuDataError(f"primer3 NUM_RETURN inconsistent, don't know how to deal with this: "
                                   f"PRIMER_LEFT_NUM_RETURNED ({p3_output.PRIMER_LEFT_NUM_RETURNED}) != "
                                   f"PRIMER_PAIR_NUM_RETURNED ({p3_output.PRIMER_PAIR_NUM_RETURNED}) != "
                                   f"PRIMER_RIGHT_NUM_RETURNED ({p3_output.PRIMER_RIGHT_NUM_RETURNED})")
    for result_num in range(p3_output.PRIMER_PAIR_NUM_RETURNED):
        results.pairs.append(_parse_result(p3_output, result_num))

    return results

def make_ok_region(target_transcript):
    # for SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
    # positive transcript strand, (0-based, length)
    region_5p = target_transcript.region_5p.trans.abs()
    region_3p = target_transcript.region_3p.trans.abs()
    return [[region_5p.start, len(region_5p),
             region_3p.start, len(region_3p)]]

def _build_junction_overlap(features):
    # For SEQUENCE_OVERLAP_JUNCTION_LIST, the junctions are specified as a position
    # to the right on the 1-based position in the sequence.  This is true even
    # if PRIMER_FIRST_BASE_INDEX=0
    assert len(features) in (1, 3)
    if len(features) == 3:
        return (features[1].trans.abs().start + 1,)
    else:
        return ()

def _build_strong_match_str(config):
    mstr = config.num_5_prime_strong_match * 'S'  # Strong (G or C)
    if config.num_5_prime_strong_match < config.PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION:
        # must cover the 5' junction overlap, so extend with `N' to match any
        mstr += (config.MIN_5_PRIME_OVERLAP_OF_JUNCTION - config.num_5_prime_strong_match) * 'N'
    return mstr

def _build_seq_args(config, target_transcript):
    ok_regions = make_ok_region(target_transcript)
    seq_args = ObjDict(SEQUENCE_ID=target_transcript.trans_id.name,
                       SEQUENCE_TEMPLATE=target_transcript.rna,
                       SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=ok_regions)
    junction_overlaps = _build_junction_overlap(target_transcript.features_5p) + _build_junction_overlap(target_transcript.features_3p)
    if len(junction_overlaps) > 0:
        seq_args.SEQUENCE_OVERLAP_JUNCTION_LIST = junction_overlaps

    return seq_args

def _compute_primer_product_size_range(target_transcript, global_args):
    # must be at least twice PRIMER_MAX_SIZE
    min_size = 2 * global_args.PRIMER_MAX_SIZE

    # inside of primer regions and including primer regions, uses abs() and sort() to ignore strand.
    size_ranges = [
        max(abs(target_transcript.region_5p.trans.start - target_transcript.region_3p.trans.end),
            min_size),
        max(abs(target_transcript.region_5p.trans.end - target_transcript.region_3p.trans.start),
            min_size)]
    size_ranges.sort()
    return [size_ranges]

def _copy_config_global_args(config, global_args):
    for fld in dir(config):
        if fld.startswith("PRIMER_"):
            val = getattr(config, fld)
            if val is not None:
                global_args[fld] = val

def _build_global_args(config, target_transcript):
    # set PRIMER_PRODUCT_SIZE_RANGE for specific sequence
    global_args = copy.copy(_global_args_defaults)
    _copy_config_global_args(config, global_args)

    # likes oligos 5' with G or C will anneal with 3' of sequence
    # S = Strong (G or C)
    if config.num_5_prime_strong_match > 0:
        global_args.PRIMER_MUST_MATCH_FIVE_PRIME = _build_strong_match_str(config)

    global_args.PRIMER_PRODUCT_SIZE_RANGE = _compute_primer_product_size_range(target_transcript, global_args)

    return global_args

def _check_common_errors(target_transcript, seq_args, global_args):
    """check for common errors in specification that would generate primer3
    errors and generate more useful error message."""

    def _check_region(region_features, size_desc, min_size):
        exon_size = len(region_features.bounds.trans)
        if exon_size < min_size:
            genome_region = region_features.bounds.genome
            raise PrimersJuJuDataError(f"transcript {target_transcript} region {genome_region} exon length {exon_size} is less than {size_desc} {min_size}")

    for region_features in (target_transcript.features_5p, target_transcript.features_3p):
        _check_region(region_features, "PRIMER_MIN_SIZE", global_args.PRIMER_MIN_SIZE)
        _check_region(region_features, "PRIMER_MAX_SIZE", global_args.PRIMER_MAX_SIZE)

def primer3_design(primer3_config, target_transcript, *, debug=False):
    """main entry to run primer3
    global_args defined here:
    https://www.primer3plus.com/primer3plusHelp.html#globalTags

    global PRIMER_FIRST_BASE_INDEX must be zero.
    """

    global_args = _build_global_args(primer3_config, target_transcript)
    seq_args = _build_seq_args(primer3_config, target_transcript)
    if debug:
        print(">>>>> primer3 debug:", target_transcript, file=sys.stderr)
        primer3_dump_args(sys.stderr, primer3_config, target_transcript,
                          global_args=global_args, seq_args=seq_args)

    _check_common_errors(target_transcript, seq_args, global_args)

    # FIXME: dict() is a hack around primer3 give type error on ObjDict
    p3_output = primer3.bindings.design_primers(dict(seq_args), dict(global_args),
                                                misprime_lib=primer3_config.misprime_lib,
                                                mishyb_lib=primer3_config.mishyb_lib)
    results = primer3_parse_output(p3_output)
    if "PRIMER_ERRORS" in results:
        raise PrimersJuJuError(f"primer3 errors: {results.PRIMER_ERRORS}")
    if "PRIMER_WARNINGS" in results:
        raise PrimersJuJuError(f"primer3 warnings: {results.PRIMER_WARNINGS}")
    return results

def primer3_dump_args(fh, primer3_config, target_transcript, *, global_args=None, seq_args=None):
    "print the arguments that will be used for this design"
    if global_args is None:
        global_args = _build_global_args(primer3_config, target_transcript)
    if seq_args is None:
        seq_args = _build_seq_args(primer3_config, target_transcript)
    pp = pprint.PrettyPrinter(stream=fh, sort_dicts=False, indent=4)
    print(">>> primer3_args <<<", file=fh)
    print("global_args:", file=fh)
    pp.pprint(global_args)
    print("seq_args:", file=fh)
    pp.pprint(seq_args)

def _make_point_char_inserts(points, mark_char):
    # primer3 point 0-based index position before junction
    return tuple((p, mark_char) for p in points)

def _make_region_char_inserts(region_range, start_char, end_char):
    # region_range is (start, end)
    return ((region_range[0], start_char),
            (region_range[1], end_char))

def primer3_annotate_amplicon(primer3_config, target_transcript):
    """Generate primer3 web annotated amplicon, for debugging purposes"""
    # Internally we use SEQUENCE_PRIMER_PAIR_OK_REGION_LIST to define regions to
    # design primers.  However, this isn't support by primer3 annotations, instead we
    # define:
    #   Targets: [ ] Region inside of specified primer regions
    # Base on SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
    # Targeted junctions are marked with: '-'
    #
    # None of the approaches tried to excluding the regions before and after the primer
    # target regions work, so we output the amplicaton and exclude the before and after
    # regions
    rna = target_transcript.rna
    seq_args = _build_seq_args(primer3_config, target_transcript)

    # build regions in primer3 style, first convert to [start, end]
    ok_region_ranges = []
    for ok_regions in seq_args.SEQUENCE_PRIMER_PAIR_OK_REGION_LIST:
        ok_region_ranges += [(ok_regions[i], ok_regions[i] + ok_regions[i + 1])
                             for i in range(0, len(ok_regions) - 1, 2)]
    assert len(ok_region_ranges) == 2
    # regions we consider
    before = (0, ok_region_ranges[0][0])
    target = (ok_region_ranges[0][1], ok_region_ranges[1][0])
    after = (ok_region_ranges[1][1], target_transcript.region_5p.trans.size)

    # Build insert list [(before_index, char), ...], sort for insertion order.  transcript must cover
    # the whole range
    inserts = list(sorted(_make_region_char_inserts(target, '[', ']') +
                          _make_point_char_inserts(seq_args.get("SEQUENCE_OVERLAP_JUNCTION_LIST", ()), '-')))
    # generate annotated sequence
    amplicon_parts = []
    prev_end = before[1]
    # add target and junctions
    for idx, char in inserts:
        amplicon_parts.append(rna[prev_end:idx])
        amplicon_parts.append(char)
        prev_end = idx
    amplicon_parts.append(rna[prev_end:after[0]])
    return "".join(amplicon_parts)
