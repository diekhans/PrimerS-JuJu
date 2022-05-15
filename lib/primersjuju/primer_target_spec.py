"""
Primer targets specification from user.  Including Parsing and validated
primer_targets TSV input file.

This does
"""
import re
from pycbio.tsv import TsvReader
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuDataError

# for validation
REGION_COLS = frozenset(["region_5p", "region_3p"])
TRANSCRIPT_COLS = frozenset(["trans_track", "trans_id"])
REQUIRED_COLS = frozenset(frozenset(["target_id"]) | REGION_COLS | TRANSCRIPT_COLS)

class TargetTranscriptSpec:
    """a target transcript"""
    def __init__(self, trans_track, trans_id, user_attrs):
        self.trans_track = trans_track
        self.trans_id = trans_id
        self.user_attrs = user_attrs

    def __str__(self):
        ua = [(k, self.user_attrs[k]) for k in sorted(self.user_attrs.keys())]
        return f"({self.trans_track}, {self.trans_id}) {ua}"

class PrimerTargetSpec:
    """
    A given primer target regions and associate target transcripts
    """
    def __init__(self, target_id, region_5p, region_3p, user_attrs):
        self.target_id = target_id
        self.region_5p, self.region_3p = region_5p, region_3p
        self.user_attrs = user_attrs
        self.tracks = {}  # by [trans_track][trans_id]

    def add_transcript(self, trans_track, trans_id, user_attrs):
        track = self.tracks.get(trans_track)
        if track is None:
            track = self.tracks[trans_track] = {}
        if trans_id in track:
            raise PrimersJuJuDataError(f"duplicate transcript for primer: ({trans_track}, {trans_id})")
        trans = TargetTranscriptSpec(trans_track, trans_id, user_attrs)
        self.tracks[trans_track][trans_id] = trans
        return trans

    def access_transcript(self, trans_track, trans_id):
        "transcript or error"
        try:
            return self.tracks[trans_track][trans_id]
        except KeyError:
            raise PrimersJuJuDataError(f"unknown transcript ({trans_track}, {trans_id})")

    def get_tracks_trans(self):
        "returns list of tuples of (trans_track, trans_id)"
        return [(tr, tid)
                for tr in self.tracks.keys()
                for tid in self.tracks[tr].keys()]


class PrimerTargetSpecs:
    """All specified primer target transcripts, each with a unique id"""
    def __init__(self):
        self.targets = {}

    def add_target(self, target_id, region_5p, region_3p, user_attrs):
        if target_id in self.targets:
            raise PrimersJuJuDataError(f"duplicate primer target_id '{target_id}'")
        target = PrimerTargetSpec(target_id, region_5p, region_3p, user_attrs)
        self.targets[target_id] = target
        return target

    def get_target_ids(self):
        return list(self.targets.keys())

    def get_target(self, target_id):
        "target or None"
        return self.targets.get(target_id)

    def access_target(self, target_id):
        "target or error"
        target = self.targets.get(target_id)
        if target is None:
            raise PrimersJuJuDataError(f"unknown primer target_id references '{target_id}")
        return target

def _check_target_id(target_id):
    if not re.match("^[A-Za-z][-_.=%+A-Za-z0-9]*$", target_id):
        raise PrimersJuJuDataError(f"invalid target_id string '{target_id}', see documentation")

def _must_not_be_empty(columns, row):
    for col in columns:
        if row[col] == "":
            raise PrimersJuJuDataError(f"row column {col} must not be empty")

def _must_be_empty(columns, row):
    for col in columns:
        if row[col] != "":
            raise PrimersJuJuDataError(f"row column {col} must be empty")

def _parse_coords(coord_str):
    coords = Coords.parse(coord_str, oneBased=True)
    if len(coords) > 10000000:
        raise PrimersJuJuDataError(f"coordinates seem absurdly long: '{coord_str}'")
    return coords

def _get_user_cols(rows):
    row0 = rows[0]
    target_user_cols = set()
    transcript_user_cols = set()
    for col in row0._columns_:
        if col not in REQUIRED_COLS:
            if col.startswith("trans_"):
                transcript_user_cols.add(col)
            else:
                target_user_cols.add(col)
    return target_user_cols, transcript_user_cols

def _build_column_dict(columns, row):
    return {col: row[col] for col in columns}

def _do_add_primary_row(primer_target_specs, target_user_cols, transcript_user_cols, row):
    _check_target_id(row.target_id)
    _must_not_be_empty(REQUIRED_COLS, row)

    region_5p = _parse_coords(row.region_5p)
    region_3p = _parse_coords(row.region_3p)
    if region_3p.overlaps(region_5p):
        raise PrimersJuJuDataError(f"region_5p ({region_5p}) overlaps region_3p ({region_3p})")
    if region_3p.name != region_5p.name:
        raise PrimersJuJuDataError(f"region_5p ({region_5p}) is on a different chromosome than region_3p ({region_3p})")

    target = primer_target_specs.add_target(row.target_id, region_5p, region_3p,
                                            _build_column_dict(target_user_cols, row))
    target.add_transcript(row.trans_track, row.trans_id,
                          _build_column_dict(transcript_user_cols, row))

def _add_primary_row(primer_target_specs, target_user_cols, transcript_user_cols, row):
    try:
        _do_add_primary_row(primer_target_specs, target_user_cols, transcript_user_cols, row)
    except Exception as ex:
        raise PrimersJuJuDataError(f"error parsing primary row: {str(row)}") from ex

def _do_add_continue_row(primer_target_specs, transcript_user_cols, row):
    _check_target_id(row.target_id)
    _must_be_empty(REGION_COLS, row)
    _must_not_be_empty(TRANSCRIPT_COLS, row)
    target = primer_target_specs.access_target(row.target_id)
    target.add_transcript(row.trans_track, row.trans_id,
                          _build_column_dict(transcript_user_cols, row))

def _add_continue_row(primer_target_specs, transcript_user_cols, row):
    try:
        _do_add_continue_row(primer_target_specs, transcript_user_cols, row)
    except Exception as ex:
        raise PrimersJuJuDataError(f"error parsing continuation row: {str(row)}") from ex

def _primer_target_specs_build(rows):
    # since no order required, two passes; one to for primary data and the
    # other for continuation
    primer_target_specs = PrimerTargetSpecs()
    target_user_cols, transcript_user_cols = _get_user_cols(rows)
    for row in rows:
        if row.region_5p != "":
            _add_primary_row(primer_target_specs, target_user_cols, transcript_user_cols, row)
    for row in rows:
        if row.region_5p == "":
            _add_continue_row(primer_target_specs, transcript_user_cols, row)
    return primer_target_specs

def _check_required_columns(rows):
    if len(rows) == 0:
        raise PrimersJuJuDataError("no data in TSV")
    row0 = rows[0]
    for col in REQUIRED_COLS:
        if col not in row0:
            raise PrimersJuJuDataError(f"required column is missing: '{col}'")

def primer_targets_specs_read(primer_targets_tsv, in_fh=None):
    """read all primer targets into PrimerTargets object"""
    try:
        rows = [row for row in TsvReader(primer_targets_tsv, inFh=in_fh)]
        _check_required_columns(rows)
        return _primer_target_specs_build(rows)
    except Exception as ex:
        raise PrimersJuJuDataError(f"error parsing primary target specification TSV: '{primer_targets_tsv}'") from ex
