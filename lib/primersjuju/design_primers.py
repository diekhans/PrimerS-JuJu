"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuError
from .primer3_interface import Primer3Results, Primer3Pair, primer3_design
from .primer_targets import PrimerTargets, TargetTranscript
from .transcript_features import Features, transcript_range_to_features, features_to_transcript_coords, features_to_genomic_coords
from .primer_uniqueness import PrimerUniqueness, primer_uniqueness_query, primer_uniqueness_none

class DesignStatus(SymEnum):
    """Status of the design for a given target.  Smaller is better"""
    GOOD = 0
    NOT_GENOME_UNIQUE = 1
    NOT_TRANSCRIPTOME_UNIQUE = 2
    NO_PRIMERS = 3

@dataclass
class PrimerDesign:
    """information collected on one primer pair"""
    ppair_id: str
    primer3_pair: Primer3Pair
    features_5p: Features  # in transcription order, positive genomic strand
    features_3p: Features
    amplicon_coords: Coords
    uniqueness: PrimerUniqueness
    # ranking by stability and uniqueness
    priority: int = None

    def spans_splice_juncs(self):
        return (len(self.features_5p) > 1) or (len(self.features_3p) > 1)

    @property
    def amplicon_length(self):
        return len(self.amplicon_coords)

    def dump(self, fh):
        def _lfmt(l):
            if l is None:
                return str(l)
            else:
                return '[' + ",\n\t".join([str(v) for v in l]) + '\n    ]'

        def _print_p3_attr(attr):
            print("   ", attr, self.primer3_pair[attr], file=fh)
        print(">>> PrimerDesign <<<", file=fh)
        print("    ppair_id", self.ppair_id, file=fh)
        print("    features_5p", self.features_5p, file=fh)
        print("    features_3p", self.features_3p, file=fh)
        print("    priority", self.priority, file=fh)
        print("    amplicon_coords", self.amplicon_coords, file=fh)
        print("    amplicon_length", self.amplicon_length, file=fh)
        _print_p3_attr("PRIMER_LEFT")
        _print_p3_attr("PRIMER_RIGHT")
        _print_p3_attr("PRIMER_LEFT_SEQUENCE")
        _print_p3_attr("PRIMER_RIGHT_SEQUENCE")
        _print_p3_attr("PRIMER_PAIR_PRODUCT_SIZE")
        print("    genome_on_targets", _lfmt(self.uniqueness.genome_on_targets), file=fh)
        print("    genome_off_targets", _lfmt(self.uniqueness.genome_off_targets), file=fh)
        print("    genome_non_targets", _lfmt(self.uniqueness.genome_non_targets), file=fh)
        print("    transcriptome_on_targets", _lfmt(self.uniqueness.transcriptome_on_targets), file=fh)
        print("    transcriptome_off_targets", _lfmt(self.uniqueness.transcriptome_off_targets), file=fh)
        print("    transcriptome_non_targets", _lfmt(self.uniqueness.transcriptome_non_targets), file=fh)


@dataclass
class PrimerDesigns:
    """information collected on all primers for a target"""
    target_id: str
    primer_targets: PrimerTargets
    target_transcript: TargetTranscript
    primer3_results: Primer3Results
    uniqueness_checked: bool
    designs: Sequence[PrimerDesign]
    status: DesignStatus

    def dump(self, fh):
        print(f">>> PrimerDesigns {len(self.designs)} <<<", file=fh)
        print("    target_id", self.target_id, file=fh)
        print("    tran_id", self.target_transcript.trans_id, file=fh)
        for design in self.designs:
            design.dump(fh)

def _get_exon_features(target_transcript, start, end):
    # start/end are transcript forward direction
    feat0 = target_transcript.features[0]
    primer_tcoords = Coords(feat0.trans.name, start, end,
                            strand='+', size=feat0.trans.size)
    return transcript_range_to_features(target_transcript.features, primer_tcoords)

def _get_exon_left_features(target_transcript, primer3_tcoords):
    start = primer3_tcoords[0]
    end = start + primer3_tcoords[1]
    return _get_exon_features(target_transcript, start, end)

def _get_exon_right_features(target_transcript, primer3_tcoords):
    # primer3 has the zero-based end position of reverse primer in RNA
    end = primer3_tcoords[0] + 1
    start = end - primer3_tcoords[1]
    return _get_exon_features(target_transcript, start, end)

def _validate_primer_features(features_5p, features_3p):
    for feature_5p in features_5p:
        for feature_3p in features_3p:
            if feature_5p.trans.overlaps(feature_3p.trans):
                raise PrimersJuJuError(f"primer3 pairs overlap in transcript space {features_5p} and {features_3p}")
            if feature_5p.genome.overlaps(feature_3p.genome):
                raise PrimersJuJuError(f"primer3 pairs overlap in genome space {features_5p} and {features_3p}")

def _build_primer_design(target_transcript, target_id, result_num, primer3_pair, uniqueness_query):
    ppair_id = "{}_pp{}".format(target_id, result_num)
    features_5p = _get_exon_left_features(target_transcript, primer3_pair.PRIMER_LEFT)
    features_3p = _get_exon_right_features(target_transcript, primer3_pair.PRIMER_RIGHT)
    _validate_primer_features(features_5p, features_3p)
    amplicon_coords = features_to_transcript_coords(features_5p + features_3p)
    assert len(amplicon_coords) == primer3_pair.PRIMER_PAIR_PRODUCT_SIZE

    if uniqueness_query is not None:
        uniqueness = primer_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair)
    else:
        uniqueness = primer_uniqueness_none()

    return PrimerDesign(ppair_id, primer3_pair, features_5p, features_3p, amplicon_coords, uniqueness)

def _calc_design_status(primer_design) -> DesignStatus:
    if primer_design.uniqueness.transcriptome_off_target_cnt > 0:
        return DesignStatus.NOT_TRANSCRIPTOME_UNIQUE
    elif primer_design.uniqueness.genome_off_target_cnt > 0:
        return DesignStatus.NOT_GENOME_UNIQUE
    else:
        return DesignStatus.GOOD

def _get_design_status(primer_design_list) -> DesignStatus:
    "determine status base on primer3 and uniqueness"
    design_status = DesignStatus.NO_PRIMERS  # worst
    for primer_design in primer_design_list:
        design_status = min(design_status, _calc_design_status(primer_design))
    return design_status

def _primer_design_target_score(primer_design):
    # order of prority:
    #   on_target only or no alignments in it has introns (might not be in transcriptome)
    #   no off-target, but can have non-target or no alignments
    #   off-target
    on_cnt = off_cnt = non_cnt = 0
    uniqueness = primer_design.uniqueness
    on_cnt = uniqueness.genome_on_target_cnt + uniqueness.transcriptome_on_target_cnt
    off_cnt = uniqueness.genome_off_target_cnt + uniqueness.transcriptome_off_target_cnt
    non_cnt = uniqueness.genome_non_target_cnt + uniqueness.transcriptome_non_target_cnt
    if ((off_cnt == 0) and (non_cnt == 0)
        and ((primer_design.spans_splice_juncs() and on_cnt == 0) or (on_cnt > 0))):
        return 1
    elif off_cnt == 0:
        return 2
    else:
        return 3

def _primer_design_sort_key(primer_design):
    """lowest values are best"""
    # higher delta-G is better, take total
    return (_primer_design_target_score(primer_design),
            -(primer_design.primer3_pair.PRIMER_LEFT_END_STABILITY +
              primer_design.primer3_pair.PRIMER_RIGHT_END_STABILITY))

def _sort_primer_designs(primer_design_list):
    "sort and assign priorities"
    primer_design_list = sorted(primer_design_list, key=_primer_design_sort_key)
    for i in range(len(primer_design_list)):
        primer_design_list[i].priority = i + 1
    return primer_design_list

def _build_primer_designs(primer_targets, target_transcript, primer3_results, uniqueness_query):
    primer_design_list = [_build_primer_design(target_transcript, primer_targets.target_id, i + 1, pair, uniqueness_query)
                          for i, pair in enumerate(primer3_results.pairs)]
    primer_design_list = _sort_primer_designs(primer_design_list)
    return PrimerDesigns(primer_targets.target_id, primer_targets, target_transcript, primer3_results,
                         uniqueness_query is not None,
                         primer_design_list, _get_design_status(primer_design_list))

def design_primers(genome_data, primer_targets, *, uniqueness_query=None, primer3_debug=False,
                   amplicon_isoform_query=None):
    """design transcripts """
    target_transcript = primer_targets.transcripts[0]
    primer3_results = primer3_design(target_transcript, debug=primer3_debug)

    return _build_primer_designs(primer_targets, target_transcript, primer3_results, uniqueness_query)

def primer_design_amplicon_coords(primer_design, target_transcript):
    """return amplicon coordinates for an RNA and primer design """

    gcoords = features_to_genomic_coords(primer_design.features_5p + primer_design.features_3p)
    amp_features = target_transcript.features.intersect_genome(gcoords)
    return features_to_transcript_coords(amp_features)

def primer_design_amplicon(primer_design, target_transcript):
    """return amplicon for an RNA and primer design """
    tcoords = primer_design_amplicon_coords(primer_design, target_transcript)
    return target_transcript.rna[tcoords.start:tcoords.end]
