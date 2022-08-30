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
from .transcript_features import Features, ExonFeature, transcript_range_to_features
from .uniqueness_query import GenomeHit, TranscriptomeHit

class DesignStatus(SymEnum):
    """Status of the design for a given target.  Smaller is better"""
    GOOD = 0
    NOT_GENOME_UNIQUE = 1
    NOT_TRANSCRIPTOME_UNIQUE = 2
    NO_PRIMERS = 3

def _len_none(l):
    "0 if list is None, else length"
    return 0 if l is None else len(l)

@dataclass
class PrimerDesign:
    """information collected on one primer pair"""
    ppair_id: str
    primer3_pair: Primer3Pair
    features_5p: Features  # in transcription order, positive genomic strand
    features_3p: Features
    # results from uniqueness query, non_target are ones on sequences liked fixes and alts
    # None if uniqueness check was not done.
    genome_on_targets: Sequence[GenomeHit]
    genome_off_targets: Sequence[GenomeHit]
    genome_non_targets: Sequence[GenomeHit]
    # results from transcriptome query
    transcriptome_on_targets: Sequence[TranscriptomeHit]
    transcriptome_off_targets: Sequence[TranscriptomeHit]
    transcriptome_non_targets: Sequence[TranscriptomeHit]
    priority: int = None

    def amplicon_trans_coords(self) -> Coords:
        """amplicon region, in positive transcript coordinates """
        trans_5p = self.features_5p[0].trans
        trans_3p = self.features_3p[0].trans
        if trans_5p.strand == '+':
            return Coords(trans_5p.name, trans_5p.start, trans_3p.end, '+', trans_5p.size)
        else:
            return Coords(trans_5p.name, trans_5p.size - trans_5p.end, trans_5p.size - trans_3p.start, '+', trans_5p.size)

    def spans_splice_juncs(self):
        return (len(self.features_5p) > 1) or (len(self.features_3p) > 1)

    @property
    def amplicon_length(self):
        # ignore strand with abs
        return abs(self.features_3p[-1].trans.end - self.features_5p[0].trans.start)

    @property
    def genome_on_target_cnt(self):
        return _len_none(self.genome_on_targets)

    @property
    def genome_off_target_cnt(self):
        return _len_none(self.genome_off_targets)

    @property
    def genome_non_target_cnt(self):
        return _len_none(self.genome_non_targets)

    @property
    def transcriptome_on_target_cnt(self):
        return _len_none(self.transcriptome_on_targets)

    @property
    def transcriptome_off_target_cnt(self):
        return _len_none(self.transcriptome_off_targets)

    @property
    def transcriptome_non_target_cnt(self):
        return _len_none(self.transcriptome_non_targets)

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
        print("    amplicon_trans_coords", self.amplicon_trans_coords(), file=fh)
        print("    amplicon_length", self.amplicon_length, file=fh)
        _print_p3_attr("PRIMER_LEFT")
        _print_p3_attr("PRIMER_RIGHT")
        _print_p3_attr("PRIMER_LEFT_SEQUENCE")
        _print_p3_attr("PRIMER_RIGHT_SEQUENCE")
        _print_p3_attr("PRIMER_PAIR_PRODUCT_SIZE")
        print("    genome_on_targets", _lfmt(self.genome_on_targets), file=fh)
        print("    genome_off_targets", _lfmt(self.genome_off_targets), file=fh)
        print("    genome_non_targets", _lfmt(self.genome_non_targets), file=fh)
        print("    transcriptome_on_targets", _lfmt(self.transcriptome_on_targets), file=fh)
        print("    transcriptome_off_targets", _lfmt(self.transcriptome_off_targets), file=fh)
        print("    transcriptome_non_targets", _lfmt(self.transcriptome_non_targets), file=fh)


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
        print("    tran_id", self.target_transcript.track_name, self.target_transcript.trans_id, file=fh)
        for design in self.designs:
            design.dump(fh)

def _get_exon_features(target_transcript, start, end):
    # start/end are transcript forward direction
    feat0 = target_transcript.features[0]
    primer_trans_coords = Coords(feat0.trans.name, start, end,
                                 strand='+', size=feat0.trans.size)
    return transcript_range_to_features(target_transcript.features, primer_trans_coords)

def _get_exon_left_features(target_transcript, primer3_coords):
    start = primer3_coords[0]
    end = start + primer3_coords[1]
    return _get_exon_features(target_transcript, start, end)

def _get_exon_right_features(target_transcript, primer3_coords):
    start = primer3_coords[0] - primer3_coords[1]
    end = start + primer3_coords[1]
    return _get_exon_features(target_transcript, start, end)

def _validate_primer_features(features_5p, features_3p):
    for feature_5p in features_5p:
        for feature_3p in features_3p:
            if feature_5p.trans.overlaps(feature_3p.trans):
                raise PrimersJuJuError(f"primer3 pairs overlap in transcript space {features_5p} and {features_3p}")
            if feature_5p.genome.overlaps(feature_3p.genome):
                raise PrimersJuJuError(f"primer3 pairs overlap in genome space {features_5p} and {features_3p}")

def _is_target_chrom(genome_data, chrom_name):
    """Should a chromosome be consider at target for off-target/on target
    check.  This excludes patches and alts.  For genomes without this information,
    all sequences are conserved a target.
    """
    if genome_data.assembly_info is None:
        # no information so consider this a possible target
        return True
    chrom_info = genome_data.assembly_info.getByName(chrom_name)
    return ((chrom_info.sequenceRole == "assembled-molecule") and
            (chrom_info.assemblyUnit == "Primary Assembly"))

def _check_coords_overlap(features, coordses):
    for f in features.iter_type(ExonFeature):
        for coords in coordses:
            if f.genome.overlaps(coords):
                return True
    return False

def _check_hit_overlap(target_transcript, left_coordses, right_coordses):
    """check overlap based on list of genomic coordinates"""
    features_first, features_last = target_transcript.get_genome_ordered_features()
    return (_check_coords_overlap(features_first, left_coordses) and
            _check_coords_overlap(features_last, right_coordses))

def _check_genome_hit_overlap(target_transcript, hit):
    """does a genome uniqueness hit correspond to the target regions"""
    # GENOME hits will not span introns, but may partially align one
    # of the exons
    return _check_hit_overlap(target_transcript, [hit.left_coords], [hit.right_coords])

def _genome_uniqueness_classify(genome_data, target_transcript, hits):
    on_targets = []
    off_targets = []
    non_targets = []
    for hit in hits:
        if not _is_target_chrom(genome_data, hit.left_coords.name):
            non_targets.append(hit)
        elif _check_genome_hit_overlap(target_transcript, hit):
            on_targets.append(hit)
        else:
            off_targets.append(hit)
    return on_targets, off_targets, non_targets

def _genome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair):
    """guery to find genomic on and off target hits via an alignment method"""
    max_size = 1_000_000  # arbitrary
    hits = uniqueness_query.query_genome(ppair_id, primer3_pair.PRIMER_LEFT_SEQUENCE, primer3_pair.PRIMER_RIGHT_SEQUENCE, max_size)
    return _genome_uniqueness_classify(uniqueness_query.genome_data, target_transcript, hits)

def _check_transcriptome_hit_overlap(target_transcript, hit):
    """does a transcriptome uniqueness hit correspond to the target regions"""
    return _check_hit_overlap(target_transcript,
                              hit.left_features.genome_coords_type(ExonFeature),
                              hit.right_features.genome_coords_type(ExonFeature))

def _transcriptome_uniqueness_classify(genome_data, target_transcript, hits):
    on_targets = []
    off_targets = []
    non_targets = []
    for hit in hits:
        if not _is_target_chrom(genome_data, hit.left_features[0].genome.name):
            non_targets.append(hit)
        elif _check_transcriptome_hit_overlap(target_transcript, hit):
            on_targets.append(hit)
        else:
            off_targets.append(hit)
    return on_targets, off_targets, non_targets

def _transcriptome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair):
    """guery to find transcriptome on and off target hits via an alignment method"""
    max_size = 500_000  # arbitrary
    hits = uniqueness_query.query_transcriptome(ppair_id, primer3_pair.PRIMER_LEFT_SEQUENCE, primer3_pair.PRIMER_RIGHT_SEQUENCE, max_size)
    return _transcriptome_uniqueness_classify(uniqueness_query.genome_data, target_transcript, hits)

def _build_primer_design(target_transcript, target_id, result_num, primer3_pair, uniqueness_query):
    ppair_id = "{}_pp{}".format(target_id, result_num)
    features_5p = _get_exon_left_features(target_transcript, primer3_pair.PRIMER_LEFT)
    features_3p = _get_exon_right_features(target_transcript, primer3_pair.PRIMER_RIGHT)
    _validate_primer_features(features_5p, features_3p)

    if uniqueness_query is None:
        genome_on_targets = genome_off_targets = genome_non_targets = None
        transcriptome_on_targets = transcriptome_off_targets = transcriptome_non_targets = None
    else:
        genome_on_targets, genome_off_targets, genome_non_targets = _genome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair)
        transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets = _transcriptome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair)

    return PrimerDesign(ppair_id, primer3_pair, features_5p, features_3p,
                        genome_on_targets, genome_off_targets, genome_non_targets,
                        transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets)

def _calc_design_status(primer_design) -> DesignStatus:
    if primer_design.transcriptome_off_target_cnt > 0:
        return DesignStatus.NOT_TRANSCRIPTOME_UNIQUE
    elif primer_design.genome_off_target_cnt > 0:
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
    on_cnt = primer_design.genome_on_target_cnt + primer_design.transcriptome_on_target_cnt
    off_cnt = primer_design.genome_off_target_cnt + primer_design.transcriptome_off_target_cnt
    non_cnt = primer_design.genome_non_target_cnt + primer_design.transcriptome_non_target_cnt
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

def design_primers(genome_data, primer_targets, *, uniqueness_query=None, primer3_debug=False):
    """design transcripts """
    target_transcript = primer_targets.transcripts[0]
    primer3_results = primer3_design(target_transcript, debug=primer3_debug)

    return _build_primer_designs(primer_targets, target_transcript, primer3_results, uniqueness_query)

def primer_design_amplicon(primer_design, target_transcript):
    """return amplicon for and RNA and primer"""
    rna = target_transcript.rna
    amplicon_coords = primer_design.amplicon_trans_coords()
    assert amplicon_coords.end <= len(rna)
    return rna[amplicon_coords.start:amplicon_coords.start + amplicon_coords.size]
