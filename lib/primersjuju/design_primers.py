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
    features_5p: Features
    features_3p: Features
    # results from uniqueness query, non_target are ones on sequences liked fixes and alts
    genome_on_targets: Sequence[GenomeHit]
    genome_off_targets: Sequence[GenomeHit]
    genome_non_targets: Sequence[GenomeHit]
    # results from transcriptome query
    transcriptome_on_targets: Sequence[TranscriptomeHit]
    transcriptome_off_targets: Sequence[TranscriptomeHit]
    transcriptome_non_targets: Sequence[TranscriptomeHit]

    def genome_on_target_cnt(self):
        return _len_none(self.genome_on_targets)

    def genome_off_target_cnt(self):
        return _len_none(self.genome_off_targets)

    def genome_non_target_cnt(self):
        return _len_none(self.genome_non_targets)

    def transcriptome_on_target_cnt(self):
        return _len_none(self.transcriptome_on_targets)

    def transcriptome_off_target_cnt(self):
        return _len_none(self.transcriptome_off_targets)

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
    designs: Sequence[PrimerDesign]
    status: DesignStatus

    def dump(self, fh):
        print(f">>> PrimerDesigns {len(self.designs)} <<<", file=fh)
        print("    target_id", self.target_id, file=fh)
        print("    tran_id", self.target_transcript.track_name, self.target_transcript.trans_id, file=fh)
        for design in self.designs:
            design.dump(fh)

def _get_exon_features(target_transcript, primer3_coords):
    feat0 = target_transcript.features_5p[0]
    start = primer3_coords[0] - 1  # one-based
    end = start + primer3_coords[1]  # has length
    trans_coords = Coords(feat0.trans.name, start, end,
                          strand='+', size=feat0.trans.size)
    feats = transcript_range_to_features(target_transcript.features, trans_coords)
    print("primer3_coords", primer3_coords)
    for i, feat in enumerate(feats):
        print("   F", i, str(feat), str(feat.trans.reverse()))
    return feats
    #return transcript_range_to_features(target_transcript.features, trans_coords)

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
        # no information so consider in a possible target
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
    features_5p, features_3p = target_transcript.get_genome_ordered_features()
    return (_check_coords_overlap(features_5p, left_coordses) and
            _check_coords_overlap(features_3p, right_coordses))

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

def _transcriptome_uniqueness_check(genome_data, target_transcript, hits):
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
    return _transcriptome_uniqueness_check(uniqueness_query.genome_data, target_transcript, hits)

def _build_primer_design(target_transcript, target_id, result_num, primer3_pair, uniqueness_query):
    ppair_id = "{}+pp{}".format(target_id, result_num)
    print("\n@@", ppair_id)
    features_5p = _get_exon_features(target_transcript, primer3_pair.PRIMER_LEFT)
    features_3p = _get_exon_features(target_transcript, primer3_pair.PRIMER_RIGHT)
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
    if primer_design.transcriptome_off_target_cnt() > 0:
        return DesignStatus.NOT_TRANSCRIPTOME_UNIQUE
    elif primer_design.genome_off_target_cnt() > 0:
        return DesignStatus.NOT_GENOME_UNIQUE
    else:
        return DesignStatus.GOOD

def get_design_status(primer_design_list) -> DesignStatus:
    "determine status base on primer3 and uniqueness"
    design_status = DesignStatus.NO_PRIMERS  # worst
    for primer_design in primer_design_list:
        design_status = min(design_status, _calc_design_status(primer_design))
    return design_status

def _build_primer_designs(primer_targets, target_transcript, primer3_results, uniqueness_query):
    primer_design_list = [_build_primer_design(target_transcript, primer_targets.target_id, i + 1, pair, uniqueness_query)
                          for i, pair in enumerate(primer3_results.pairs)]
    return PrimerDesigns(primer_targets.target_id,
                         primer_targets, target_transcript, primer3_results,
                         primer_design_list, get_design_status(primer_design_list))

def design_primers(genome_data, primer_targets, *, uniqueness_query=None):
    """design transcripts """
    target_transcript = primer_targets.transcripts[0]
    primer3_results = primer3_design(target_transcript)

    return _build_primer_designs(primer_targets, target_transcript, primer3_results, uniqueness_query)
