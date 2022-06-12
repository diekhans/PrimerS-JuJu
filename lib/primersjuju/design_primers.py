"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from pycbio.hgdata.bed import Bed
from . import PrimersJuJuError
from .primer3_interface import Primer3Results, Primer3Pair, primer3_design
from .target_transcripts import TargetTranscripts, TargetTranscript, ExonFeature
from .transcript_features import Features, transcript_range_to_features, features_to_genomic_coords
from .uniqueness_query import GenomeHit, TranscriptomeHit

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

    def dump(self, fh):
        print(">>> PrimerDesign <<<", file=fh)
        print("    ppair_id", self.ppair_id, file=fh)
        print("    features_5p", self.ppair_id, file=fh)
        print("    features_3p", self.ppair_id, file=fh)
        print("    genome_on_targets", self.genome_on_targets, file=fh)
        print("    genome_off_targets", self.genome_off_targets, file=fh)
        print("    genome_non_targets", self.genome_non_targets, file=fh)
        print("    transcriptome_on_targets", self.transcriptome_on_targets, file=fh)
        print("    transcriptome_off_targets", self.transcriptome_off_targets, file=fh)
        print("    transcriptome_non_targets", self.transcriptome_non_targets, file=fh)


@dataclass
class PrimerDesigns:
    """information collected on all primers for a target"""
    target_id: str
    target_transcripts: TargetTranscripts
    target_transcript: TargetTranscript
    primer3_results: Primer3Results
    designs: Sequence[PrimerDesign]

    def dump(self, fh):
        print(f">>> PrimerDesigns {len(self.designs)} <<<", file=fh)
        print("    target_id", self.target_id, file=fh)
        print("    tran_id", self.target_transcript.track_name, self.target_transcript.trans_id, file=fh)
        for design in self.designs:
            design.dump(fh)

def _get_exon_features(target_transcript, primer3_coords):
    some_feat = target_transcript.features_5p[0]
    start = primer3_coords[0] - 1  # one-based
    end = start + primer3_coords[1]  # has length
    trans_coords = Coords(some_feat.trans.name, start, end,
                          strand=some_feat.trans.strand, size=some_feat.trans.size)
    return transcript_range_to_features(target_transcript.features, trans_coords)

def _validate_primer_features(features_5p, features_3p):
    for feature_5p in features_5p:
        for feature_3p in features_3p:
            if feature_5p.trans.overlaps(feature_3p.trans):
                raise PrimersJuJuError(f"primer3 pairs overlap in transcript space {features_5p} and {features_3p}")
            if feature_5p.genome.overlaps(feature_3p.genome):
                raise PrimersJuJuError(f"primer3 pairs overlap in genome space {features_5p} and {features_3p}")

def _is_target_chrom(uniqueness_query, chrom_name):
    """Should a chromosome be consider at target for off-target/on target
    check.  This excludes patches and alts.  For genomes without this information,
    all sequences are conserved a target.
    """
    if uniqueness_query.genome_data.assembly_info is None:
        # no information so consider in a possible target
        return True
    chrom_info = uniqueness_query.genome_data.assembly_info.getByName(chrom_name)
    return ((chrom_info.sequenceRole == "assembled-molecule") and
            (chrom_info.assemblyUnit == "Primary Assembly"))

def _check_genome_hit_overlap(target_transcript, hit):
    """does a genome uniqueness hit correspond to the target regions"""
    # genome hits can't span splice sites
    if ((len(target_transcript.features_5p) > 1) or
        (len(target_transcript.features_3p) > 1)):
        return False
    return (target_transcript.features_5p[0].genome.contains(hit.left_coords) and
            target_transcript.features_3p[0].genome.contains(hit.right_coords))

def _is_genome_off_target(uniqueness_query, target_transcript, hit):
    return (_is_target_chrom(uniqueness_query, hit.left_coords.name) and
            (not _check_genome_hit_overlap(target_transcript, hit)))

def _genome_uniqueness_classify(uniqueness_query, target_transcript, ppair_id, primer3_pair):
    """guery to find genomic on and off target hits via an alignment method"""
    max_size = 1_000_000  # arbitrary
    hits = uniqueness_query.query_genome(ppair_id, primer3_pair.PRIMER_LEFT_SEQUENCE, primer3_pair.PRIMER_RIGHT_SEQUENCE, max_size)
    on_targets = []
    off_targets = []
    non_targets = []
    for hit in hits:
        if not _is_target_chrom(uniqueness_query, hit.left_coords.name):
            non_targets.append(hit)
        elif _check_genome_hit_overlap(target_transcript, hit):
            on_targets.append(hit)
        else:
            off_targets.append(hit)
    return on_targets, off_targets, non_targets

def _check_transcriptome_hit_overlap(target_transcript, hit):
    """does a transcriptome uniqueness hit correspond to the target regions"""
    def _check_one(trans_features, hit_features):
        # check both hit and regions must be junction spanning or not junction spanning
        if len(trans_features) != len(hit_features):
            return False
        check_idxs = (0,) if len(trans_features) == 1 else (0, 2)  # exon or exon-intron-exon
        for i in check_idxs:
            if not trans_features[i].genome.contains(hit_features[i].genome):
                return False
        return True

    return (_check_one(target_transcript.features_5p, hit.left_features) and
            _check_one(target_transcript.features_3p, hit.right_features))

def _transcriptome_uniqueness_classify(uniqueness_query, target_transcript, ppair_id, primer3_pair):
    """guery to find transcriptome on and off target hits via an alignment method"""
    max_size = 500_000  # arbitrary
    hits = uniqueness_query.query_transcriptome(ppair_id, primer3_pair.PRIMER_LEFT_SEQUENCE, primer3_pair.PRIMER_RIGHT_SEQUENCE, max_size)
    on_targets = []
    off_targets = []
    non_targets = []
    for hit in hits:
        if not _is_target_chrom(uniqueness_query, hit.left_features[0].genome.name):
            non_targets.append(hit)
        elif _check_transcriptome_hit_overlap(target_transcript, hit):
            on_targets.append(hit)
        else:
            off_targets.append(hit)
    return on_targets, off_targets, non_targets

def _build_primer_design(target_transcript, target_id, result_num, primer3_pair, uniqueness_query):
    ppair_id = "{}+pp{}".format(target_id, result_num)
    features_5p = _get_exon_features(target_transcript, primer3_pair.PRIMER_LEFT)
    features_3p = _get_exon_features(target_transcript, primer3_pair.PRIMER_RIGHT)
    _validate_primer_features(features_5p, features_3p)

    if uniqueness_query is None:
        genome_on_targets = genome_off_targets = genome_non_targets = None
        transcriptome_on_targets = transcriptome_off_targets = transcriptome_non_targets = None
    else:
        genome_on_targets, genome_off_targets, genome_non_targets = _genome_uniqueness_classify(uniqueness_query, target_transcript, ppair_id, primer3_pair)
        transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets = _transcriptome_uniqueness_classify(uniqueness_query, target_transcript, ppair_id, primer3_pair)

    return PrimerDesign(ppair_id, primer3_pair, features_5p, features_3p,
                        genome_on_targets, genome_off_targets, genome_non_targets,
                        transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets)

def _build_primer_designs(target_transcripts, target_transcript, primer3_results, uniqueness_query):
    return PrimerDesigns(target_transcripts.target_id,
                         target_transcripts, target_transcript, primer3_results,
                         [_build_primer_design(target_transcript, target_transcripts.target_id, i + 1, pair, uniqueness_query)
                          for i, pair in enumerate(primer3_results.pairs)])

def design_primers(genome_data, target_transcripts, *, uniqueness_query=None):
    """design transcripts """
    target_transcript = target_transcripts.transcripts[0]
    primer3_results = primer3_design(target_transcript)

    return _build_primer_designs(target_transcripts, target_transcript, primer3_results, uniqueness_query)

_primer_bed_columns = (
    'PRIMER_PAIR_PRODUCT_SIZE',    # 1133
    'PRIMER_PAIR_COMPL_ANY_TH',    # 12.498152433966595
    'PRIMER_PAIR_COMPL_END_TH',    # 5.033689592343535
    'PRIMER_PAIR_PENALTY',         # 0.42111006752037383
    'PRIMER_LEFT',                 # (27, 20)
    'PRIMER_LEFT_SEQUENCE',        # 'GGTTCTTCTGCGCTACTGCT'
    'PRIMER_LEFT_END_STABILITY',   # 4.24
    'PRIMER_LEFT_GC_PERCENT',      # 55.0
    'PRIMER_LEFT_HAIRPIN_TH',      # 36.97289132676252
    'PRIMER_LEFT_PENALTY',         # 0.3904716385882807
    'PRIMER_LEFT_SELF_ANY_TH',     # 9.732052967213463
    'PRIMER_LEFT_SELF_END_TH',     # 0.0
    'PRIMER_LEFT_TM',              # 60.39047163858828
    'PRIMER_RIGHT',                # (1159, 20)
    'PRIMER_RIGHT_SEQUENCE',       # 'CAAAAACCCACGCAGACAGG'
    'PRIMER_RIGHT_END_STABILITY',  # 4.0
    'PRIMER_RIGHT_GC_PERCENT',     # 55.0
    'PRIMER_RIGHT_HAIRPIN_TH',     # 0.0
    'PRIMER_RIGHT_PENALTY',        # 0.03063842893209312
    'PRIMER_RIGHT_SELF_ANY_TH',    # 0.0
    'PRIMER_RIGHT_SELF_END_TH',    # 0.0
    'PRIMER_RIGHT_TM',             # 59.96936157106791
)

def _get_extra_cols(primer_design):
    """get extra BED columns from primer design"""
    extra_cols = []
    for col_name in _primer_bed_columns:
        col = getattr(primer_design.primer3_pair, col_name)
        if isinstance(col, float):
            col = "{:.2f}".format(col)
        extra_cols.append(col)
    return extra_cols

def _coords_to_bed(name, color, coords_list, extra_cols=None):
    coords_list = sorted(coords_list, key=lambda c: (c.name, c.start, c.end))
    first = coords_list[0]
    last = coords_list[-1]

    bed = Bed(first.name, first.start, last.end, name,
              strand=first.strand, thickStart=first.start, thickEnd=last.end,
              itemRgb=color.toRgb8Str(), numStdCols=12, blocks=[],
              extraCols=extra_cols)
    for coords in coords_list:
        bed.addBlock(coords.start, coords.end)
    return bed

def _primer_to_bed(primer_design, color):
    coords_list = (features_to_genomic_coords(primer_design.features_5p, ExonFeature) +
                   features_to_genomic_coords(primer_design.features_3p, ExonFeature))
    return _coords_to_bed(primer_design.ppair_id, color, coords_list,
                          _get_extra_cols(primer_design))

def build_primer_beds(primer_designs, color):
    return sorted([_primer_to_bed(pd, color) for pd in primer_designs.designs],
                  key=Bed.genome_sort_key)

def _genome_hit_to_bed(hit, name, color):
    return _coords_to_bed(name, color, (hit.left_coords, hit.right_coords))

def _genome_hits_to_bed(hits, name, color):
    if hits is None:
        return []  # no data
    return [_genome_hit_to_bed(hit, name, color) for hit in hits]

def _transcriptome_hit_to_bed(hit, name_pre, color):
    return _coords_to_bed(name_pre + '|' + hit.trans_id + '|' + hit.gene_name, color,
                          features_to_genomic_coords(hit.left_features, ExonFeature) +
                          features_to_genomic_coords(hit.right_features, ExonFeature))

def _transcriptome_hits_to_bed(hits, name_pre, color):
    if hits is None:
        return []  # no data
    return [_transcriptome_hit_to_bed(hit, name_pre, color) for hit in hits]

def _build_design_uniqueness_hits_beds(primer_design, on_color, off_color, non_color):
    beds = []
    beds += _genome_hits_to_bed(primer_design.genome_on_targets, "on=" + primer_design.ppair_id, on_color)
    beds += _genome_hits_to_bed(primer_design.genome_off_targets, "off=" + primer_design.ppair_id, off_color)
    beds += _genome_hits_to_bed(primer_design.genome_off_targets, "non=" + primer_design.ppair_id, non_color)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_on_targets, "on=" + primer_design.ppair_id, on_color)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_off_targets, "off=" + primer_design.ppair_id, off_color)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_non_targets, "non=" + primer_design.ppair_id, non_color)
    return beds

def build_uniqueness_hits_beds(primer_designs, on_color, off_color, non_color):
    beds = []
    for primer_design in primer_designs.designs:
        beds += _build_design_uniqueness_hits_beds(primer_design, on_color, off_color, non_color)
    return sorted(beds, key=Bed.genome_sort_key)
