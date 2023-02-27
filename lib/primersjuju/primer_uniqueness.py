"""
Primer uniqueness query and results
"""

from typing import Sequence
from dataclasses import dataclass
from pycbio.ncbi.assembly import AssemblyReportNotFound
from .transcript_features import ExonFeature
from .uniqueness_query import GenomeHit, TranscriptomeHit

def _len_none(l):
    "0 if list is None, else length"
    return 0 if l is None else len(l)

@dataclass
class PrimerUniqueness:
    "Results or querying a primer design for genome and transcriptone uniqueness."

    # results from uniqueness query, non_target are ones on sequences liked fixes and alts
    # None if uniqueness check was not done.
    genome_on_targets: Sequence[GenomeHit]
    genome_off_targets: Sequence[GenomeHit]
    genome_non_targets: Sequence[GenomeHit]
    # results from transcriptome query
    transcriptome_on_targets: Sequence[TranscriptomeHit]
    transcriptome_off_targets: Sequence[TranscriptomeHit]
    transcriptome_non_targets: Sequence[TranscriptomeHit]

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

def _is_target_chrom(genome_data, chrom_name):
    """Should a chromosome be consider at target for off-target/on target
    check.  This excludes patches and alts.  For genomes without this information,
    all sequences are conserved a target.
    """
    if genome_data.assembly_info is None:
        # no information so consider this a possible target
        return True
    try:
        chrom_info = genome_data.assembly_info.getByName(chrom_name)
    except AssemblyReportNotFound:
        # can be due to UCSC keep sequences that were removed
        return False
    return ((chrom_info.sequenceRole == "assembled-molecule") and
            (chrom_info.assemblyUnit == "Primary Assembly"))

def _check_gcoords_overlap(features, gcoordss):
    for f in features.iter_type(ExonFeature):
        for gcoords in gcoordss:
            if f.genome.overlaps(gcoords):
                return True
    return False

def _check_hit_overlap(target_transcript, left_gcoordss, right_gcoordss):
    """check overlap based on list of genomic coordinates"""
    features_first, features_last = target_transcript.get_genome_ordered_features()
    return (_check_gcoords_overlap(features_first, left_gcoordss) and
            _check_gcoords_overlap(features_last, right_gcoordss))

def _check_genome_hit_overlap(target_transcript, hit):
    """does a genome uniqueness hit correspond to the target regions"""
    # Genome hits will not span introns, but may partially align one
    # of the exons
    return _check_hit_overlap(target_transcript, [hit.left_gcoords], [hit.right_gcoords])

def _genome_uniqueness_classify(genome_data, target_transcript, hits):
    on_targets = []
    off_targets = []
    non_targets = []
    for hit in hits:
        if not _is_target_chrom(genome_data, hit.left_gcoords.name):
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

def primer_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair):
    "Run genome and transcriptome queries and collate results"
    genome_on_targets, genome_off_targets, genome_non_targets = _genome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair)
    transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets = _transcriptome_uniqueness_query(uniqueness_query, target_transcript, ppair_id, primer3_pair)

    return PrimerUniqueness(genome_on_targets, genome_off_targets, genome_non_targets,
                            transcriptome_on_targets, transcriptome_off_targets, transcriptome_non_targets)

def primer_uniqueness_none():
    "When uniqueness_query is not available"
    return PrimerUniqueness(None, None, None, None, None, None)
