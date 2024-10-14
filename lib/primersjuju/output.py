# common output functions

import os.path as osp
import urllib.parse
from pycbio.sys import fileOps
from pycbio.sys.svgcolors import SvgColors
from pycbio.hgdata.bed import Bed
from primersjuju.transcript_features import features_to_genomic_coords_list, features_sort_genome
from primersjuju.primer_targets import ExonFeature
from primersjuju.design_primers import primer_design_amplicon, primer_design_amplicon_features, primer_design_amplicon_coords
from primersjuju.primer3_interface import primer3_dump_args, primer3_annotate_amplicon

# Document these in README

# target tracks
TARGET_SPEC_COLOR = SvgColors.blue
TARGET_FEAT_COLOR = SvgColors.green

# amplicon tracks
AMP_TARGET_COLOR = SvgColors.orange
AMP_PRIMER_COLOR = SvgColors.indigo

# primer tracks
PRIMERS_UNSTABLE_COLOR = SvgColors.fuchsia
PRIMERS_ON_COLOR = SvgColors.green
PRIMERS_ON_OFF_COLOR = SvgColors.orange
PRIMERS_OFF_COLOR = SvgColors.red
PRIMERS_NONE_COLOR = SvgColors.darkorange
PRIMERS_NON_COLOR = SvgColors.purple

# uniqueness track
UNIQ_ON_COLOR = SvgColors.green
UNIQ_OFF_COLOR = SvgColors.red
UNIQ_NON_COLOR = SvgColors.purple


#  Delta G: secondary structures of any self-dimers, hairpins, and
#  heterodimers that should be weaker (more positive) thanÂ -9.0 kcal/mole.
STABILITY_THRSEHOLD = -9.0

def _coords_list_to_range(coords_list):
    return coords_list[0].adjrange(None, end=coords_list[-1].end)

def _gcoords_to_bed(name, color, gcoords_list, *, strand=None, extra_cols=None, thick_gcoords=None):
    gcoords_list = sorted(gcoords_list, key=lambda c: (c.name, c.start, c.end))
    first = gcoords_list[0]
    last = gcoords_list[-1]
    assert first < last
    if strand is None:
        strand = first.strand
    assert strand is not None
    if thick_gcoords is None:
        thick_gcoords = first.adjrange(first.start, last.end)

    bed = Bed(first.name, first.start, last.end, name, strand=strand,
              thickStart=thick_gcoords.start, thickEnd=thick_gcoords.end,
              itemRgb=color.toRgb8Str(), numStdCols=12,
              extraCols=extra_cols)
    for gcoords in gcoords_list:
        bed.addBlock(gcoords.start, gcoords.end)
    return bed

def make_transcript_bed(primer_targets, trans, name, start, end, color):
    """Generate a transcript BED with an amplication regions annotation and a custom name."""
    thick_gcoords = trans.bounds.genome.adjrange(start, end)
    return _gcoords_to_bed(name, color,
                           features_to_genomic_coords_list(trans.features, ExonFeature),
                           strand=primer_targets.strand, thick_gcoords=thick_gcoords)

def build_target_transcript_bed(primer_targets, trans):
    "transcript, with maximum possbile amplicon as thick"
    features_first, features_last = trans.get_genome_ordered_features()
    return make_transcript_bed(primer_targets, trans, trans.trans_id.name,
                               features_first[0].genome.start, features_last[-1].genome.end,
                               TARGET_FEAT_COLOR)

def build_target_beds(primer_targets):
    # specified target regions
    target_bed = [_gcoords_to_bed(primer_targets.target_id, TARGET_SPEC_COLOR,
                                  [primer_targets.region_5p, primer_targets.region_3p],
                                  strand=primer_targets.strand,
                                  thick_gcoords=primer_targets.region_5p)]
    trans_beds = [build_target_transcript_bed(primer_targets, trans)
                  for trans in primer_targets.transcripts]

    return target_bed + trans_beds

def build_target_amplicon_bed(primer_targets, trans, primer_design):
    primer_feats = primer_design.features_5p + primer_design.features_3p
    features_sort_genome(primer_feats)
    thick_gcoords = trans.bounds.genome.adjrange(primer_feats[0].genome.start, primer_feats[-1].genome.end)
    name = trans.trans_id.name + "^" + primer_design.ppair_id
    return make_transcript_bed(primer_targets, trans, name,
                               thick_gcoords.start, thick_gcoords.end, AMP_TARGET_COLOR)

def build_amplicon_beds(primer_targets, primer_designs):
    if len(primer_designs.designs) == 0:
        return []
    primer_design = primer_designs.designs[0]
    assert primer_design.priority == 1

    primer_bed = [_primer_to_bed(primer_designs, primer_design, color=AMP_PRIMER_COLOR, add_extra=False)]
    trans_beds = [build_target_amplicon_bed(primer_targets, trans, primer_design)
                  for trans in primer_targets.transcripts]
    return primer_bed + trans_beds

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

def _get_extra_cols(primer_designs, primer_design):
    """get extra BED columns from primer design"""
    uniqueness = primer_design.uniqueness
    extra_cols = [primer_design.priority,
                  primer_design.amplicon_length,
                  _count_amplicon_exons(primer_design, primer_designs.target_transcript),
                  uniqueness.transcriptome_on_target_cnt,
                  uniqueness.transcriptome_off_target_cnt,
                  uniqueness.transcriptome_non_target_cnt,
                  uniqueness.genome_on_target_cnt,
                  uniqueness.genome_off_target_cnt,
                  uniqueness.genome_non_target_cnt]
    # directly from primer3
    for col_name in _primer_bed_columns:
        col = getattr(primer_design.primer3_pair, col_name)
        if isinstance(col, float):
            col = "{:.2f}".format(col)
        extra_cols.append(col)
    return extra_cols

def _primer_color(primer_design):
    if ((primer_design.primer3_pair.PRIMER_LEFT_END_STABILITY <= STABILITY_THRSEHOLD) or
        (primer_design.primer3_pair.PRIMER_RIGHT_END_STABILITY <= STABILITY_THRSEHOLD)):
        return PRIMERS_UNSTABLE_COLOR
    uniqueness = primer_design.uniqueness
    if (((uniqueness.transcriptome_on_target_cnt + uniqueness.genome_on_target_cnt) > 0) and
        ((uniqueness.transcriptome_off_target_cnt + uniqueness.genome_off_target_cnt) > 0)):
        return PRIMERS_ON_OFF_COLOR
    if ((uniqueness.transcriptome_on_target_cnt + uniqueness.genome_on_target_cnt) > 0):
        return PRIMERS_ON_COLOR
    if ((uniqueness.transcriptome_off_target_cnt + uniqueness.genome_off_target_cnt) > 0):
        return PRIMERS_OFF_COLOR
    if ((uniqueness.transcriptome_non_target_cnt + uniqueness.genome_non_target_cnt) > 0):
        return PRIMERS_NON_COLOR
    return PRIMERS_NONE_COLOR

def _primer_to_bed(primer_designs, primer_design, *, color=None, add_extra=True):
    if color is None:
        color = _primer_color(primer_design)
    gcoords_5p_list = features_to_genomic_coords_list(primer_design.features_5p, ExonFeature)
    gcoords_3p_list = features_to_genomic_coords_list(primer_design.features_3p, ExonFeature)
    extra_cols = _get_extra_cols(primer_designs, primer_design) if add_extra else None
    return _gcoords_to_bed(primer_design.ppair_id, color,
                           gcoords_5p_list + gcoords_3p_list,
                           strand=primer_designs.primer_targets.strand,
                           thick_gcoords=_coords_list_to_range(gcoords_5p_list),
                           extra_cols=extra_cols)

def build_primer_beds(primer_designs):
    return [_primer_to_bed(primer_designs, pd) for pd in primer_designs.designs]

def _genome_hit_to_bed(hit, name, color):
    return _gcoords_to_bed(name, color, (hit.left_gcoords, hit.right_gcoords),
                           thick_gcoords=hit.left_gcoords)

def _genome_hits_to_bed(hits, name, color):
    if hits is None:
        return []  # no data
    return [_genome_hit_to_bed(hit, name, color) for hit in hits]

def _build_genome_uniqueness_hits_beds(primer_design):
    beds = []
    beds += _genome_hits_to_bed(primer_design.uniqueness.genome_on_targets, "on:" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _genome_hits_to_bed(primer_design.uniqueness.genome_off_targets, "off:" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _genome_hits_to_bed(primer_design.uniqueness.genome_non_targets, "non:" + primer_design.ppair_id, UNIQ_NON_COLOR)
    return beds

def build_genome_uniqueness_hits_beds(primer_designs):
    beds = []
    for primer_design in primer_designs.designs:
        beds.extend(_build_genome_uniqueness_hits_beds(primer_design))
    return beds

def _transcriptome_hit_to_bed(hit, name_pre, color):
    name = name_pre + '|' + hit.trans_id
    if hit.gene_name is not None:
        name += '|' + hit.gene_name
    left_gcoords = features_to_genomic_coords_list(hit.left_features, ExonFeature)
    right_gcoords = features_to_genomic_coords_list(hit.right_features, ExonFeature)
    return _gcoords_to_bed(name, color, left_gcoords + right_gcoords,
                           thick_gcoords=_coords_list_to_range(left_gcoords))

def _transcriptome_hits_to_bed(hits, name_pre, color):
    if hits is None:
        return []  # no data
    return [_transcriptome_hit_to_bed(hit, name_pre, color) for hit in hits]

def _build_transcriptome_uniqueness_hits_beds(primer_design):
    beds = []
    beds += _transcriptome_hits_to_bed(primer_design.uniqueness.transcriptome_on_targets, "on:" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.uniqueness.transcriptome_off_targets, "off:" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.uniqueness.transcriptome_non_targets, "non:" + primer_design.ppair_id, UNIQ_NON_COLOR)
    return beds

def build_transcriptome_uniqueness_hits_beds(primer_designs):
    beds = []
    for primer_design in primer_designs.designs:
        beds.extend(_build_transcriptome_uniqueness_hits_beds(primer_design))
    return beds

def _write_beds(beds, bed_file):
    beds = sorted(beds, key=Bed.genome_sort_key)
    with open(bed_file, "w") as fh:
        for bed in beds:
            bed.write(fh)

def _get_out_path(outdir, target_id, suffix):
    return osp.join(outdir, target_id + "." + suffix)

def output_target_design_file(outdir, target_id):
    """get path to target design TSV file for an experiment.  The existence of this file indicates design
    for this target was complete.
    """
    return _get_out_path(outdir, target_id, "designs.tsv")

def _make_excel_link(url, position):
    return f'=HYPERLINK("{url}", "{str(position)}")'

def _make_browser_link(genome_name, position, hub_urls=None):
    browser_url = "https://genome.ucsc.edu/cgi-bin/hgTracks"

    # browser doesn't allow entire string to be quotes, only arguments
    cgi_args = ["position=" + urllib.parse.quote(str(position)),
                "db=" + urllib.parse.quote(genome_name),
                "genome=" + urllib.parse.quote(genome_name)]
    if hub_urls is not None:
        for hub_url in hub_urls:
            cgi_args.append("hubUrl=" + urllib.parse.quote(hub_url))
    url = browser_url + "?" + "&".join(cgi_args)
    return _make_excel_link(url, position)

def _make_design_browser_link(primer_designs, hub_urls):
    return _make_browser_link(primer_designs.primer_targets.genome_name,
                              primer_designs.target_transcript.bounds.genome, hub_urls)

def _make_uniqeness_hits_browser_gcoords(hits):
    """this makes a list of coordinates, there isn't a way to add multiple
    hyperlinks to a cell via a TSV"""
    # drop duplicates
    if hits is None:
        return ""
    else:
        gcoords_list = sorted(set([h.get_genome_range() for h in hits]))
        return ", ".join([str(c) for c in gcoords_list])

def _count_amplicon_exons(primer_design, target_transcript):
    amp_features = primer_design_amplicon_features(primer_design, target_transcript)
    return len([exon for exon in amp_features.iter_type(ExonFeature)])

def output_target_debug(config, outdir, primer_targets, primer_designs):
    target_transcript = primer_targets.transcripts[0]
    fileOps.ensureDir(outdir)
    with open(_get_out_path(outdir, primer_targets.target_id, "debug.txt"), 'w') as fh:
        primer_targets.dump(fh)
        print("annotated_amplicon:", primer3_annotate_amplicon(config.primer3, primer_targets.transcripts[0]), file=fh)
        primer3_dump_args(fh, config.primer3, target_transcript)
        print(file=fh)
        primer_designs.primer3_results.dump(fh)
        print(file=fh)
        primer_designs.dump(fh)
        print(file=fh)

def output_target_beds(outdir, primer_targets, primer_designs):
    fileOps.ensureDir(outdir)
    _write_beds(build_target_beds(primer_targets),
                _get_out_path(outdir, primer_targets.target_id, "target.bed"))
    _write_beds(build_amplicon_beds(primer_targets, primer_designs),
                _get_out_path(outdir, primer_targets.target_id, "amplicon.bed"))
    _write_beds(build_primer_beds(primer_designs),
                _get_out_path(outdir, primer_targets.target_id, "primers.bed"))
    if primer_designs.uniqueness_checked:
        _write_beds(build_genome_uniqueness_hits_beds(primer_designs),
                    _get_out_path(outdir, primer_targets.target_id, "genome-uniqueness.bed"))
        _write_beds(build_transcriptome_uniqueness_hits_beds(primer_designs),
                    _get_out_path(outdir, primer_targets.target_id, "transcriptome-uniqueness.bed"))


_design_tsv_header = ("target_id", "transcript_id", "design_status", "position", "browser",
                      "primer_id", "left_primer", "right_primer", "pri",
                      "amplicon_len", "amplicon_exons", "left_delta_G", "right_delta_G",
                      "on_target_trans", "off_target_trans",
                      "on_target_genome", "off_target_genome",
                      "amplicon")

def _write_primer_pair_design_trans(fh, primer_designs, primer_design, trans, first, hub_urls):
    "write one design, if primer_design is None, it means there were no primers found"
    row = [primer_designs.primer_targets.target_id, trans.trans_id.name,
           primer_designs.status]
    if first:
        row.extend([primer_designs.target_transcript.bounds.genome,
                    _make_design_browser_link(primer_designs, hub_urls)])
    else:
        row.extend(2 * [''])
    if primer_design is None:
        row += 13 * ['']
    else:
        amp_seq = primer_design_amplicon(primer_design, trans)
        row += [primer_design.ppair_id,
                primer_design.primer3_pair.PRIMER_LEFT_SEQUENCE,
                primer_design.primer3_pair.PRIMER_RIGHT_SEQUENCE,
                primer_design.priority,
                len(amp_seq),
                _count_amplicon_exons(primer_design, trans),
                primer_design.primer3_pair.PRIMER_LEFT_END_STABILITY,
                primer_design.primer3_pair.PRIMER_RIGHT_END_STABILITY,
                _make_uniqeness_hits_browser_gcoords(primer_design.uniqueness.transcriptome_on_targets),
                _make_uniqeness_hits_browser_gcoords(primer_design.uniqueness.transcriptome_off_targets),
                _make_uniqeness_hits_browser_gcoords(primer_design.uniqueness.genome_on_targets),
                _make_uniqeness_hits_browser_gcoords(primer_design.uniqueness.genome_off_targets),
                amp_seq]
    fileOps.prRow(fh, row)

def _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls):
    for trans in primer_designs.primer_targets.transcripts:
        _write_primer_pair_design_trans(fh, primer_designs, primer_design, trans, first, hub_urls)

def _write_primer_designs(fh, primer_designs, hub_urls):
    fileOps.prRow(fh, _design_tsv_header)
    first = True
    if len(primer_designs.designs) == 0:
        _write_primer_pair_design(fh, primer_designs, None, True, hub_urls)
    else:
        for primer_design in sorted(primer_designs.designs, key=lambda pd: pd.priority):
            _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls)
            first = False

def output_primer_designs(outdir, primer_targets, primer_designs, hub_urls):
    fileOps.ensureDir(outdir)
    with fileOps.AtomicFileCreate(_get_out_path(outdir, primer_targets.target_id, "designs.tsv")) as tmp_tsv:
        with open(tmp_tsv, "w") as fh:
            _write_primer_designs(fh, primer_designs, hub_urls)

_isoform_tsv_header = ("target_id", "primer_id", "pri", "track", "transcript_id",
                       "amplicon_coords", "amplicon_len", "amplicon_exons", "amplicon")

def _write_primer_pair_isoform(fh, primer_targets, target_transcript, primer_design):
    "write isoform information for one primer design"

    # map primer regions to this transcript
    amp_tcoords = primer_design_amplicon_coords(primer_design, target_transcript)

    fileOps.prRowv(fh, primer_targets.target_id, primer_design.ppair_id, primer_design.priority,
                   target_transcript.trans_id.track, target_transcript.trans_id.name,
                   amp_tcoords, len(amp_tcoords),
                   _count_amplicon_exons(primer_design, target_transcript),
                   primer_design_amplicon(primer_design, target_transcript))

def _write_primers_isoforms(fh, primer_targets, primer_designs):
    fileOps.prRow(fh, _isoform_tsv_header)
    for primer_design in primer_designs.designs:
        for target_transcript in primer_targets.transcripts:
            _write_primer_pair_isoform(fh, primer_targets, target_transcript, primer_design)

def output_primers_isoforms(outdir, primer_targets, primer_designs):
    fileOps.ensureDir(outdir)
    with fileOps.AtomicFileCreate(_get_out_path(outdir, primer_targets.target_id, "isoforms.tsv")) as tmp_tsv:
        with open(tmp_tsv, "w") as fh:
            _write_primers_isoforms(fh, primer_targets, primer_designs)

def output_target_designs(config, outdir, primer_targets, primer_designs, hub_urls=None):
    """output primer TSV, BEDs, and debug information for one target.  The
    $target_id.design.tsv file is created atomically, so it can be used as a
    marker that this design is complete.  If hub_urls is specified, the
    design.tsv file will contain links to the target region in the UCSC
    browser.
    """

    output_target_debug(config, outdir, primer_targets, primer_designs)
    output_target_beds(outdir, primer_targets, primer_designs)
    output_primers_isoforms(outdir, primer_targets, primer_designs)
    output_primer_designs(outdir, primer_targets, primer_designs, hub_urls)
