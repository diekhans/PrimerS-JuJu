# common output functions

import os.path as osp
import urllib.parse
from pycbio.sys import fileOps
from pycbio.sys.svgcolors import SvgColors
from pycbio.hgdata.bed import Bed
from primersjuju.transcript_features import features_to_genomic_coords
from primersjuju.primer_targets import ExonFeature
from primersjuju.design_primers import primer_design_amplicon
from primersjuju.primer3_interface import primer3_dump_args, primer3_annotate_amplicon

# Document these in README

# target track
TARGET_SPEC_COLOR = SvgColors.blue
TARGET_FEAT_COLOR = SvgColors.green

# primer tracks
PRIMERS_UNSTABLE_COLOR = SvgColors.fuchsia
PRIMERS_ON_COLOR = SvgColors.green
PRIMERS_ON_OFF_COLOR = SvgColors.orange
PRIMERS_OFF_COLOR = SvgColors.red
PRIMERS_NONE_COLOR = SvgColors.darkorange
PRIMERS_NON_COLOR = SvgColors.purple

# uniquness track
UNIQ_ON_COLOR = SvgColors.green
UNIQ_OFF_COLOR = SvgColors.red
UNIQ_NON_COLOR = SvgColors.purple


#  Delta G: secondary structures of any self-dimers, hairpins, and
#  heterodimers that should be weaker (more positive) thanÂ -9.0 kcal/mole.
STABILITY_THRSEHOLD = -9.0

def _coords_list_to_range(coords_list):
    return coords_list[0].adjrange(None, end=coords_list[-1].end)

def _coords_to_bed(name, color, coords_list, *, strand=None, extra_cols=None, thick_coords=None):
    coords_list = sorted(coords_list, key=lambda c: (c.name, c.start, c.end))
    first = coords_list[0]
    last = coords_list[-1]
    assert first < last
    if strand is None:
        strand = first.strand
    assert strand is not None
    if thick_coords is None:
        thick_coords = first.adjrange(first.start, last.end)

    bed = Bed(first.name, first.start, last.end, name, strand=strand,
              thickStart=thick_coords.start, thickEnd=thick_coords.end,
              itemRgb=color.toRgb8Str(), numStdCols=12,
              extraCols=extra_cols)
    for coords in coords_list:
        bed.addBlock(coords.start, coords.end)
    return bed

def build_target_beds(primer_targets):
    # specified target regions
    target_beds = [_coords_to_bed(primer_targets.target_id, TARGET_SPEC_COLOR,
                                  [primer_targets.region_5p, primer_targets.region_3p],
                                  strand=primer_targets.strand,
                                  thick_coords=primer_targets.region_5p)]

    # transcript, with amplicon as thick
    trans0 = primer_targets.transcripts[0]

    features_first, features_last = trans0.get_genome_ordered_features()
    thick_coords = trans0.bounds.genome.adjrange(features_first[0].genome.start,
                                                 features_last[-1].genome.end)
    feat_beds = [_coords_to_bed(trans0.trans_id, TARGET_FEAT_COLOR,
                                features_to_genomic_coords(trans0.features, ExonFeature),
                                strand=primer_targets.strand, thick_coords=thick_coords)]
    return target_beds + feat_beds

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
    extra_cols = [primer_design.priority,
                  primer_design.amplicon_length,
                  _count_amplicon_exons(primer_design, primer_designs.target_transcript),
                  primer_design.transcriptome_on_target_cnt,
                  primer_design.transcriptome_off_target_cnt,
                  primer_design.transcriptome_non_target_cnt,
                  primer_design.genome_on_target_cnt,
                  primer_design.genome_off_target_cnt,
                  primer_design.genome_non_target_cnt]
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
    elif (((primer_design.transcriptome_on_target_cnt + primer_design.genome_on_target_cnt) > 0) and
          ((primer_design.transcriptome_off_target_cnt + primer_design.genome_off_target_cnt) > 0)):
        return PRIMERS_ON_OFF_COLOR
    elif ((primer_design.transcriptome_on_target_cnt + primer_design.genome_on_target_cnt) > 0):
        return PRIMERS_ON_COLOR
    elif ((primer_design.transcriptome_off_target_cnt + primer_design.genome_off_target_cnt) > 0):
        return PRIMERS_OFF_COLOR
    elif ((primer_design.transcriptome_non_target_cnt + primer_design.genome_non_target_cnt) > 0):
        return PRIMERS_NON_COLOR
    else:
        return PRIMERS_NONE_COLOR

def _primer_to_bed(primer_designs, primer_design):
    coords_5p_list = features_to_genomic_coords(primer_design.features_5p, ExonFeature)
    coords_3p_list = features_to_genomic_coords(primer_design.features_3p, ExonFeature)
    return _coords_to_bed(primer_design.ppair_id, _primer_color(primer_design),
                          coords_5p_list + coords_3p_list,
                          strand=primer_designs.primer_targets.strand,
                          thick_coords=_coords_list_to_range(coords_5p_list),
                          extra_cols=_get_extra_cols(primer_designs, primer_design))

def build_primer_beds(primer_designs):
    return [_primer_to_bed(primer_designs, pd) for pd in primer_designs.designs]

def _genome_hit_to_bed(hit, name, color):
    return _coords_to_bed(name, color, (hit.left_coords, hit.right_coords),
                          thick_coords=hit.left_coords)

def _genome_hits_to_bed(hits, name, color):
    if hits is None:
        return []  # no data
    return [_genome_hit_to_bed(hit, name, color) for hit in hits]

def _build_genome_uniqueness_hits_beds(primer_design):
    beds = []
    beds += _genome_hits_to_bed(primer_design.genome_on_targets, "on:" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _genome_hits_to_bed(primer_design.genome_off_targets, "off:" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _genome_hits_to_bed(primer_design.genome_non_targets, "non:" + primer_design.ppair_id, UNIQ_NON_COLOR)
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
    left_features = features_to_genomic_coords(hit.left_features, ExonFeature)
    right_features = features_to_genomic_coords(hit.right_features, ExonFeature)
    return _coords_to_bed(name, color, left_features + right_features,
                          thick_coords=_coords_list_to_range(left_features))

def _transcriptome_hits_to_bed(hits, name_pre, color):
    if hits is None:
        return []  # no data
    return [_transcriptome_hit_to_bed(hit, name_pre, color) for hit in hits]

def _build_transcriptome_uniqueness_hits_beds(primer_design):
    beds = []
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_on_targets, "on:" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_off_targets, "off:" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_non_targets, "non:" + primer_design.ppair_id, UNIQ_NON_COLOR)
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

def _make_uniqeness_hits_browser_coords(hits):
    """this makes a list of coordinates, there isn't a way to add multiple
    hyperlinks to a cell via a TSV"""
    # drop duplicates
    if hits is None:
        return ""
    else:
        coords_list = sorted(set([h.get_genome_range() for h in hits]))
        return ", ".join([str(c) for c in coords_list])

def _count_amplicon_exons(primer_design, target_transcript):
    amp_coords = primer_design.amplicon_trans_coords()
    cnt = 0
    for exon in target_transcript.features.iter_type(ExonFeature):
        if exon.trans.overlaps(amp_coords):
            cnt += 1
    return cnt

def output_target_debug(outdir, primer_targets, primer_designs):
    target_transcript = primer_targets.transcripts[0]
    fileOps.ensureDir(outdir)
    with open(_get_out_path(outdir, primer_targets.target_id, "debug.txt"), 'w') as fh:
        primer_targets.dump(fh)
        print("annotated_amplicon:", primer3_annotate_amplicon(primer_targets.transcripts[0]), file=fh)
        primer3_dump_args(fh, target_transcript)
        print(file=fh)
        primer_designs.primer3_results.dump(fh)
        print(file=fh)
        primer_designs.dump(fh)
        print(file=fh)

def output_target_beds(outdir, primer_targets, primer_designs):
    fileOps.ensureDir(outdir)
    _write_beds(build_target_beds(primer_targets),
                _get_out_path(outdir, primer_targets.target_id, "target.bed"))
    _write_beds(build_primer_beds(primer_designs),
                _get_out_path(outdir, primer_targets.target_id, "primers.bed"))
    if primer_designs.uniqueness_checked:
        _write_beds(build_genome_uniqueness_hits_beds(primer_designs),
                    _get_out_path(outdir, primer_targets.target_id, "genome-uniqueness.bed"))
        _write_beds(build_transcriptome_uniqueness_hits_beds(primer_designs),
                    _get_out_path(outdir, primer_targets.target_id, "transcriptome-uniqueness.bed"))

_design_tsv_header = ("target_id", "design_status", "transcript_id", "browser",
                      "primer_id", "left_primer", "right_primer", "pri",
                      "amplicon_len", "amplicon_exons", "left_delta_G", "right_delta_G",
                      "on_target_trans", "off_target_trans",
                      "on_target_genome", "off_target_genome",
                      "amplicon")

def _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls):
    "write one design, if primer_design is None, it means there were no primers found"
    row = [primer_designs.primer_targets.target_id,
           primer_designs.status, primer_designs.primer_targets.transcripts[0].trans_id]
    if first:
        row.append(_make_design_browser_link(primer_designs, hub_urls))
    else:
        row.append('')
    if primer_design is None:
        row += 12 * ['']
    else:
        row += [primer_design.ppair_id,
                primer_design.primer3_pair.PRIMER_LEFT_SEQUENCE,
                primer_design.primer3_pair.PRIMER_RIGHT_SEQUENCE,
                primer_design.priority,
                primer_design.amplicon_length,
                _count_amplicon_exons(primer_design, primer_designs.target_transcript),
                primer_design.primer3_pair.PRIMER_LEFT_END_STABILITY,
                primer_design.primer3_pair.PRIMER_RIGHT_END_STABILITY,
                _make_uniqeness_hits_browser_coords(primer_design.transcriptome_on_targets),
                _make_uniqeness_hits_browser_coords(primer_design.transcriptome_off_targets),
                _make_uniqeness_hits_browser_coords(primer_design.genome_on_targets),
                _make_uniqeness_hits_browser_coords(primer_design.genome_off_targets),
                primer_design_amplicon(primer_design, primer_designs.primer_targets.transcripts[0])]
    fileOps.prRow(fh, row)

def _write_primer_designs(fh, primer_designs, hub_urls):
    fileOps.prRow(fh, _design_tsv_header)
    first = True
    if len(primer_designs.designs) == 0:
        _write_primer_pair_design(fh, primer_designs, None, first, hub_urls)
    else:
        for primer_design in primer_designs.designs:
            _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls)
            first = False

def output_primer_designs(outdir, primer_targets, primer_designs, hub_urls):
    fileOps.ensureDir(outdir)
    with fileOps.AtomicFileCreate(_get_out_path(outdir, primer_targets.target_id, "designs.tsv")) as tmp_tsv:
        with open(tmp_tsv, "w") as fh:
            _write_primer_designs(fh, primer_designs, hub_urls)

def output_target_designs(outdir, primer_targets, primer_designs, hub_urls=None):
    """output primer TSV, BEDs, and debug information for one target.  The
    $target_id.design.tsv file is created atomically, so it can be used as a
    marker that this design is complete.  If hub_urls is specified, the
    design.tsv file will contain links to the target region in the UCSC
    browser.
    """

    output_target_debug(outdir, primer_targets, primer_designs)
    output_target_beds(outdir, primer_targets, primer_designs)
    output_primer_designs(outdir, primer_targets, primer_designs, hub_urls)
