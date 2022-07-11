
## output BEDs

## design TSV and spreadsheet

* target_id - target id specified in input
* design_status - indicates the overall results of the design:
  * GOOD - primer design was successful with no off-target alignments
  * NOT_GENOME_UNIQUE - has off-target genomic alignments
  * NOT_TRANSCRIPTOME_UNIQUE - has off-target transcriptome alignments
  * NO_PRIMERS - primer3 did not find any primers
* transcript_id - transcript that was used for the design
* browser - coordinates of the target transcript, as a link to the UCSC Genome Browser that includes the generated track hub
* primer_id - unique id assign to the primer; the target_id plus _ppN.
* left_primer - left primer sequence returned by primer3
* right_primer - right primer sequence returned by primer3
* pri - priority of the primer, with 1 being the best.  Primers having no off-target alignments and the higher combined delta-G have a better priority
* amplicon_len - length of the amplicon in bases
* amplicon_exons - number of exons in the amplicon
* left_delta_G - primer3 computed stability of the left primer, higher is better
* right_delta_G - primer3 computed stability of the right primer, higher is better
* on_target_trans - on-target transcriptome alignments of the primer pair from isPcr.  These should match the designed primers
* off_target_trans - list of off-target transcriptome alignments locations of the primer pair from isPcr
* on_target_genome - on-target genome alignments of the primer pair from isPcr.  These should match the designed primers.
* off_target_genome - list of off-target genome alignments locations of the primer pair from isPcr
* annotated_amplicon - amplicon sequence with annotation for feeding manually to primer3 web.  This is for debugging and will be drop and may not be correct.


### target track

BED file, $target_id.target.bed, has the target regions that were specified as input.

Color coding:
* blue - input target regions
* green - targeted transcript, amplicion is drawn thick

# primer track

BED file, $target_id.primers.bed, of the primers returned by primer3 mapped
back to the genome using the target transcript.  These are colored based
uniqueness mappings.

Color coding:
* fuchsia - strong delta-G, (<= -9.0 kcal/mole).
* green - a primer that is only has on-target mappings when querying both the genome and transcriptome.
* orange  - a primer that has both on-target and off-target mappings when querying both the genome and transcriptome.
* red - a primer that only has off-target mappings when querying both the genome and transcriptome.
* purple - primer that is only has mappings to non-target chromosomes when querying both the genome and transcriptome.  A non-target chromosome is an alternate or patch sequence
* darkorange - primer that is no mappings when querying both the genome and transcriptome.

Primers pair has one end drawn thick, the other thin.

# uniqueness track

BED files, $target_id.transcriptome-uniqueness.bed and $target_id.genome-uniqueness.bed,
are the alignments produced by querying the transcriptome and genome with isPCR.

Color coding:
* green - an on-target alignment of a primer
* red - an off-target alignment of a primer
* purple - a alignment to a non-target chromosome (alt or patch) of a primer

Primers pair has one end drawn thick, the other thin.
