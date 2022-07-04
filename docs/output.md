
## output BEDs

### target track

BED file, $target_id.target.bed, has the target regions that were specified as input.

Color coding:
* blue - input target regions
* green - targeted transcript, ampliction is drawn thick

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
