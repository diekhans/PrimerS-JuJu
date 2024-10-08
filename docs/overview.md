# Overview of PrimersJuJu workflow

Primers-JuJu takes manually selected regions and isoforms and use
primer3 with other validation to produce a set of primer pairs.

The command-line program primers-juju takes user-defined input sets
and produces a set of prioritized primer sequence pairs.  The results are both
a TSV file with information about the primers and tracks in a UCSC browser
assembly hub.

## Input

For each desired primer pair, the user provides as input:

* a unique symbolic name to identify the primer set being designed
* the genome assembly
* list of transcript models to use
* transcript set to check for off-target primers, currently this must be part of the target model set.
* coordinates ranges of two regions from which to pick the primer pairs
* parameters to use for the primer design

The genome assembly is specified as the genome assembly in UCSC two-bit
file.  The transcript models are specified as one or more tracks and the
transcript ids within those tracks.

The coordinate range pairs must meet the primer3 requirement that at
least one of the two sequences is less than 60 bp in length.  If a region
spans an intron, then. 

The input file specification is described in
[Input and output data formats used by PrimerS JuJu](data_formats.md).

## Processing

* Validate input to make sure it is consistent with the requirements:
** Allow for slight fuzziness in region boundaries and adjust to match exons
** Check intron-spanning requirements.
* Call primer3 to produce a list of possible primer pairs.
* Double-check delta-G for secondary structures of any self-dimers, hairpins, and heterodimers that should be weaker (more positive) than -9.0 kcal/mole.
* Cs and Gs (three hydrogen bonds) are preferred at the primer terminus (3’) to enhance specific binding.
* check for off-target mappings using isPcr to check for the uniqueness of each primer pair in both the transcriptome and genome
* quality control to check against other isoforms

## Output

A TSV file is created that contains the prioritized list of primer pairs.
The file includes the genomic coordinates, all attributes reported by primer3,
the results of isPcr, and the sequences.

BED files are produced for:
** browser track of original regions
** browser track with primer pairs
** browser track with amplicons

Color coding is used to identify the highest priority primer pairs.

A track hub is created.

## Possible extensions
* if a region exceeds the primer3 maximum, split it up into multiple calls to primer3

