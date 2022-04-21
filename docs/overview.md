# Overview of PrimersJuJu workflow

PrimersJuJu takes manually select regions and isoforms and use
primer3 to produce a set of primer pairs.

The command line program pjuju is used take user-defined input sets and
produce a set of prioritized primer sequence pairs.

  ** generate track of region being designed
# Overview of PrimersJuJu workflow

PrimersJuJu takes manually select regions and isoforms and use
primer3 to produce a set of primer pairs.

The command line program pjuju is used take user-defined input sets and
produce a set of prioritized primer sequence pairs.

  ** generate track of region being designed

## Input

For a given locus, the user provides as input:

* a unique name to identify the primer set being designed.
* genome assembly
* a list of transcript models to use
* coordinates ranges of two regions from which to pick the primer pairs
* parameters of uses for the primer design

The genome assembly is specifics as the genome assembly in UCSC two-bit
format.  The transcript models are specified as one or more tracks and the
transcript ids within those tracks.

The coordinate range pairs must meet the primer3 requirement that that at
least one of the two sequences is less than 60 bp in length.  

Session
inteserct range with trans fuzzy

HOW TO DO THIS:
If a coordinate
range cross a splice junction, then it will be a required splice junction for
primer3.  A given region must not contain more than one splice junction.

All of the primer3 parameters maybe changed from the defaults, either by the
command line or a configuration file.

## Processing

* validate input to make sure it is consistent with the requirements
* call primer3 to produce a list of possible primer pairs.
* use isPcr to check for uniqueness of each primer pair.
* output results

## Output

A TSV file is created that contains the prioritized list of primer pairs.
This includes the genomic cooridnates, all attributes report by primer3,
the results of isPcr and the sequences.

A BED version annotating the primers in the genome.  This includes all
attributes, suitable for turning create a bigBed for a track hub.  Color
coding is used to identify the highest priority primer pairs.


## Possible extensions

* pre-filter regions with mappability tracks
* if a region exceeds the primer3 maximum, split it up into multiple calls to primer3
* split regions exceeding primer3 maximum size into overlapping regions and make multiple calls to primer3
* create separate primer pair lists when a region crosses multiple splice junctions

## Input

For a given locus, the user provides as input:

* a unique name to identify the primer set being designed.
* genome assembly
* a list of transcript models to use
* coordinates ranges of two regions from which to pick the primer pairs
* parameters of uses for the primer design

The genome assembly is specifics as the genome assembly in UCSC two-bit
format.  The transcript models are specified as one or more tracks and the
transcript ids within those tracks.

The coordinate range pairs must meet the primer3 requirement that that at
least one of the two sequences is less than 60 bp in length.  

Session
inteserct range with trans fuzzy

HOW TO DO THIS:
If a coordinate
range cross a splice junction, then it will be a required splice junction for
primer3.  A given region must not contain more than one splice junction.

All of the primer3 parameters maybe changed from the defaults, either by the
command line or a configuration file.

## Processing

* validate input to make sure it is consistent with the requirements
* call primer3 to produce a list of possible primer pairs.
* use isPcr to check for uniqueness of each primer pair.
* output results

## Output

A TSV file is created that contains the prioritized list of primer pairs.
This includes the genomic cooridnates, all attributes report by primer3,
the results of isPcr and the sequences.

A BED version annotating the primers in the genome.  This includes all
attributes, suitable for turning create a bigBed for a track hub.  Color
coding is used to identify the highest priority primer pairs.


## Possible extensions

* pre-filter regions with mappability tracks
* if a region exceeds the primer3 maximum, split it up into multiple calls to primer3
* split regions exceeding primer3 maximum size into overlapping regions and make multiple calls to primer3
* create separate primer pair lists when a region crosses multiple splice junctions
