# Input and output data formats used by Primers JuJu

## Definitions and conventions

* *symbolic id* - string used to identify starting with an alphabetic character,
  with remaining characters being alpha-numeric, or one '_', '-', '.', '=', '%', or '+'
* column names are all valid Python identifiers to allow easy access from code; however, sometimes clunky names.  These start with an alphabetic character, followed by an alpha-numeric or '_" character

## Primer target specification

The *primer targets* file is a TSV file, typically exported from a spreadsheet
that specifies one or more primer targets.  Multiple targets maybe are included
in a single primer targets file.

The file consists of rows that define the genomic attributes of the primer pair
along with one target transcript.  Continuation rows are used to specify more
target transcripts.

* *target_id* - a unique symbolic id for this primer design.  It should be
  unique across all primers designed for the project.  It is suggested that this
  be based on gene symbols.  If multiple primer targets are being defined for a gene,
  it is suggested to use '+' to separate the name gene name from another symbol to
  identify the particular target.
* *region_5p* - genomic coordinates of the 5' region, in the form of valid UCSC browser
  display coordinates, such as chrX:15,547,624-15,602,400.  Note what the browser
  displays are one-based coordinates.  Commas are optional.
* *region_3p* - genomic coordinates of the 3' region.
* *trans_target_track* - symbolic name of track that contains the target transcript. A list of valid
  tracks names and their corresponding file will be defined for the project.
* *trans_target_id* - target transcript id, which is the BED name.  This must be
  unique within the track.


Continuation rows consist of the same *target_id*, with the *region_5p* and
*region_3p* cells empty.  They containing *trans_target_track* and
*trans_target_id* for each additional transcript.  No order is required, as the
*target_id* connects the rows.  However, having the continuation row follow
The primary row is good for maintaining sanity.

All other columns in TSV are passed through.  Other column names starting with
*trans_* will be associated with that specific transcripts, while others are
associated with the entire primer design.  This is significant because metadata
for browser tracks are on a per transcript basis.  The *trans_* values will
only be displayed for that transcript, while all other columns are displayed for all
transcripts.


# TODO: need examples
