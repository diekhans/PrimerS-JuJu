# Input and output data formats used by Primers JuJu

## Definitions and conventions

* *symbolic id* - string used to identify starting with an alphabetic character,
  with remaining characters being alpha-numeric, or one '_', '-', '.', '=', '%', or '+'
* column names - these are all valid Python identifiers, to allow easy access from code, however sometimes clunky names. These start with an alphabetic character, followed by a alpha-numeric or '_" character

## Primer target specification

The *primer targets* file is a TSV file, normally exported from a spreadsheet,
that specified one or more primer targets.  Multiple targets maybe be included
in a single primer targets file.

The file consists of rows that define the genomic attributes of the primer pair
along with one target transcript.  Continuation rows are used to specify more
target transcripts.

* *target_id* - a unique symbolic id for this primer design.  It should be
  unique across all primers design for the project.  It is suggested that this
  be based on gene symbol.  If multiple primer targets are being defined for a gene,
  it is suggest using '+' to separate the name gene name from another symbol to
  identify the particular target.
* *region_5p* - genomic coordinates of the 5' region, it the form of valid UCSC browser
  display coordinates, such as chrX:15,547,624-15,602,400.  Note what the browser
  displays are one-based coordinates.  Commas are optional.
* *region_3p* - genomic coordinates of the 3' region.
* *target_trans_track* - symbolic name of track that contains the target transcript. A list of valid
  tracks names and their corresponding file will be defined for the project.
* *target_trans_id* - target transcript id, which is the BED name.  This must be
  unique withing the track.


Continuation rows consist of the same *target_id*, with the *region_5p* and
*region_3p* cells empty.  They containing *target_trans_track* and
*target_trans_id* for each additional transcript.  No order is required, as the
*target_id* connects the rows.  However having the continuation row follow
the primary row is good for maintained sanity.
