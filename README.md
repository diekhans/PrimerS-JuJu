# PrimerS-JuJu

Tools to assist with designing multiple RT-PCR primers for full-length
isoform amplification.

Package is currently alpha-quality and requires a bit manual setup to run.

This package builds on [primer3](https://github.com/primer3-org/primer3)
and
[Primer3-py](https://libnano.github.io/primer3-py/index.html).

> Untergasser, Andreas, et al. "Primer3—new capabilities and interfaces."
> Nucleic acids research 40.15 (2012): e115-e115.
> doi: 10.1093/nar/gks596

Licensed under GPL v2 to comply with the primer3 license.

# Install in a python virtualenv

    CFLAGS=-std=c99 pip install git+https://github.com/diekhans/PrimerS-JuJu.git

primer3-py needs -std=c99, at least on some systems.

# Addition requirements

* available from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
  * bigBedNamedItems
  * bigBedToBed
  
* available from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/
  * gfPcr
  * gfServer  (if not using an assembly in the UCSC browser)



