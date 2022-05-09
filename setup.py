#!/usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

requirements = [
    "primer3-py>=0.6.1",
    "twobitreader>=3.1.7",
    "pycbio @ git+https://github.com/diekhans/pycbio.git",
]

setuptools.setup(
    name = 'PrimersJuJu',
    version = '0.1.0',
    description = "PrimersJuJu",
    long_description = "PrimersJuJu",
    author = "Mark Diekhans",
    author_email = 'markd@ucsc.edu',
    url = 'https://github.com/diekhans/PrimerS-JuJu',
    scripts=[
    ],
    packages = [
        'primersjuju',
    ],
    package_dir = {'': 'lib'},
    include_package_data = True,  # MUST update MANIFEST.in to include files
    install_requires = requirements,
    license = "GPL2",
    zip_safe = True,
    keywords = ['Bioinformatics', 'genomics', 'transcriptomics'],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    python_requires='>=3.7',

)
