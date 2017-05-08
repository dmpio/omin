"""Omin: Omincs Modeling Integrating Normalization.

Omin aims to be a general omics tool for Python. Curently functionality is
limited to the normalization of searched LC-MSMS TMT Proteomics data. However
we aim to gradually to broaden the scope of this software as time, interest and
funding permit.

"""
from __future__ import division, print_function
# import os
import sys
# import subprocess
# import textwrap

# FIXME: Create standard issue python setup install
# FIXME: Make base compatible with Python 2.7 and up.

DOCLINES = (__doc__ or '').split("\n")

if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[:2] < (3, 4):
    raise RuntimeError("Python version 2.7 or >= 3.4 required.")

# Version information
# Multiple versions of the code base have existed since September of 2016
# However Versioning began with 0.0.0 on April 12 2017
MAJOR               = 0
MINOR               = 0
MICRO               = 3
ISRELEASED          = False
VERSION             = '{}.{}.{}'.format(MAJOR, MINOR, MICRO)
# print(VERSION)

def get_version_info():
    """Return Omin version as string

    Returns
    -------
    FULLVERSION : str
    """
    FULLVERSION = VERSION
    return FULLVERSION
