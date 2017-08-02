"""Omin: Omincs Modeling Integrating Normalization.

Omin aims to be a general omics tool for Python. Curently functionality is
limited to the normalization of searched LC-MSMS TMT Proteomics data. However
we aim to gradually to broaden the scope of this software as time, interest and
funding permit.

"""

# LICENSE
# -------

# Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
# Blair Chesnut, and Elizabeth Hauser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, (the software)), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions: The above copyright
# notice and this permission notice shall be included in all copies or
# substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from __future__ import division, print_function
# import os
import sys
import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# FIXME: Create standard python setup install
# FIXME: Attempt to make base compatible with Python 2.7 and up.

DOCLINES = (__doc__ or '').split("\n")

if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[:2] < (3, 4):
    raise RuntimeError("Python version 2.7 or >= 3.4 required.")

# Version information
# Multiple versions of the code base have existed since September of 2016
# However Versioning began with 0.0.0 on April 12 2017
MAJOR = 0
MINOR = 0
MICRO = 11
ISRELEASED = False
VERSION = '{}.{}.{}'.format(MAJOR, MINOR, MICRO)
# print(VERSION)

def version():
    """Return Omin version as string

    Returns
    -------
    FULLVERSION : str
    """
    FULLVERSION = VERSION
    return FULLVERSION
