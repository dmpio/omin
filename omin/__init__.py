# -*- coding: utf-8 -*-
"""
Omin.

====

Provides
  1. Rapid normalization of proteomics data.
  2. State of the art visualiazion tools for proteomics data.
  3. Interactive web ready objects for investigation of omics data.

How to use the documentation
----------------------------
FIXME: ADD BRIEF OVERVIEW

Available subpackages
---------------------

utils
    Utilites for manaing omics data.
core
    Provides core classes.
databases
    Locally stored databases and tools to read and convert them.
visualize
    Data graphing tools.
export
    Tools for explorting processed data.
stats
    Statistical methods.


Utilities
---------
__version__
    omin version string

DEVELOPMENT PROTIPS
===================

PROTIP: Define variables used at the lowest possible class level.

"""

# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach,
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

# BOILERPLATE

import os
import sys
import warnings

# # PYTHON COMPATIBILITY
# from __future__ import division, absolute_import, print_function
# if sys.version_info[0] >= 3:
#     from builtins import bool, int, float, complex, object, str
#     unicode = str
# else:
#     from __builtin__ import bool, int, float, complex, object, unicode, str

# UTILS IMPORTS
from .utils import StringTools
from .utils import SelectionTools
from .utils import IOTools
from .utils import UniProtTools

# CORE IMPORTS
from . import core
from .core.handles import *

# STATS IMPORTS
from . import stats
from .stats import Compare

# VISUALIZE IMPORTS
from . import visualize

# EXPORT IMPORTS
from . import export

# DATABASE IMPORTS
from . import databases
from .databases import MitoCarta
from .databases import MitoCartaTwo

# GUI IMPORTS
from . import gui
from .gui import OminNotebookController as nb

# CLI IMPORTS
from . import cli

# FIXME: DEPRECATE THE FOLLOWING
from . import pathfinder
from . import normalize

here = os.getcwd()

with open(os.path.join(here, 'omin', '__version__')) as f:
    __version__ = f.read()


__docformat__ = 'restructuredtext'
