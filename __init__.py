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
FIXME: Add content.

Available subpackages
---------------------
core
    Provides core classes.
databases
    Locally stored databases and tools to read and convert them.
ExampleData
    Examples of input databases.
export
    Tools for explorting processed data.
stats
    Statistical methods.

Utilities
---------
__version__
    omin version string
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

from __future__ import division, absolute_import, print_function
import os
import sys
import warnings
from .utils import StringTools
from .utils import SelectionTools
from .utils import IOTools
from . import core
from .core.handles import *
from . import normalize
from . import stats
from .stats import Compare
from . import visualize
from . import export
from . import databases
from .databases import MitoCarta
from . import gui
from .gui import OminNotebook as nb
from . import pathfinder

if sys.version_info[0] >= 3:
    from builtins import bool, int, float, complex, object, str
    unicode = str
else:
    from __builtin__ import bool, int, float, complex, object, unicode, str

__docformat__ = 'restructuredtext'

# __dir_path__ = os.path.dirname(os.path.realpath(__file__))
