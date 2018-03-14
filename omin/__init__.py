# -*- coding: utf-8 -*-
"""
Omin.

====

Provides:
    Rapid normalization of proteomics data. Visualiazion tools of proteomics
    data.

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
stats
    Statistical methods.


Utilities
---------
__version__
    omin version string

"""
# -------
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

# ----------
# TO DO LIST
# ----------
# FIXME: DOCUMENT OR DIE #DOD
# FIXME: ADD TRACEBACK TO ALL TRY AND EXCEPT LOOPS
# FIXME: Make a .omin config file class that points to a file in %userprofile%
# FIXME: Make a .omin/Database file class that could hold MitoCartaTwo and perhaps read and process others
# FIXME: ADD KFW EXAMPLE FILES
# FIXME: *** update the readme ** [ ] How to install,
# FIXME: *** Document dependencies *** [ ] guipyter, [ ] panomics
# FIXME: Investigate a SQLite/Json file method in APPDATA or linux eqv.
# FIXME: Store each handle class as SQLite database in same parent dir.
# FIXME: Include paragraph description of the types of filtering.

# -------
# PROTIPS
# -------

# PROTIP: Define variables used at the lowest possible class level.
# PROTIP: Also ways set your returns in a function.
# PROTIP: Add more type checking.
# PROTIP: Add more try and excepts but try to put them on function level

# ----------------
# EXTERNAL IMPORTS
# ----------------
import os, sys, warnings

# -------------
# UTILS IMPORTS
# -------------
from .utils import StringTools
from .utils import SelectionTools
from .utils import IOTools
from .utils import UniProtTools

# ------------
# CORE IMPORTS
# ------------
from . import core
from .core.handles import Process
from .core.operations import Operate

# -------------
# STATS IMPORTS
# -------------
from . import stats
from .stats import Compare

# -----------------
# VISUALIZE IMPORTS
# -----------------
from . import visualize

# ----------------
# DATABASE IMPORTS
# ----------------
from .databases import MitoCartaTwo

# -----------
# CLI IMPORTS
# -----------
from . import cli

# -------
# VERSION
# -------
def find_path():
    """Find the location of omin package in any given file system."""
    __dir_path__ = os.path.dirname(os.path.realpath(__file__))
    return __dir_path__

with open(os.path.join(find_path(), '__version__')) as f:
    __version__ = f.read().strip()

__docformat__ = 'numpy'
