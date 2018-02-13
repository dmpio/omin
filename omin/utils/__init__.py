# -*- coding: utf-8 -*-
"""Utilites for the omin module."""

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

# FIXME: Add try and except to most if not all functions. NO MORE QUITE FAILS!
# FIXME: DEPRECATE OLD FUNCTIONS
# FIXME: Create tools to search against against databases the user specifies.

# from .utils import *
from .filter_tools import FilterTools
from .fasta_tools import FastaTools
from .fasta_tools import UniProtTools
from .string_tools import StringTools
from .selection_tools import SelectionTools
from .io_tools import IOTools
from .pd_tools import PDStudyTools
from .object_tools import inspectObject
from .object_tools import objectWalker
from .warnings import deprecated
