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

import os
import sys
# from __future__ import division, absolute_import, print_function
# import warnings
from . import utils
from .utils import StringTools
from .utils import SelectionTools
from .utils import IOTools
from . import core
from .core.handles import *
from . import normalize
from . import stats
from . import visualize
from . import export
from . import databases
from . import gui

if sys.version_info[0] >= 3:
    from builtins import bool, int, float, complex, object, str
    unicode = str
else:
    from __builtin__ import bool, int, float, complex, object, unicode, str

__docformat__ = 'restructuredtext'


__dir_path__ = os.path.dirname(os.path.realpath(__file__))
