# -*- coding: utf-8 -*-
"""Utilities for visualization."""
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

from ..utils import IOTools

import os

def save_fig(path, parent_file=None, dpi=300, ftype="png"):
    """Save your figures.
    """
    from matplotlib import pyplot as plt

    here = os.path.abspath('.')

    fn = sanitize_file_path(path)
    fn = '.'.join([fn, ftype])
    fn = os.path.join(here, fn)

    if parent_file is not None:
        fn = os.path.join(parent_file, fn)

    plt.savefig(fn, dpi=dpi)

    print("Your file has been saved: \n", fn, "\n@", dpi, "dpi")

    return

# if '__name__' == '__main__':
#     print(IOTools.sanitize_file_path('imports are working'))
