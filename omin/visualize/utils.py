# -*- coding: utf-8 -*-
"""Utilities for visualization."""
# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import numpy as np
from ..utils import IOTools

try:
    from matplotlib import pyplot as plt

except Exception as err:
    print(err)


def save_fig(path, parent_file=None, dpi=300, ftype="png"):
    """Save your figures.
    """

    here = os.path.abspath('.')

    fn = IOTools.sanitize_file_path(path)
    fn = '.'.join([fn, ftype])
    fn = os.path.join(here, fn)

    if parent_file is not None:
        fn = os.path.join(parent_file, fn)

    plt.savefig(fn, dpi=dpi)

    print("Your file has been saved: \n", fn, "\n@", dpi, "dpi")

    return


def boolean_color(mask, true_color='black', false_color='white',
                  name="face_color"):
    """Return a Series with colors replacing True and False."""
    _mask = mask.copy()
    result = _mask.apply(lambda x: true_color if x == True else false_color)
    result.name = name
    return result


def scatter_scale(target, max_size=50.0):
    """Return a scaled list of sizes."""
    result = target.apply(lambda x:max_size*(1-x))
    return result
