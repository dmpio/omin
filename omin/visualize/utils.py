# -*- coding: utf-8 -*-
"""Utilities for visualization."""
# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import numpy as np

from ..utils import IOTools

try:
    import webcolors
except Exception as err:
    print(err)

def save_fig(path, parent_file=None, dpi=300, ftype="png"):
    """Save your figures.
    """
    from matplotlib import pyplot as plt

    here = os.path.abspath('.')

    fn = IOTools.sanitize_file_path(path)
    fn = '.'.join([fn, ftype])
    fn = os.path.join(here, fn)

    if parent_file is not None:
        fn = os.path.join(parent_file, fn)

    plt.savefig(fn, dpi=dpi)

    print("Your file has been saved: \n", fn, "\n@", dpi, "dpi")

    return


def web2rgb(color):
    """Turn your webcolor names and hex strings into rgb format."""
    try:
        return np.array(webcolors.name_to_rgb(color))/255.0
    except Exception:
        pass
    try:
        return np.array(webcolors.hex_to_rgb(color))/255.0
    except Exception:
        pass
