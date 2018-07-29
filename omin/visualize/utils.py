# -*- coding: utf-8 -*-
"""Utilities for visualization."""
# LICENSE
# -------

# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import numpy as np
from pandomics import pandas as pd

from ..utils import IOTools

try:
    from matplotlib import pyplot as plt

except Exception as err:
    print(err)


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


def volcano_latex(title=None, ax=None, *args, **kwargs):
    """Returns matplotlib plot with latex labels."""

    if ax is None:
        ax = plt.gca()

    # f, ax = plt.subplots(1, 1, figsize=(10,5))

    latex_title = '$\mathrm{{{}}}$'.format(title)
    # plt.title(latex_title, fontsize=18)
    ax.set_title(latex_title, fontsize=18)

    latex_xlabel = '$\mathrm{Log_{2}FC}$'
    # plt.xlabel(latex_xlabel, fontsize=16)
    ax.set_xlabel(latex_xlabel, fontsize=16)

    latex_ylabel = '$\mathrm{Log_{10}(p value)}$'
    # plt.ylabel(latex_ylabel, fontsize=16)
    ax.set_ylabel(latex_ylabel, fontsize=16)

    return ax.plot(*args, **kwargs)


def compartmentalize(boolean_series,
                     true_face_color=None,
                     true_edge_color=None,
                     false_face_color=None,
                     false_edge_color=None,
                     fc_name="face_color",
                     ec_name="edge_color",
                     axis=1):

    """Returns a DataFrame with the folloing columns ..."""

    if true_face_color is None:
        true_face_color = [.19, .15, .11, .9]

    if true_edge_color is None:
        true_edge_color = [.6, .25, .25, .25]

    if false_face_color is None:
        false_face_color = [.97, .981, .9, .50]

    if false_edge_color is None:
        false_edge_color = [.35, .35, .35, .8]

    face_color_series = boolean_color(boolean_series, true_face_color, false_face_color, name=fc_name)
    edge_color_series = boolean_color(boolean_series, true_edge_color, false_edge_color, name=ec_name)

    results = concat([face_color_series, edge_color_series], axis=axis)

    return results

def compartmentalize_and_scale(boolean_series,
                               scale_series,
                               max_size=50,
                               axis=1,
                               **kwargs):
    """Returns DataFrame for volcano plots."""

    compartment = compartmentalize(boolean_series=boolean_series, **kwargs)
    scalar = scatter_scale(scale_series, max_size=max_size)
    results = pd.concat([compartment, scalar.rename('scalar')], axis=axis)
    return results


def mitocarta_scale(dataframe):
    return compartmentalize_and_scale(dataframe["MitoCarta2_List"], dataframe["p_adjusted"])
