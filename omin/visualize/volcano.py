# -*- coding: utf-8 -*-
"""Volcano Plotting Tools
"""

# import omin
# import pandas as pd
import re
import numpy as np

try:
    from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
except Exception as err:
    print(err)

try:
    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except Exception as err:
    print(err)

swatch = {"gray": (0.35, 0.3, 0.3),
          "opti-black": (0.10, 0.05, 0.15)}

cmap = LinearSegmentedColormap.from_list('mycmap',
                                         ["beige", swatch["opti-black"]],
                                         N=2)


def by_compartment(lfc, pval, mask, aspect=None, title=None):

    # LOAD VARIBLES
    p1 = pval.loc[mask]
    p1 = p1.reindex(index=mask.index)

    p2 = pval.loc[~mask]
    p2 = p2.reindex(index=mask.index)

    y1 = -np.log10(p1)
    y2 = -np.log10(p2)

    x1 = lfc.loc[mask]
    x1 = x1.reindex(index=mask.index)

    x2 = lfc.loc[~mask]
    x2 = x2.reindex(index=mask.index)

    # On matplotlib titles; http://matplotlib.org/users/text_intro.html
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)

    if title is not None:
        ax.set_title(title)

    # PLOT 1
    plt.scatter(x1, y1,
                c=swatch["gray"],
                edgecolors=swatch["opti-black"],
                zorder=10)
    # PLOT2
    plt.scatter(x2, y2,
                c="w",
                edgecolors=swatch["opti-black"],
                zorder=9)

    plt.ylim([0, 5])
    plt.xlim([-3, 3])
    # Show grid
    plt.grid(True)
    # Set aspect
    if aspect is not None:
        plt.axes().set_aspect(aspect)

    # Label the axes
    plt.xlabel("Log2 Fold Change", fontname="arial")
    plt.ylabel('-Log10(p-value)', fontname="arial")

    return

###

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

def scale(lfc=None, pval=None, scalar=None, compartment_mask=None,
                          true_color=None, false_color=None, ax=None):
    """Return a scaled scatter subplot.
    """

    color_matrix = boolean_color(mask=compartment_mask,
                                 true_color=true_color,
                                 false_color=false_color)

    if ax is None:
        ax = plt.gca()

    ax.scatter(lfc,
               -np.log10(pval),
               sizes=scatter_scale(scalar),
               edgecolors=[.5, .5, .5 , .7],
               color=color_matrix)

    # # Label the axes
    ax.set_xlabel("Log2 Fold Change", fontname="arial")
    ax.set_ylabel('-Log10(p-value)', fontname="arial")

    return ax
