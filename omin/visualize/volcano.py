# -*- coding: utf-8 -*-
"""Volcano Plotting Tools

Copyright 2017 James Draper, Paul Grimsrud, Deborah Muoio

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files, Omics Modeling Integrating
Normalization (OMIN), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM.
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# import omin
# import pandas as pd
import re
import matplotlib
from matplotlib import pyplot as plt
# from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap

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
