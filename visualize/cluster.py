"""Visualization tools for clustering.

LICENSE:

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

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt


def dendro_heatmap(dataframe, figsize=None, cmap=None, fig_title=None,
                   labels=None):
    """Create a vertical heatmap with a dendrogram for data.

    Dendrogram is created with scipy.cluster.hierarchy.linkage and
    scipy.cluster.hierarchy.dendrogram. Heatmap is created with
    matplotlib.pyplot.pcolor.

    Parameters
    ----------
    dataframe : Dataframe
        Comparisons will be made across columns.
    figsize : tuple
        Defaults to (5,7)
    cmap : str
        Can be any of the standard matplotlib colormaps.
    fig_title : str
        Tiltle the figure.
    labels : list
        If None will attempt to use the DataFrame column headers.

    Returns
    -------
    fig : figure
        A matplotlib figure.

    See Also
    --------
    scipy.cluster.hierarchy.linkage
    scipy.cluster.hierarchy.dendrogram
    matplotlib.pyplot.pcolor
    """
    # Set the heatmap if None
    cmap = cmap or "jet"

    # Create the figure and set it's size.
#     figsize = figsize or (5, 7)
    figsize = figsize or (8.5, 11)
    fig = plt.figure(figsize=figsize)

    # Create the subplot for the dendrogram
    dendro_subplot = plt.subplot2grid((10,10), (0,1), rowspan = 1, colspan = 9)

    # Create a linkage object. The transposed dataframe needs to be use for this part
    linkmat = linkage(dataframe.T)

    # Make a dendrogram from the linkage object.
    dendro_subplot = dendrogram(linkmat)

    if fig_title is not None:
        fig.suptitle(fig_title)

    # Remove all axes for dendrogram
    plt.axis("off")

    # ax2 = plt.subplot2grid((10, 10), (1, 0), rowspan= 4, colspan=1)
    heatmap_subplot = plt.subplot2grid((10, 10), (1, 1), rowspan = 9, colspan=9)

    # Generate heatmap
    # plt.pcolor(dataframe, cmap = cmap)
    plt.pcolormesh(dataframe, cmap = cmap)

    # Set dataframe column names as the axis labels
    if labels is not None:
        plt.xticks(np.arange(0.5, dataframe.shape[1], 1))
        plt.gca().set_xticklabels(labels,rotation=270)
    else:
        try:
            plt.xticks(np.arange(0.5, dataframe.shape[1], 1))
            plt.gca().set_xticklabels(list(dataframe.columns), rotation=270)
        except Exception:
            print("Not possible to add labels.")
            pass

    # Make the vertical distance between plots equal to zero
    plt.subplots_adjust(hspace=0)

    # Insert the colorbar for scale.
    # plt.colorbar()
    # cb = plt.colorbar(heatmap_subplot)
    # cb.ax.set_visible(False)
    return fig
