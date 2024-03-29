"""Visualization tools for clustering.

"""

import numpy as np

try:
    from matplotlib import pyplot as plt
except Exception as err:
    print(err)

try:
    from scipy.spatial import distance
    from scipy.cluster.hierarchy import linkage
    from scipy.cluster.hierarchy import dendrogram
except Exception as err:
    print(err)



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
    dendro_subplot = plt.subplot2grid((10, 10), (0, 1),
                                      rowspan=1, colspan=9)

    # Create a linkage object
    # Dataframe most be transposed for this.
    linkmat = linkage(distance.pdist(dataframe.T))

    # Make a dendrogram from the linkage object.
    den0 = dendrogram(linkmat)


    if fig_title is not None:
        fig.suptitle(fig_title)

    # Remove all axes for dendrogram
    plt.axis("off")

    # ax2 = plt.subplot2grid((10, 10), (1, 0), rowspan= 4, colspan=1)
    heatmap_subplot = plt.subplot2grid((10, 10), (1, 1),
                                       rowspan=9, colspan=9)

    # Generate heatmap
    # plt.pcolor(dataframe, cmap = cmap)
    heatmap_subplot.pcolormesh(dataframe.values[:, den0["leaves"]], cmap=cmap)

    # Set dataframe column names as the axis labels
    if labels is not None:
        plt.xticks(np.arange(0.5, dataframe.shape[1], 1))
        plt.gca().set_xticklabels(labels,rotation=270)
    else:
        try:
            plt.xticks(np.arange(0.5, dataframe.shape[1], 1))
            plt.gca().set_xticklabels(list(dataframe.columns[den0["leaves"]]),
                                      rotation=270)
        except Exception:
            print("Not possible to add labels.")
            pass

    # Insert the colorbar for scale.
    # plt.colorbar()
    # cb = plt.colorbar(ax=dendro_subplot)
    # cb.ax.set_visible(False)

    # Make the vertical distance between plots equal to zero
    plt.subplots_adjust(hspace=0)
