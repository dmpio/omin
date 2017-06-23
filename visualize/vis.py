# -*- coding: utf-8 -*-
"""
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
import pandas as pd
import re
import matplotlib
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import numpy as np


class Volcano(object):

    swatch = {"gray": (0.35, 0.3, 0.3),
              "opti-black": (0.10, 0.05, 0.15)}

    @classmethod
    def simple(cls, FC, pvals, cutoff=-np.log10(.05), aspect=None,
               c=None, edgecolors=None):

        """Takes p-values and log2 fold changes and returns a basic volcano plot
        figure.

        Parameters
        ----------
        pvals : DataFrame
        FC : DataFrame
        aspect : int
        cutoff : float

        Returns
        -------
        [Matplotlib Figure]

        """
        if c is None:
            # c = (.3, .3, .3)
            c = cls.swatch["gray"]
        if edgecolors is None:
            # edgecolors = (.1, .1, .1)
            edgecolors = cls.swatch["opti-black"]
        y1 = -np.log10(pvals)
        x1 = FC
        # PLOT
        plt.scatter(x1, y1, c=c, edgecolors=edgecolors, zorder=10)
        plt.axhspan(cutoff, -np.log10(pvals.min()) + 2, facecolor='y', alpha=0.25)
        plt.ylim([0, int(y1.max()) + 2])
        plt.xlim([-5, 5])
        plt.grid(True)
        # LABELS
        plt.xlabel("Log2 Fold Change", fontname="arial")
        plt.ylabel('-Log10 of P-Value', fontname="arial")

        if aspect is not None:
            plt.axes().set_aspect(aspect)
        return plt

    @classmethod
    def by_compartment(cls, lfc, pval, bdex, wdex, cutoff=-np.log10(.05),
                       aspect=None):
        """Volcano plot mitochondrial peptides(black) non-mitochondrial (white).

        Parameters
        ----------
        lfc : DataFrame
        pval : DataFrame
        bdex : DataFrame
        wdex : DataFrame
        cutoff : float
        aspect : float

        Returns
        -------
        matplotlib plot

        """
        # LOAD VARIBLES
        p1 = pval.ix[bdex.index]
        p2 = pval.ix[wdex.index]
        y1 = -np.log10(p1)
        x1 = lfc.ix[bdex.index]
        y2 = -np.log10(p2)
        x2 = lfc.ix[wdex.index]

        # PLOT 1
        plt.scatter(x1, y1,
                    c=cls.swatch["gray"],
                    edgecolors=cls.swatch["opti-black"],
                    # color=swatch["opti-black"],
                    zorder=10)
        # PLOT2
        plt.scatter(x2, y2,
                    c="w",
                    edgecolors=cls.swatch["opti-black"],
                    # color=swatch["opti-black"],
                    zorder=9)

        # Show cutoff
        plt.axhspan(cutoff,
                    -np.log10(pval.min()) + 10,
                    facecolor='y',
                    alpha=0.25)

        # Set plot limits
        plt.ylim([0, 5])
        plt.xlim([-3, 3])
        # Show grid
        plt.grid(True)
        # Set aspect
        if aspect is not None:
            plt.axes().set_aspect(aspect)
        # Label the axes
        plt.xlabel("Log2 Fold Change", fontname="arial")
        plt.ylabel('-Log10 of P-Value', fontname="arial")
        # return plt
        return

    @classmethod
    def by_significance(cls, lfc, pvals, padj,
                        aspect=None, c=None, edgecolors=None,
                        y_axis=None, x_axis=None):

        if c is None:
            c = cls.swatch["gray"]

        if edgecolors is None:
            edgecolors = cls.swatch["opti-black"]

        y1 = pvals[padj.reject]
        y1 = -np.log10(y1)
        x1 = lfc[padj.reject]

        y2 = pvals[~padj.reject]
        y2 = -np.log10(y2)
        x2 = lfc[~padj.reject]
        if y_axis is None:
            y_axis = [0, int(y1.max()) + 2]

        if x_axis is None:
            x_axis = [-5, 5]

        # PLOT
        plt.scatter(x1, y1, c=c, edgecolors=edgecolors, zorder=10)
        plt.scatter(x2, y2, c=(.5, .5, .5), edgecolors=edgecolors, alpha=.8)

        plt.ylim(y_axis)
        plt.xlim(x_axis)
        plt.grid(True)
        # LABELS
        plt.xlabel("Log2 Fold Change", fontname="arial")
        plt.ylabel('-Log10 of P-Value', fontname="arial")

        if aspect is not None:
            plt.axes().set_aspect(aspect)
        return



def volcan(FC, pvals, cutoff=-np.log10(.05), aspect=None):
    """Takes p-values and log2 fold changes and returns a basic volcano plot
    figure.

    Parameters
    ----------
    pvals : DataFrame
    FC : DataFrame
    aspect : int
    cutoff : float

    Returns
    -------
    [Matplotlib Figure]

    """
    y1 = -np.log10(pvals)
    x1 = FC
    # PLOT
    plt.scatter(x1, y1, c=(.3, .3, .3), edgecolors=(.1, .1, .1), zorder=10)
    plt.axhspan(cutoff, -np.log10(pvals.min()) + 2, facecolor='y', alpha=0.25)
    plt.ylim([0, int(y1.max()) + 2])
    plt.xlim([-5, 5])
    plt.grid(True)
    # LABELS
    plt.xlabel("Log2 Fold Change", fontname="arial")
    plt.ylabel('-Log10 of P-Value', fontname="arial")

    if aspect is not None:
        plt.axes().set_aspect(aspect)
    return plt


def plotByMito(lfc, pval, bdex, wdex, cutoff=-np.log10(.05), aspect=None):
    """Creates a volcano plot that plots mitochondrial proteins as black and the
     non-mitochondrial proteins as white.

    Parameters
    ----------
    lfc : DataFrame
    pval : DataFrame
    bdex : DataFrame
    wdex : DataFrame
    cutoff : float
    aspect : float

    Returns
    -------
    matplotlib plot

    """
    # LOAD VARIBLES
    p1 = pval.ix[bdex.index]
    p2 = pval.ix[wdex.index]

    y1 = -np.log10(p1)
    x1 = lfc.ix[bdex.index]

    y2 = -np.log10(p2)
    x2 = lfc.ix[wdex.index]

    # DEFINE RGB COLORS
    swatch = {"gray": (0.35, 0.3, 0.3),
              "opti-black": (0.10, 0.05, 0.15)}

    # PLOT 1
    plt.scatter(x1, y1,
                c=swatch["gray"],
                edgecolors=swatch["opti-black"],
                # color=swatch["opti-black"],
                zorder=10)
    # PLOT2
    plt.scatter(x2, y2,
                c="w",
                edgecolors=swatch["opti-black"],
                # color=swatch["opti-black"],
                zorder=9)

    # Show cutoff
    plt.axhspan(cutoff, -np.log10(pval.min()) + 2, facecolor='y', alpha=0.25)
    # Set plot limits
    plt.ylim([0, 5])
    plt.xlim([-3, 3])
    # Show grid
    plt.grid(True)
    # Set aspect
    if aspect is not None:
        plt.axes().set_aspect(aspect)
    # return plt
    return


def trePlot(trunob):
    """Takes trun of replicate object and returns between genotype comparisons.

    Notes
    -----
        FIXME : This function needs to be generalized

    Parameters
    ----------
    trunob : (:obj)

    Returns
    -------
    matplotlib plot
    """
    bdex = trunob.mitopep[trunob.mitopep.MitoCarta2_List == 1.0]
    wdex = trunob.mitopep[trunob.mitopep.MitoCarta2_List != 1.0]

    xtix = ytix = 7
    # Subplot1
    ax1 = plt.subplot(131)
    plotByMito(trunob.lfc.ix[:, 0], trunob.pval.ix[:, 0], bdex, wdex)
    plt.ylabel('-Log10 of P-Value', fontname="arial")
    plt.xlabel('Log2 Fold Change', fontname="arial")
    plt.xticks(size=xtix)
    plt.yticks(size=ytix)
    # Subplot2
    ax2 = plt.subplot(132, sharex=ax1, sharey=ax1)
    plotByMito(trunob.lfc.ix[:, 1], trunob.pval.ix[:, 1], bdex, wdex)
    plt.xlabel('Log2 Fold Change', fontname="arial")
    plt.xticks(size=xtix)
    plt.setp(ax2.get_yticklabels(), visible=False)

    # Subplot3
    ax3 = plt.subplot(133, sharex=ax1, sharey=ax1)
    plotByMito(trunob.lfc.ix[:, 2], trunob.pval.ix[:, 2], bdex, wdex)
    plt.xlabel('Log2 Fold Change', fontname="arial")
    plt.xticks(size=xtix)
    plt.setp(ax3.get_yticklabels(), visible=False)

    # set limits
    plt.xlim(-3.0, 3.0)
    plt.ylim(0, 5)
    # adjust width
    plt.subplots_adjust(wspace=0.08)
    # plt.show()
    return plt


def quadPlot(trunob):
    """For within genotype comparisons.

    Notes
    -----
        FIXME: This function needs be generalized.

    Parameters
    ----------
    trunob : (:obj)

    Returns
    -------
    matplotlib plot

    """
    bdex = trunob.mitopep[trunob.mitopep.MitoCarta2_List == 1.0]
    wdex = trunob.mitopep[trunob.mitopep.MitoCarta2_List != 1.0]
    xtix = ytix = 7
    # Subplot1 Upper Left
    ax1 = plt.subplot(221)
    plotByMito(trunob.elfc_ko.ix[:, 0], trunob.epval_ko.ix[:, 0], bdex, wdex)
    plt.ylabel('-Log10 of P-Value', fontname="arial")
    plt.xticks(size=xtix)
    plt.yticks(size=ytix)
    # Subplot2 Upper Right
    ax2 = plt.subplot(222, sharex=ax1, sharey=ax1)
    plotByMito(trunob.elfc_ko.ix[:, 1], trunob.epval_ko.ix[:, 1], bdex, wdex)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.xticks(size=xtix)
    plt.yticks(size=ytix)
    # Subplot3 Lower Left
    ax3 = plt.subplot(223, sharex=ax1, sharey=ax1)
    plotByMito(trunob.elfc_wt.ix[:, 0], trunob.epval_wt.ix[:, 0], bdex, wdex)
    plt.ylabel('-Log10 of P-Value', fontname="arial")
    plt.xlabel('Log2 Fold Change', fontname="arial")
    plt.xticks(size=xtix)
    plt.yticks(size=ytix)
    # Subplot4 Lower Right
    ax4 = plt.subplot(224, sharex=ax1, sharey=ax1)
    plotByMito(trunob.elfc_wt.ix[:, 1], trunob.epval_wt.ix[:, 1], bdex, wdex)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.xlabel('Log2 Fold Change', fontname="arial")
    plt.xticks(size=xtix)
    plt.yticks(size=ytix)
    # set limits
    plt.xlim(-3.0, 3.0)
    plt.ylim(0, 5)
    # adjust width
    plt.subplots_adjust(wspace=0.08, hspace=0.25)
    return plt


def saveFig(fig_title, parent_file=None, dpi=300, ftype=".png"):
    """Takes a figure and figure title and saves it with a file system
    compatable name.

    Parameters
    ----------
    figtitle : str
    parent_file : str or None
        Defaults to None
    dpi : int
        Defaults to 300
    ftype : str
        Defaults to '.png'

    """
    from matplotlib import pyplot as plt

    if parent_file is None:
        # Make fig file name
        fn = re.sub('[^0-9a-zA-Z]+',  # regex pattern
                    '_',  # replacement
                    fig_title)  # input string
        fn = fn + ftype
        plt.savefig(fn, dpi=dpi)
        print("Your file has been saved in this directory with the title: ",
              fn, "@", dpi, "dpi")
    else:
        fn = re.sub('[^0-9a-zA-Z]+',  # regex pattern
                    '_',  # replacement
                    fig_title)  # input string
        fn = fn + ftype
        fpath = parent_file + "/" + fn
        plt.savefig(fpath, dpi=dpi)
        print("Your figure has been saved in " + parent_file +
              " directory with the title: ", fn, "@", dpi, "dpi")

# VENN DIAGRAMS


def mitoAboveCut(ligob, bdex, cond, cut=None):
    """Finds the interesting peptides in a DataFrame.

    Parameters
    ----------
    ligob: DataFrame
    bdex : DataFrame
        DataFrame index of mitocarta calls inside mitochondria
    cond : str
        Experimental condition as string.
    cut : str or None
        Can be 'Q1' for quadrant 1, 'Q4' for quadrant 4 or None which will
        return everything above the

    Returns:
        Dataframe of peptides located in mitochondria at specified cutoff.
    """
    if cut is None:
        pv = omin.sep(ligob.pval, cond).ix[:, 0].sort_values()
        pv = pv[pv < .05]
        mitocom = pv[bdex.index]
        return mitocom.dropna()
    if cut == 'Q1':
        pv = omin.sep(ligob.pval, cond).ix[:, 0].sort_values()
        lf = omin.sep(ligob.lfc, cond).ix[:, 0].sort_values(ascending=False)
        com = pd.concat([lf[lf > 0], pv[pv < .05]], axis=1)
        mitocom = com.ix[bdex.index]
        return mitocom.dropna()
    if cut == 'Q4':
        pv = omin.sep(ligob.pval, cond).ix[:, 0].sort_values()
        lf = omin.sep(ligob.lfc, cond).ix[:, 0].sort_values(ascending=False)
        com = pd.concat([lf[lf < 0], pv[pv < .05]], axis=1)
        mitocom = com.ix[bdex.index]
        return mitocom.dropna()


def grabTre(ligob, bdex, condlist):
    """
    Parameters
    ----------
    ligob : (:obj)
    bdex : DataFrame
    condlist : list

    Returns
    -------
    grab_list : list

    """
    grab_list = [mitoAboveCut(ligob, bdex, i) for i in condlist]
    return grab_list


def vennTre(ligob, bdex, condlist):
    """
    Notes
    -----
        Needs better documentation.

    Parameters
    ----------
    ligob : (:obj)
    bdex : DataFrame
    condlist : list

    Returns
    -------
    matplotlib_venn
    """
    gt = grabTre(ligob, bdex, condlist)
    return venn3([set(i.index) for i in gt], condlist)


def vennDuo(ob1, ob2, bdex, cond, cut=None):
    """
    Parameters
    ----------
    ob1 : (:obj)
    ob2 : (:obj)
    bdex : DataFrame
    cond : str
    cut : int

    Returns
    -------
    matplotlib_venn

    """
    midex1 = mitoAboveCut(ob1, bdex, cond, cut)
    midex2 = mitoAboveCut(ob2, bdex, cond, cut)
    return venn2([set(midex1.index), set(midex2.index)])


def vennRepVsRep(ob1, ob2, bdex, condlist, cut=None):
    """

    Notes
    -----
        Needs better documentation.

    Parameters
    ----------
    ob1 : (:obj)
    ob2 : (:obj)
    bdex : DataFrame
    condlist : list
    cut : int or None

    Returns
    -------
    matplotlib plot
    """
    ax1 = plt.subplot(131)
    v1 = vennDuo(ob1, ob2, bdex, condlist[0], cut)
    v1.get_label_by_id('A').set_text('Rep1')
    v1.get_label_by_id('B').set_text('Rep2')
    plt.title("Rest")

    ax2 = plt.subplot(132)
    v2 = vennDuo(ob1, ob2, bdex, condlist[1], cut)
    v2.get_label_by_id('A').set_text('Rep1')
    v2.get_label_by_id('B').set_text('Rep2')
    plt.title("10 Post")

    ax3 = plt.subplot(133)
    v3 = vennDuo(ob1, ob2, bdex, condlist[2], cut)
    v3.get_label_by_id('A').set_text('Rep1')
    v3.get_label_by_id('B').set_text('Rep2')
    plt.title("60 Post")
    return plt


def vennConditionsWithInGeno(trunch_object,
                             modification_object,
                             condition_list):
    """
    Parameters
    ----------
    trunch_object : (:obj)
    modification_object : (:obj)
    condition_list : list

    Returns
    -------
    matplotlib_venn

    """
    comp_set = [set(omin.allComp(trunch_object, modification_object,
                                 i).index.tolist()) for i in condition_list]
    return venn2([comp_set[0], comp_set[1]])
