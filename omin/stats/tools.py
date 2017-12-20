# -*- coding: utf-8 -*-
"""Statistical tools"""


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

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy.stats import f_oneway
import functools
from statsmodels.sandbox.stats.multicomp import multipletests


class Compare(object):
    """Tools for comparisons."""

    @staticmethod
    def log2FC(numer, denom, new_column_name=""):
        """Take the log2 fold change of normalized DataFrames of simillar size.

        Parameters
        ----------
        numer : DataFrame
            The numerator DataFrame
        denom : DataFrame
            The denominator DataFrame
        new_column_name : str
            Include your the name of comparison. Defaults to a blank string.

        Returns
        -------
        lfc : DataFrame

        """
        if len(new_column_name) > 0:
            new_column_name = " "+new_column_name
        lfc = numer.mean(axis=1) - denom.mean(axis=1)
        lfc = pd.DataFrame(lfc,
                           columns=["LFC"+new_column_name],
                           index=numer.index)
        return lfc

    @staticmethod
    def ttester(numer, denom, new_column_name=""):
        """For pvalue comparision of types of DataFrames of similar shape.

        Notes
        -----
        Make sure that your DataFrames are the same size. Future versions
        Should rely on this method as a base for pvalue comparison.

        Parameters
        ----------
        numer : DataFrame
            The numerator DataFrame
        denom : DataFrame
            The denominator DataFrame

        Returns
        -------
        pvals : DataFrame

        Examples
        --------
        Comparing two DataFrames of similar size.
        >>>omin.ttester(KO_DataFrame,WT_DataFrame)

        """
        if len(new_column_name) > 0:
            new_column_name = " "+new_column_name
        # The loop below suppresses an irrelevent error message.
        # For more details on this see:
        # http://stackoverflow.com/questions/40452765/invalid-value-in-less-when-comparing-np-nan-in-an-array
        with np.errstate(invalid='ignore'):
            np.less([np.nan, 0], 1)
            # ttest_ind implemented
            pvals = ttest_ind(numer, denom, axis=1).pvalue
        pvals = pd.DataFrame(pvals,
                             columns=["pval"+new_column_name],
                             index=numer.index)
        return pvals

    @staticmethod
    def anova_oneway(numer, denom, new_column_name=""):
        """For pvalue comparision of types of DataFrames of similar shape.

        Parameters
        ----------
        numer : DataFrame
            The numerator DataFrame
        denom : DataFrame
            The denominator DataFrame

        Returns
        -------
        pvals : DataFrame

        """
        if len(new_column_name) > 0:
            new_column_name = " "+new_column_name
        # The loop below suppresses an irrelevent error message.
        with np.errstate(invalid='ignore'):
            np.less([np.nan, 0], 1)
            pvals = f_oneway(numer.T, denom.T).pvalue
        pvals = pd.DataFrame(pvals,
                             columns=["pval"+new_column_name],
                             index=numer.index)
        return pvals

    @staticmethod
    def bh_fdr(p_val, alpha=None):
        """Return adjusted p-values.

        statsmodels.sandbox.stats.multicomp.multipletests function with
        method = "fdr_bh" default and alpha=.05

        Parameters
        ----------
        p_val : DataFrame

        Returns
        -------
        p_adj : DataFrame
        """
        alpha = alpha or .05
        # FIXME : Add somekind of dataprovenonce measure here.
        bh_funct = functools.partial(multipletests,
                                     method="fdr_bh",
                                     alpha=alpha)

        p_adj = bh_funct(pvals=p_val.dropna().values.T[0])
        # p_adj = bh_funct(pvals=p_val.values.T[0])
        p_adj = pd.DataFrame([p_adj[0], p_adj[1]]).T
        p_adj.index = p_val.dropna().index
        p_adj = p_adj.reindex(index=p_val.index)
        p_adj.columns = ["reject", "p_adjusted"]
        # write over NaNs will boolean False
        p_adj.reject = p_adj.reject.fillna(False)
        return p_adj

    @classmethod
    def comparison_summary_stats(cls, numerator, denominator, cut_off=None,
                                 mitodex=None, nonmitodex=None):
        """Summarize the stats of a comparison."""
        # Define the significance cut off
        cut_off = cut_off or -np.log10(.05)
        pval = cls.ttester(numerator, denominator)
        lfc = cls.log2FC(numerator, denominator)
        over_cutoff = -np.log10(pval).dropna() > cut_off
        over_cutoff = over_cutoff.reindex(index=numerator.index).fillna(False)
        total_sig = over_cutoff.values.sum()
        print("Total peptides with -Log10(p-values) > {:.3}:".format(cut_off),
              total_sig)
        # significant LFCs.
        sig_lfc = lfc.ix[over_cutoff.ix[:, 0]]
        # Find all negative significant fold changes.
        negative_sig_lfc = sig_lfc < 0
        negative_sig_lfc = sig_lfc.ix[negative_sig_lfc.ix[:, 0]]
        print("\t Total negative LFCs:", negative_sig_lfc.shape[0])

        # Find all positive significant fold changes.
        positive_sig_lfc = sig_lfc > 0
        positive_sig_lfc = sig_lfc.ix[positive_sig_lfc.ix[:, 0]]
        print("\t Total positive LFCs:", positive_sig_lfc.shape[0])

        # number of significant mitochondrial peptides
        mito_over_cutoff = -pval.reindex(index=numerator.index).ix[mitodex.index].apply(np.log10)
        mito_over_cutoff = mito_over_cutoff.dropna() > cut_off
        mito_sig = mito_over_cutoff.values.sum()
        print("Total mitochondrial peptides with -Log10(p-values) > {:.3}:".format(cut_off),
              mito_sig)

        # number of significant nonmitochondrial peptides
        nonmito_over_cutoff = -pval.reindex(index=numerator.index).ix[nonmitodex.index].apply(np.log10)
        nonmito_over_cutoff = nonmito_over_cutoff.dropna() > cut_off
        nonmito_sig = nonmito_over_cutoff.values.sum()
        print("Total of non-mitochondrial peptides with -Log10(p-values) > Log10(.05):",
              nonmito_sig)

    @classmethod
    def annotated_comparison(cls, numerator, denominator, master_index,
                             fdr_alpha=None):
        """Return an annotated comparison as a DataFrame."""
        pvl = cls.ttester(numerator, denominator)
        lfc = cls.log2FC(numerator, denominator)
        qvl = cls.bh_fdr(pvl, alpha=fdr_alpha)
        output = pd.concat([master_index, lfc, pvl, qvl], axis=1)
        return output
