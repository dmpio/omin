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

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

# FIXME: Add built in python bhfdr method in here.

def log2FC(numer, denom, new_column_name=""):
    """Takes the log2 fold change of normalized DataFrames of simillar size.

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
    lfc = pd.DataFrame(lfc, columns=["LFC"+new_column_name], index=numer.index)
    return lfc


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
        #ttest_ind implemented
        pvals = ttest_ind(numer, denom, axis=1).pvalue
    pvals = pd.DataFrame(pvals, columns=["pval"+new_column_name], index=numer.index)
    return pvals