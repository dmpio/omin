# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re
import omin

def logNormToAve(pepdf):
    """Takes a DataFrame composed of a fraction of peptide abundances and then subtracts each element in each row by the sum of it's row.

    Notes
    -----
    FIXME : Make sure this method could possibly handle protein data aswell.

    Parameters
    ----------
    pepdf : DataFrame
        Limit to a single fraction of peptide abundance data.

    """
    log2 = pepdf.copy().apply(np.log2)
    ave = log2.mean(axis=1)
    log2_div_ave = log2.sub(ave, axis=0)
    log2_div_ave.columns = "Log2-AVE: " + log2_div_ave.columns
    return log2_div_ave

def normToPool(log2_div_ave):
    """Normalizes a DataFrame composed of a single fraction of peptide data that contains a single 'Control' or 'Pool' column.

    Notes
    -----
    See also: logNormToAve

    Parameters
    ----------
    log2_div_ave : DataFrame

    """
    pool = omin.betSep(log2_div_ave, "control", "pool")[0]
    lda_div_pool = log2_div_ave.sub(pool.ix[:, 0], axis=0)
    # Rename the columns to reflect the operations on them.
    lda_div_pool.columns = [re.sub("Log2-AVE", "Log2-AVE-Pool", i) for i in lda_div_pool.columns]
    return lda_div_pool

###COMPARISON TOOLS###
#======================================================================================================================

def log2FC(numer, denom, new_column_name = ""):
    """Takes the log2 fold change of normalized DataFrames of simillar size.

    Parameters
    ----------
    numer : DataFrame
        The numerator DataFrame
    denom : DataFrame
        The denominator DataFrame
    new_column_name : str
        Include your the name of your comparison here. Defaults to a blank string.

    Returns
    -------
    lfc : DataFrame

    """
    if len(new_column_name) > 0:
        new_column_name = " "+new_column_name
    lfc = numer.mean(axis=1) - denom.mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["LFC"+new_column_name],index=numer.index)
    return lfc

def ttester(numer, denom,new_column_name = ""):
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

    """
    from scipy.stats import ttest_ind
    if len(new_column_name) > 0:
        new_column_name = " "+new_column_name
    pvals = ttest_ind(numer, denom, axis=1).pvalue
    pvals = pd.DataFrame(pvals, columns=["pval"+new_column_name], index=numer.index)
    return pvals
