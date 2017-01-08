# -*- coding: utf-8 -*-
import re
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
# import omin.ceive


def logger(d):
    """ Takes a DataFrame and returns a logspace object.
    Notes
    -----
    This function is deprecated in favor of the class Logger

    Parameters
    ----------
    d : DataFrame

    Returns
    -------
    gout : (:obj)
        logspace object.

    """
    dout = d.apply(np.log2)
    dout.columns = "Log2 " + dout.columns
    ave = pd.DataFrame(dout.mean(axis=1), columns=["AVE"])
    gout = logspace(dout, ave)
    return gout

def logdiv(lso):
    """
    Parameters
    ----------
    lso : (:obj)
        logspace object.

    Returns
    -------
    dout : DataFrame
        log2-AVE
    """
    dout = lso.log.as_matrix() - lso.ave.as_matrix()
    dout = pd.DataFrame(dout, index=lso.log.index)
    colnm = [re.sub("Log2", "Log2-AVE", i) for i in lso.log.columns]
    dout.columns = colnm
    return dout

def dfMin(d1, d2):
    """ Take the difference of two DataFrames.


    """
    return pd.DataFrame(d1.as_matrix() - d2.as_matrix(), index=d1.index)

def treatAve(dataframe, treat):
    return pd.DataFrame(sep(dataframe, treat).mean(axis=1),
                        columns=["AVE " + treat],
                        index=dataframe.index)

def treatSTD(d, treat):
    return pd.DataFrame(sep(d, treat).std(axis=1),
                        columns=["STDEV " + treat],
                        index=d.index)

def treatAveSTD(d, numer, denom):
    avet1 = treatAve(d, numer)
    avet2 = treatAve(d, denom)
    stdt1 = treatSTD(d, numer)
    stdt2 = treatSTD(d, denom)
    return pd.concat([avet1, avet2, stdt1, stdt2], axis=1)

def logFC(d, num, dem):
    """
    Takes:
    d = DataFrame
    num = treatment
    """
    dout = treatAve(d, num).as_matrix() - treatAve(d, dem).as_matrix()
    dout = pd.DataFrame(dout,
                        columns=["LogFC(" + num + "/" + dem + ")"],
                        index=d.index)
    return dout

def logFolder(ko, wt):
    """

    Parameters
    ----------
    ko : DataFrame
    wt : DataFrame

    Returns
    -------
    lfc : DataFrame

    """
    lfc = ko.mean(axis=1) - wt.mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["lfc"])
    return lfc

def pvalr(d, A, B):
    """
    d = DataFrame
    A = term 1
    B = term 2
    """
    dout = ttest_ind(sep(d, A), sep(d, B), axis=1).pvalue
    dout = pd.DataFrame(dout, columns=["pval " + A + " vs. " + B], index=d.index)
    return dout

def pvaln(numer, denom):
    """For pvalue comparision of types of DataFrames of similar shape.

    Note: Make sure that your DataFrames are the same size. Future versions
    Should rely on this method as a base for pvalue comparison.
    Args:
        numer (:obj): DataFrame numer stands for numerator
        denom(:obj): DataFrame denom stands for denominator

    return:
        pvals (:obj): DataFrame of the pvalues
    """
    pvals = ttest_ind(numer, denom, axis=1).pvalue
    pvals = pd.DataFrame(pvals, columns=["pval"], index=numer.index)
    return pvals


def overPooler(dataframe):
    """Devides DataFrame by pool.

    Parameters
    ----------
    dataframe : DataFrame

    Returns
    -------
    dopool : DataFrame
        You dataframe devided by the pool
    """
    logoverave = logdiv(logger(dataframe))
    # c,s = sepCon(logoverave,"Control")#FIXME
    c, s = betSep(logoverave, "Control", "Pool")  # FIXED!
    dopool = dfMin(logoverave, c)
    colnm = [re.sub("Log2-AVE", "Log2-AVE/pool", i) for i in logoverave.columns]
    dopool.columns = colnm
    return dopool

def tooLogFC(num, dem, search_term):
    """Takes the Log2 fold change for two DataFrames.

    Parameters
    ----------
    num : DataFrame
    dem : DataFrame
    search_term : str

    Returns
    -------
    lfc : DataFrame

    """
    numop = overPooler(num)
    demop = overPooler(dem)
    lfc = sep(numop, search_term).mean(axis=1) - sep(demop, search_term).mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["LogFC " + search_term])
    return lfc

def tooPvalr(num, dem, search_term):
    """Takes the p-value for two DataFrames.

    Parameters
    ----------
    num : DataFrame
    dem : DataFrame
    search_term : str

    Returns
    -------
    pval : DataFrame

    """
    numop = overPooler(num)
    demop = overPooler(dem)
    pval = ttest_ind(sep(numop, search_term), sep(demop, search_term), axis=1, nan_policy="omit").pvalue
    pval = pd.DataFrame(pval, index=num.index, columns=["Pval " + search_term])
    return pval

def preLogFC(numop, demop, search_term):
    """Takes the Log2 fold change If you have already divided by pool.

    Parameters
    ----------
    numop : DataFrame
    demop : DataFrame
    search_term : str

    Returns
    -------
    lfc : DataFrame

    """
    lfc = sep(numop, search_term).mean(axis=1) - sep(demop, search_term).mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["LogFC " + search_term])
    return lfc

def prePvalr(numop, demop, s):
    """Takes the p-value if you have already divided by pool.

    Parameters
    ----------
    numop : DataFrame
    demop : DataFrame
    search_term : str

    Returns
    -------
    lfc : DataFrame

    """
    pval = ttest_ind(sep(numop, s), sep(demop, s), axis=1, nan_policy="omit").pvalue
    pval = pd.DataFrame(pval, index=numop.index, columns=["Pval " + s])
    return pval
