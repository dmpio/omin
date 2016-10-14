import re
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

def sep(dataframe_in, search_term):
    """Takes DataFrame and search_term and returns a new DataFrame that contains columns that contain that search_term.

    Parameters
    ----------
    dataframe_in : DataFrame
    search_term : str
        What kind of columns do you want to find.

    Returns
    -------
    dataframe_out : DataFrame
        A DataFrame that with columns contain just search_term.

    """
    dataframe_out = None
    if dataframe_in.columns[dataframe_in.columns.str.contains(search_term, case=False)].any():
        dataframe_out = dataframe_in[dataframe_in.columns[dataframe_in.columns.str.contains(search_term, case=False)]].copy()
    else:
        print("The DataFrame has no columns that contain:", search_term)
    return dataframe_out

def modSel(dataframe_in, mod1, mod2):
    """Selects peptides with modifications

    Parameters
    ----------
    dataframe_in : DataFrame
    mod1 : str
    mod2 : str

    Returns
    -------
    dmod12 : DataFrame
        Contains peptides with both mod1 and mod2.
    dmod1 : DataFrame
        Contains peptides with mod1
    dmod2 : DataFrame
        Contains peptides with mod1

    """
    dmod1 = dataframe_in[dataframe_in.Modifications.str.contains(mod1)]
    dmod2 = dataframe_in[dataframe_in.Modifications.str.contains(mod2)]
    dmod12 = dataframe_in[dataframe_in.Modifications.str.contains(mod1) ^ dataframe_in.Modifications.str.contains(mod2, case=False)]
    return dmod12, dmod1, dmod2


def sepCon(dataframe_in, separation_term):
    """Separates dataframe_in two DataFrames one with the separation_term and one without the separation_term.

    Parameters
    ----------
    dataframe : DataFrame
    separation_term : str

    Returns
    -------
    with_term = DataFrame
    without_term = DataFrame

    """
    with_term = dataframe_in[dataframe_in.columns[dataframe_in.columns.str.contains(separation_term, case=False)]]
    without_term = dataframe_in[dataframe_in.columns[~dataframe_in.columns.str.contains(separation_term, case=False)]]
    return with_term, without_term

def betSep(dataframe_in, *args):
    """A better version of the function sepCon.

    Parameters
    ----------
    dataframe_in : DataFrame
    *args : str

    """
    sel = np.array([dataframe_in.columns.str.contains(i, case=False) for i in args])
    sel = np.any(sel, axis=0)
    with_terms = dataframe_in[dataframe_in.columns[sel]]
    without_terms = dataframe_in[dataframe_in.columns[~sel]]
    return with_terms,without_terms

class logspace:
    """
    Notes
    -----
    This class is deprecated in favor of the class Logger.

    Attributes
    ----------
    log : DataFrame
    ave : DataFrame
    com : DataFrame
    """
    def __init__(self, d, ave):
        """
        Parameters
        ----------
        d : DataFrame
        ave : DataFrame

        """
        self.log = d
        self.ave = ave
        # Combine the DataFrames
        self.com = pd.concat([self.log, self.ave], axis=1)

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

def treatAve(d, treat):
    return pd.DataFrame(sep(d, treat).mean(axis=1),
                        columns=["AVE " + treat],
                        index=d.index)

def treatSTD(d, treat):
    return pd.DataFrame(sep(d, treat).std(axis=1),
                        columns=["STDEV " + treat],
                        index=d.index)

def treatAveSTD(d, t1, t2):
    avet1 = treatAve(d, t1)
    avet2 = treatAve(d, t2)
    stdt1 = treatSTD(d, t1)
    stdt2 = treatSTD(d, t2)
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
    lfc = ko.mean(axis=1) - wt.mean(axis=1)
    return pd.DataFrame(lfc, columns=["lfc"])

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


def overPooler(d):
    logoverave = logdiv(logger(d))
    # c,s = sepCon(logoverave,"Control")#FIXME
    c, s = betSep(logoverave, "Control", "Pool")  # FIXED!
    dopool = dfMin(logoverave, c)
    colnm = [re.sub("Log2-AVE", "Log2-AVE/pool", i) for i in logoverave.columns]
    dopool.columns = colnm
    return dopool


def tooLogFC(num, dem, s):
    numop = overPooler(num)
    demop = overPooler(dem)
    lfc = sep(numop, s).mean(axis=1) - sep(demop, s).mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["LogFC " + s])
    return lfc

def tooPvalr(num, dem, s):
    numop = overPooler(num)
    demop = overPooler(dem)
    pval = ttest_ind(sep(numop, s), sep(demop, s), axis=1, nan_policy="omit").pvalue
    pval = pd.DataFrame(pval, index=num.index, columns=["Pval " + s])
    return pval

def preLogFC(numop, demop, s):
    """If you have already divided by pool.
    """
    lfc = sep(numop, s).mean(axis=1) - sep(demop, s).mean(axis=1)
    lfc = pd.DataFrame(lfc, columns=["LogFC " + s])
    return lfc

def prePvalr(numop, demop, s):
    """If you have already divided by pool.
    """
    pval = ttest_ind(sep(numop, s), sep(demop, s), axis=1, nan_policy="omit").pvalue
    pval = pd.DataFrame(pval, index=numop.index, columns=["Pval " + s])
    return pval
