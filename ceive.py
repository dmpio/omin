# -*- coding: utf-8 -*-
import omin
import pandas as pd
import numpy as np
###FILTERING FUNCTIONS###
# ----------------------------------------------------------------------------------------------------------------------
def specSel(dataframe,include_list,exclude_list,case=True):
    """Select columns whose headers contain items just the items you want.

    Parameters
    ----------
    dataframe : DataFrame
    include_list : list
    exclude_list : list

    Returns
    -------
    selected : DataFrame

    """
    allbool = np.ones(len(dataframe.columns),dtype=bool)
    for mod in include_list:
            trubool = dataframe.columns.str.contains(mod,case)
            allbool = allbool&trubool
    for ex in exclude_list:
        fakbool = dataframe.columns.str.contains(ex,case)
        allbool = allbool&~fakbool

    if allbool.any():
        selected = dataframe[dataframe.columns[allbool]]
    else:
        print("Special select failed. Try something else.")
        selected = np.nan
    return selected

def setDiff(list_A, list_B):
    """Takes difference of two lists with respect to list_A. Simillar to `list(set(list_A)-set(list_B))` however
    duplicates are not deleted and order is preserved.

    Parameters
    ----------
    list_A : list
    list_B : list

    Returns
    -------
    the_diff :  list
        List of the difference between sets.
    """
    the_diff = [i for i in list_A if i not in set(list_B)]

    return the_diff

class CompOb:
    """Make a comparison object that mirrors the output of a Venn diagram.
    Attributes
    ----------
    justA : DataFrame
        The elements exclusive to list_A.
    justB : DataFrame
        The elements exclusive to list_B.
    AandB : DataFrame
        The elements found in list_A and list_B.
    """
    def __init__(self,list_A,list_B,how=None,modification_labels = ["Acetyl","Phospho"]):
        """
        Parameters
        ----------
        list_A : list
        list_B : list
        how : str
            Can be compared by any column as a string. If no string is specified then the index is used.
        modification_labels : list
            If none is specified then it defaults to  ["Acetyl", "Phospho"]
        """
        if how == None:
            self.justA = list_A.ix[setDiff(list_A.index,list_B.index)]
            self.justB = list_B.ix[setDiff(list_B.index,list_A.index)]
            self.AandB = list_A.ix[set(list_A.index)&set(list_B.index)]
        else:
            self.justA = list_A.ix[list_A[how].isin(setDiff(list_A[how],list_B[how]))]
            self.justB = list_B.ix[list_B[how].isin(setDiff(list_B[how],list_A[how]))]
            combo = pd.concat([list_A,list_B],keys=modification_labels)
            self.AandB = combo.ix[combo[how].isin(set(list_A[how])&set(list_B[how]))]

def aboveCut(trunch_object,cond,pval_kind="pval",lfc_kind="lfc",cut=.05):
    """Selects p-values above a given cutoff in a given trunch_object and returns a DataFrame of p-values and log fold
    changes.

    Parameters
    ----------
    trunch_object : (:obj)
    cond : str
    pval_kind : str
        If none is specified it defaults to "pval"
    lfc_kind : str
        If none is specified it defaults to "lfc"
    cut : float
        If none is specified it defaults to .05

    Returns
    -------
    above_cut : DataFrame
        Contains p-values and log fold change columns.
    """
    pv = omin.sep(trunch_object.__dict__[pval_kind], cond).ix[:, 0].sort_values()
    pv = pv[pv<cut]
    lfc = omin.sep(trunch_object.__dict__[lfc_kind],cond).ix[pv.index]
    above_cut = pd.concat([pv,lfc],axis=1)
    return above_cut

def allComp(trunch_object, modification_object, cond,pval_kind="pval",lfc_kind="lfc",cut=.05):
    """
    Notes
    -----
        TODO: Include modification columns

    Parameters
    ----------
    trunch_object : (:obj)
    modification_object : (:obj)
    cond : str
    pval_kind : str
        If none is specified it defaults to "pval"
    lfc_kind : str
        If none is specified it defaults to "lfc"
    cut : float
        If none is specified it defaults to .05

    Returns
    -------
    compound : DataFrame

    """
    # Grab the P-values and LFCs
    pv_lfc = omin.aboveCut(modification_object, cond, pval_kind,lfc_kind,cut)
    # Grab the master protein accession and gene id
    mpa_gn = trunch_object.mpa.ix[pv_lfc.index]
    #Grab the Modifications columns
    mods = trunch_object.peptides.Modifications[pv_lfc.index]
    # grab the MitoCarta2.0 data
    mito = trunch_object.mitopep.ix[pv_lfc.index][["MitoCarta2_List", "Matrix", "IMS"]]
    # Concatenate DataFrames
    compound = pd.concat([mpa_gn, mods,pv_lfc, mito], axis=1)
    # Change NaNs for zeros
    compound = compound.fillna(0)
    return compound
