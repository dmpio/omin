import omin
import pandas as pd


def setDiff(list_A, list_B):
    """Takes difference of two lists with respect to list_A. Simillar to `list(set(list_A)-set(list_B))` however
    duplicates are not deleted and order is preserved.

    Args:
        list_A (list): A list of elements.
        list_B (list): A list of elements.
    Return:
    """
    return [i for i in list_A if i not in set(list_B)]

class CompOb:
    """
    """
    def __init__(self,list_A,list_B,how=None,modification_labels = ["Acetyl","Phospho"]):
        """
        """
        #If how is not specified then the index is used
        if how == None:
            self.justA = list_A.ix[setDiff(list_A.index,list_B.index)]
            self.justB = list_B.ix[setDiff(list_B.index,list_A.index)]
            self.AandB = list_A.ix[set(list_A.index)&set(list_B.index)]
        else:
            self.justA = list_A.ix[list_A[how].isin(setDiff(list_A[how],list_B[how]))]
            self.justB = list_B.ix[list_B[how].isin(setDiff(list_B[how],list_A[how]))]
            combo = pd.concat([list_A,list_B],keys=modification_labels)
            self.AandB = combo.ix[combo[how].isin(set(list_A[how])&set(list_B[how]))]




# def aboveCut(trunch_object, cond):
#     """Selects p-values above a given cutoff in a given trunch_object and returns a DataFrame of p-values and log fold
#     changes.
#
#     Parameters
#     ----------
#     trunch_object:
#
#     cond:
#
#     Returns
#     -------
#
#     """
#     pv = omin.sep(trunch_object.pval, cond).ix[:, 0].sort_values()
#     pv = pv[pv<.05]
#     lfc = omin.sep(trunch_object.lfc,cond).ix[pv.index]
#
#     return pd.concat([pv,lfc],axis=1)

def aboveCut(trunch_object,cond,pval_kind="pval",lfc_kind="lfc",cut=.05):
    """Selects p-values above a given cutoff in a given trunch_object and returns a DataFrame of p-values and log fold
    changes.

    Parameters
    ----------
    trunch_object:

    cond:

    Returns
    -------

    """
    pv = omin.sep(trunch_object.__dict__[pval_kind], cond).ix[:, 0].sort_values()
    pv = pv[pv<cut]
    lfc = omin.sep(trunch_object.__dict__[lfc_kind],cond).ix[pv.index]

    return pd.concat([pv,lfc],axis=1)


def allComp(trunch_object, modification_object, cond,pval_kind="pval",lfc_kind="lfc",cut=.05):
    """

    TODO: Include modification columns

    Parameters
    ----------
    trunch_object
    modification_object
    cond

    Returns
    -------

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
