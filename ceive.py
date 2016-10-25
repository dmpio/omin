# -*- coding: utf-8 -*-
import omin
import pandas as pd
import numpy as np
#---------------------------------------------------------------------------------------------------------------------
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

def manyModSel(pepdf, terms):
    """Returns searched peptide a tuple of searched DataFrames with a given modifications or modifications.

    Parameters
    ----------
    pepdf : DataFrame
        With peptides information.
    terms : list
        Can be any number of modifications as a string. Case does not matter and regex special characters can be
        used e.g. 'acetyl', 'Phospho',hydroxy...methyl.glutaryl,'ect'

    Returns
    -------
    selected : tuple
        Entering more than one term last element of the tuple will contain all modified peptides.

    """
    selected = ()
    for term in terms:
        term = omin.mod_dict[term]
        moddex = pepdf.Modifications.str.contains(pat=term, case=False)
        if moddex.sum() > 0:
            selected += (pepdf.ix[moddex],)
            print(moddex.sum(), "peptides with", term, "modification found.")
        else:
            print("No peptides with", term, "modification were found.")
            pass
    if len(selected) > 1:
        all_select = np.bitwise_or.reduce([df.index for df in selected])
        all_selected = pepdf.ix[all_select]
        selected = selected + (all_selected,)
    return selected

def vLook(peptides,proteins,mods):
    """Returns a tuple of selected peptides and proteins.

    Takes raw peptides and protiens returns a tuple of selected peptides and proteins. The function can also select for a sigle
    modification or many modifications.

    Parameters
    ----------
    peptides : DataFrame
    proteins : DataFrame
    mods : list

    Returns
    -------
    peptide_select : DataFrame
    protein_select : DataFrame

    Examples
    --------
    >>>mpa_pep,fdr_prot = vLook(raw.peptides,raw.proteins)
    >>>mpa_pep,fdr_prot = vLook(raw.peptides,raw.proteins,"hydroxy...methyl.glutaryl")
    >>>mpa_pep,fdr_prot = vLook(raw.peptides,raw.proteins,"Acetyl","Phospho")

    See Also
    --------
    manyModSel
    masterOne
    masterPep

    """

    fdr = omin.masterOne(proteins)
    if len(mods) == 0:
        mpa = omin.masterPep(peptides)
    else:
        mpa = omin.masterPep(omin.manyModSel(peptides,mods)[-1])
    fdrdf = pd.DataFrame(fdr.Accession,index = fdr.index)
    peptide_select = mpa.merge(fdrdf, on ="Accession",how="left",right_index=True)
    protein_select = mpa.merge(fdrdf, on ="Accession",how="left",left_index=True)
    return peptide_select,protein_select

###VENN DIAGRAM FUNCTIONS###
# ---------------------------------------------------------------------------------------------------------------------

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
