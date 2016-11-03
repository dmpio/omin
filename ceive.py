# -*- coding: utf-8 -*-
import omin
import pandas as pd
import numpy as np
import re

###Selection Functions###
# ---------------------------------------------------------------------------------------------------------------------

def sep(dataframe_in, search_term, strict=False, match=False):
    """Takes DataFrame and search_term and returns a new DataFrame that contains columns that contain that search_term.

    Parameters
    ----------
    dataframe_in : DataFrame
    search_term : str
        What kind of columns do you want to find.
    strict : bool
        Defaults to False. FIXME : Annotate this.
    match : bool
        Defaults to False. If the function will use pandas dataframe.columns.str.match which is more strict than dataframe.columns.str.search.
    Returns
    -------
    dataframe_out : DataFrame
        A DataFrame that with columns contain just search_term.

    Examples
    --------
    >>>omin.sep(mydataframe,"Search Term")
    dataframe_out
    >>>omin.sep(mydataframe,"Search Term",match=True)
    Scricter dataframe_out

    See Also
    --------
    omin.sepCon
    omin.betSep
    omin.modSel
    omin.manyModSel

    """
    dataframe_out = None
    if match:
        dataframe_out = dataframe_in[dataframe_in.columns[dataframe_in.columns.str.match(search_term)]].copy()
        return dataframe_out

    if dataframe_in.columns[dataframe_in.columns.str.contains(search_term, case=False)].any():
        if strict:
            dataframe_out = dataframe_in[dataframe_in.columns[
                np.array([search_term in set(re.sub(" ", "", i).split(",")) for i in dataframe_in.columns])]]
        else:
            dataframe_out = dataframe_in[
                dataframe_in.columns[dataframe_in.columns.str.contains(search_term, case=False)]].copy()
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

    See Also
    --------
    omin.sep
    omin.sepCon

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

    See Also
    --------
    omin.sep
    omin.betSep
    omin.modSel
    omin.manyModSel
    omin.specSel

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

    Returns
    -------
    with_terms : DataFrame
    without_terms : DataFrame

    """
    sel = np.array([dataframe_in.columns.str.contains(i, case=False) for i in args])
    sel = np.any(sel, axis=0)
    with_terms = dataframe_in[dataframe_in.columns[sel]]
    without_terms = dataframe_in[dataframe_in.columns[~sel]]
    return with_terms,without_terms

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

def sevSel(dataframe=None, term_list=None, match=False):
    """
    Parameters
    ----------
    dataframe : DataFrame
    term_list : list

    Returns
    -------
    dataframe_out : DataFrame

    """
    if type(term_list) == str:
        term_list = [term_list]
    if match:
        dataframe_out = pd.concat([omin.sep(dataframe, term, match=True) for term in term_list], axis=1)
        return dataframe_out
    else:
        dataframe_out = pd.concat([omin.sep(dataframe, term) for term in term_list], axis=1)
        return dataframe_out

def colSelPro(dataframe=None,term_list=None):
    """Returns the dataframe of the columns specified in the terms_list.

    For each term in term_list an exact match is tried before excepting when omin.sep method is used.

    Parameters
    ----------
    dataframe : DataFrame
    term_list : list

    Returns
    -------
    out_dataframe : DataFrame

    """
    df_list = []
    for term in term_list:
        try:
            selected_col = dataframe[term]
            df_list.append(selected_col)

        except:
            selected_col = omin.sep(dataframe,term)
            df_list.append(selected_col)

    out_dataframe = pd.concat(df_list,axis=1)
    return out_dataframe

def superGroup(dataframe=None,new_level=None):
    """Returns a multiindexed DataFrame with the top index named new_level.

    Parameters
    ----------
    dataframe : DataFrame
    new_level : str

    Returns
    -------
    out_df : DataFrame

    """
    if type(dataframe.columns) == pd.indexes.base.Index:
        out_df = pd.DataFrame(dataframe.values,index=dataframe.index,columns=pd.MultiIndex.from_product([[new_level],dataframe.columns]))
        return out_df
    if type(dataframe.columns) == pd.indexes.multi.MultiIndex:
        if len(dataframe.columns.levels[0])<1:
            levels = [list(i.values) for i in dataframe.columns.levels]
            levels = [[new_level]]+levels
            out_df = pd.DataFrame(dataframe.values, index = dataframe.index, columns = pd.MultiIndex.from_arrays(levels))
            return out_df
        else:
            levels = [[new_level]]+[list(i.values) for i in dataframe.columns.levels]
            labels = [list(i) for i in dataframe.columns.labels]
            new_list = list(np.linspace(0,0,len(labels[-1]),dtype=int))
            labels = [new_list]+labels
            multi = pd.MultiIndex(levels,labels)
            out_df = pd.DataFrame(dataframe.values,index=dataframe.index,columns=multi)
            return out_df

###FILTERING FUNCTIONS##
#---------------------------------------------------------------------------------------------------------------------

def masterCleanse(protein_df):
    """Filters raw protein DataFrame for master proteins.

    The raw protein data from Proteome Discoverer there is a column with the title 'Master' this funtion scans through
    that column and selects only the proteins that end with the string "IsMasterProtein"

    Parameters
    ----------
    protein_df : DataFrame
        Raw protein DataFrame

    Returns
    -------
    clean : DataFrame
        Protein DataFrame that contains only proteins with 'IsMasterProtein' in 'Master' column of protein_df
    """
    clean = protein_df.ix[protein_df.Master.str.endswith("IsMasterProtein")]
    return clean

def onePerQ(protein_df):
    """Filters raw protein DataFrame for proteins that are less than 1% the expected q-value.

    Scans through the protein DataFrame selecting only the proteins with less than 1% of the expected q-value.

    Parameters
    ----------
    protein_df : DataFrame
        Raw protein DataFrame

    Returns
    -------
    clean : DataFrame
        Protein data that contains only proteins with proteins only less than 1% of the expected q-value.
    """
    one_per = protein_df["Exp. q-value"] < .01
    one_per = protein_df.ix[one_per]
    return one_per


def masterOne(protein_df):
    """Takes a raw protein DataFrame and filters it using first the 'masterCleanse' function and 'onePerQ' function.

    Parameters
    ----------
    protein_df : DataFrame
        Raw proteins.

    Returns
    -------
    master_one : DataFrame
        Of master proteins with exp. q-value <1%
    """
    master = masterCleanse(protein_df)
    master_one = onePerQ(master)
    return master_one


def masterPep(peptide_df):
    """Takes a peptide DataFrame and returns just the first master protein accession for each peptide.

    Notes
    -----
    Assumes the first uniprot ID list is the correct one. Peptides with no master protein accession will be lost however
    the index of peptide_df will be preserved.

    Parameters
    ----------
    peptide_df : DataFrame

    Returns
    -------
    master_prot_acc : DataFrame

    """
    master_prot_acc = [i.split(';')[0] for i in peptide_df['Master Protein Accessions'].dropna()]

    master_prot_acc = pd.DataFrame(master_prot_acc,
                                   index=peptide_df['Master Protein Accessions'].dropna().index, columns=['Accession'])
    return master_prot_acc

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

    fdr = masterOne(proteins)
    if len(mods) == 0:
        mpa = masterPep(peptides)
    if len(mods) == 1:
        mpa = masterPep(manyModSel(peptides,mods)[0])
    else:
        mpa = masterPep(manyModSel(peptides,mods)[-1])
    fdrdf = pd.DataFrame(fdr.Accession,index = fdr.index)
    peptide_select = mpa.merge(fdrdf, on ="Accession",how="left",right_index=True)
    protein_select = mpa.merge(fdrdf, on ="Accession",how="left",left_index=True)
    return peptide_select,protein_select

def mitoCartaPepOut(raw_file,mods = ["Acetyl","Phospho"],dex = False):
    """
    Parameters
    ----------
    raw_file : (:obj)
        An instance of the class omin.experiment.RawData
    mods : list
        Defaults to ["Acetyl","Phospho"].
    dex : bool
        Defaults to False. When False output is mitocarta_pep if True output is a tuple containing mitodex and nonmitodex

    Returns
    -------
    (mitodex, nonmitodex) : tuple(DataFrame,DataFrame)
    mitocarta_pep : Dataframe

    Examples
    --------
    >>>mitocarta_pep = mitoCartaPepOut(raw_object)# Grab the full MitoCarta2.0 call sheet as a dataframe.
    >>>mitodex,nonmitodex = mitoCartaPepOut(raw_object,dex=True) #Grab the mito/non-mito peptides for plotting by setting dex to True

    See Also
    --------
    omin.experiment.RawData
    omin.vis.plotByMito

    """
    peptides = raw_file.peptides
    proteins = raw_file.proteins
    carta = omin.mitoCartaCall.mitoProt(proteins)
    pepsel,prosel = omin.vLook(peptides,proteins,mods)
    mitocarta_pep = pepsel.merge(carta,on="Accession",how="left")
    mitocarta_pep.index = pepsel.index
    if dex:
        nonmitodex = mitocarta_pep.ix[mitocarta_pep.MitoCarta2_List != 1]
        mitodex = mitocarta_pep.ix[mitocarta_pep.MitoCarta2_List == 1]
        return mitodex,nonmitodex
    else:
        return mitocarta_pep

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
###INTERACTIVE SELECTION###
#----------------------------------------------------------------------------------------------------------------------

def treatmentSelect(peptide_abundance):
    """Allows the user to select which columns contain treatments.

    Parameters
    ----------
    peptide_abundance : DataFrame

    Returns
    -------
    select_list : list
    """
    print(pd.DataFrame([i.split(",") for i in peptide_abundance.columns]))
    col_num = int(input("Which column has treatment data?(Enter the number)"))
    select_set = set(pd.DataFrame([i.split(",") for i in peptide_abundance.columns]).ix[:, col_num])

    select_list = [re.sub(" ", "_", i.strip()) for i in select_set]
    #print(select_list)
    [print(n, i) for n, i in enumerate(select_list)]
    remove_num = int(input("Enter the number of any element that need to be removed."))
    select_list.remove(select_list[remove_num])
    return select_list
