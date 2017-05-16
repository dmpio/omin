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

import os
import re
import pandas as pd
import numpy as np
from ..utils import SelectionTools

# Load this directory
this_dir, _ = os.path.split(__file__)

# Add the selected file.
mitocarta_local = "/MitoCarta2.0_Mouse.xlsx"

if os.name == "nt":
    mitocarta_local = mitocarta_local.replace("/", "\\")

carta_file_path = this_dir + mitocarta_local

# load MitoCarta2
mitodf = pd.read_excel(carta_file_path)

# Load MouseGeneID as a Series
mgi = mitodf.MouseGeneID.copy()

# Make copy of mitodf
cmitodf = mitodf.copy()

# Remove useless columns
cmitodf.drop(cmitodf.columns[1:5], axis=1, inplace=True)


def mitoProt(protdata):
    """Takes protein DataFrame and returns a the mitocarta2 calls

    Parameters
    ----------
    protdata : DataFrame
        Protien DataFrame

    Returns

    -------
    mitoprot : DataFrame
        Your proteins with the relevent mitocarta2 data.

    """
    uni = protdata.Accession.copy()
    uni = uni.as_matrix()
    ent = SelectionTools.sep(protdata, "Entrez")
    ent = ent.as_matrix().T[0]
    uni2ent = pd.DataFrame([ent, uni]).T
    uni2ent = uni2ent.dropna()
    # create list for DataFrame
    list_df = [np.int64(i.split(";")[0]) for i in uni2ent[uni2ent.columns[0]]]
    # Try to replace the following with ent.to_frame()
    ent = pd.DataFrame(list_df, index=uni2ent.index)
    uni2ent = pd.concat([ent, uni2ent], axis=1)
    uni2ent.columns = ["MouseGeneID", "MGI", "Accession"]
    mitoprot = uni2ent.merge(cmitodf, on="MouseGeneID", how="left").dropna()
    return mitoprot


def mitoCartaPepOut(obj=None, mods=None, dex=False):
    """
    Parameters
    ----------
    obj : (:obj)
        Any object that has a attribute that is named peptides and proteins.
    mods : list
        Defaults to ["Acetyl","Phospho"].
    dex : bool
        Defaults to False. When False output is mitocarta_pep if True
        output is a tuple containing mitodex and nonmitodex

    Returns
    -------
    (mitodex, nonmitodex) : tuple(DataFrame,DataFrame)
    mitocarta_pep : Dataframe

    Examples
    --------
    Grab all MitoCarta 2.0 calls for all peptides as a dataframe.
    >>>mitocarta_pep = mitoCartaPepOut(raw_object)

    Grab the mito/non-mito peptides for plotting by setting dex to True
    >>>mitodex,nonmitodex = mitoCartaPepOut(raw_object,dex=True)

    """
    mods = mods or ["Acetyl", "Phospho"]

    try:
        rx = re.compile("[Pp]eptides")
        peptides = list(filter(rx.findall, obj.__dict__.keys()))[0]
        rx = re.compile("[Pp]roteins")
        proteins = list(filter(rx.findall, obj.__dict__.keys()))[0]

    except Exception:
        peptides = "peptides"
        proteins = "proteins"

    peptides = obj.__dict__[peptides]
    proteins = obj.__dict__[proteins]

    carta = mitoProt(proteins)
    pepsel, prosel = SelectionTools.vLook(peptides, proteins, mods)
    mitocarta_pep = pepsel.merge(carta, on="Accession", how="left")
    mitocarta_pep.index = pepsel.index
    if dex:
        nonmitodex = mitocarta_pep.ix[mitocarta_pep.MitoCarta2_List != 1]
        mitodex = mitocarta_pep.ix[mitocarta_pep.MitoCarta2_List == 1]
        return mitodex, nonmitodex
    else:
        return mitocarta_pep
