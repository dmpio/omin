# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
from omin.utils import SelectionTools

# Load this directory
this_dir, _ = os.path.split(__file__)

# Add the selected file.
carta_file_path = this_dir + "\MitoCarta2.0_Mouse.xlsx"

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
    ent = pd.DataFrame(list_df, index=uni2ent.index)
    uni2ent = pd.concat([ent, uni2ent], axis=1)
    uni2ent.columns = ["MouseGeneID", "MGI", "Accession"]
    mitoprot = uni2ent.merge(cmitodf, on="MouseGeneID", how="left").dropna()
    return mitoprot
