# -*- coding: utf-8 -*-
from . import sep
import pandas as pd
import numpy as np
import os,sys
nb_dir = os.path.split(os.getcwd())[0]
if nb_dir not in sys.path:
    sys.path.append(nb_dir) 

def mitodf():
    mitodf = pd.read_excel(nb_dir+'\\'+"omin\Databases\MitoCarta2.0_Mouse.xlsx")
    return mitodf

mitodf = mitodf()
mgi = mitodf.MouseGeneID.copy()
#make copy of mitodf
cmitodf = mitodf.copy()
#Remove useless columns
cmitodf.drop(cmitodf.columns[1:5], axis=1, inplace=True)

def mitoProt(protdata):
    uni = protdata.Accession.copy()
    uni = uni.as_matrix()
    ent = sep(protdata,"Entrez")
    ent = ent.as_matrix().T[0]
    uni2ent = pd.DataFrame([ent,uni]).T
    uni2ent = uni2ent.dropna()
    e = pd.DataFrame([np.int64(i.split(";")[0]) for i in uni2ent[uni2ent.columns[0]]],
                      index=uni2ent.index)
    uni2ent = pd.concat([e,uni2ent],axis=1)
    uni2ent.columns = ["MouseGeneID","MGI","Accession"]
    mitoprot = uni2ent.merge(cmitodf,on="MouseGeneID",how="left").dropna()
    return mitoprot

#dex = pepacepho["Master Protein Accessions"].dropna()