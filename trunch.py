import omin
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

class modgeno:
    """Separates each group of modfied peptides into genotypes namely wild type (wt) and knockout(ko).

    Attributes
    ----------
    wt : DataFrame
    ko : DataFrame
    lfc : DataFrame
    pval : DataFrame
    epval_ko : DataFrame
    elfc_ko : DataFrame
    epval_wt : DataFrame
    elfc_wt : DataFrame

    """
    def __init__(self,wt,ko,condlist,control):
        """
        Parameters
        ----------
        wt : DataFrame
        ko : DataFrame
        condlist : list
            List of conditions.
        control : str
            A the control group amoung the condition list.

        """
        testlist = condlist.copy()
        testlist.remove(control)
        self.wt = wt
        self.ko = ko
        self.lfc = pd.concat([omin.tooLogFC(self.ko,self.wt,i) for i in condlist],axis=1)
        self.pval = pd.concat([omin.tooPvalr(self.ko,self.wt,i) for i in condlist],axis=1)
        self.epval_ko = pd.concat([omin.pvalr(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
        self.elfc_ko = pd.concat([omin.logFC(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
        self.epval_wt = pd.concat([omin.pvalr(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)
        self.elfc_wt = pd.concat([omin.logFC(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)

class Trunch:
    """Takes peptide file DataFrame and protein file DataFrame.

    Attributes
    ----------
    peptides : DataFrame
    proteins : DataFrame
    pepace : DataFrame
    peppho : DataFrame
    pepace : DataFrame
    mpa : DataFrame
    mitoprot : DataFrame
    mitopep : DataFrame
    bdex : DataFrame
    wdex : DataFrame
    allpepabun : DataFrame
    modpepabun : DataFrame
    ko : DataFrame
    wt : DataFrame
    """
    def __init__(self,pep,prot):
        """
        Notes
        -----
            Need handling for replicates and way for user to select the conditions and control from abundance columns.

        Parameters
        ----------
        pep : DataFrame
        prot : DataFrame

        """
        self.peptides = pep
        self.proteins = prot
        #separate peptides based on modifications
        self.pepacepho,self.peppho,self.pepace = omin.modSel(pep,"Phospho","Acetyl")
        #Mitocarta2.0 Calls this part could probably be done better
        dex = self.pepacepho["Master Protein Accessions"].dropna()
        self.mpa = pd.DataFrame([i.split(";")[0] for i in dex],index=dex.index)
        self.mpa.columns = ["Accession"]
        #self.mpa = self.mpa.merge(prot[["Accession","Gene ID"]],on="Accession",how="left")
        self.mpa = self.mpa.merge(prot[["Accession", "Gene ID","Description"]], on="Accession", how="left")
        self.mpa.index = dex.index
        #self.mitoprot = omin.mitoCartaCall.mitoProt(prot)
        self.mitoprot = omin.mitoProt(prot)
        self.mitopep = self.mpa.merge(self.mitoprot,on="Accession",how='left')
        self.mitopep.index = self.mpa.index
        self.bdex = self.mitopep[self.mitopep.MitoCarta2_List==1.0]
        self.wdex = self.mitopep[self.mitopep.MitoCarta2_List!=1.0]
        #sepate all peptide abundance columns
        self.allpepabun = omin.sep(pep,"Abundance:")
        #from all peptide abundance columns separate all modified proteins
        self.modpepabun = self.allpepabun.iloc[self.pepacepho.index].copy()
        self.ko,self.wt = omin.sepCon(self.modpepabun,"KO")#FIXME: KO cannot be hardcoded

    def diffGeno(self,condlist,control):
        """Separate modification columns

        Parameters
        ----------
        condlist : list
        control : str
        
        """
        wtpho,wtace = omin.sepCon(self.wt,"Phospho")
        kopho,koace = omin.sepCon(self.ko,"Phospho")
        self.ace = modgeno(wtace,koace,condlist,control)
        self.ace.mitopep = self.mitopep
        self.pho = modgeno(wtpho,kopho,condlist,control)
        self.pho.mitopep = self.mitopep
