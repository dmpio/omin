# -*- coding: utf-8 -*-
import pandas as pd
import omin

class Modgeno:
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
        control : str
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

class Nonmodgeno:
    def __init__(self,wt,ko):
        self.wt = wt
        self.ko = ko
    
class Geno:
    """
    Separates each group of modfied peptides into genotypes namely wild type
    (wt) and knockout(ko).
    """
    def __init__(self,abun,condlist,control,mitopep=None,modified=True):
        if modified == True:
            self.ko,self.wt = omin.sepCon(abun,"KO")
            wtpho,wtace = omin.sepCon(self.wt,"Phospho")
            kopho,koace = omin.sepCon(self.ko,"Phospho")
            self.ace = Modgeno(wtace,koace,condlist,control)
            self.ace.mitopep = mitopep
            self.pho = Modgeno(wtpho,kopho,condlist,control)
            self.pho.mitopep = mitopep
        if modified == False:
            self.ko,self.wt = omin.sepCon(abun,"KO")
            wtpho,wtace = omin.sepCon(self.wt,"Phospho")
            kopho,koace = omin.sepCon(self.ko,"Phospho")
            self.ace = Nonmodgeno(wtace,koace)
            self.pho = Nonmodgeno(wtpho,kopho)
            
class Rep:
    def __init__(self,abun,condlist,control,mitopep=None,modified=True):
        if modified == True:
            rep1,rep2 = omin.sepCon(abun,"1..Replicate")
            self.rep1 = Geno(rep1,condlist,control,mitopep)
            self.rep2 = Geno(rep2,condlist,control,mitopep)
        if modified == False:
            rep1,rep2 = omin.sepCon(abun,"1..Replicate")
            self.rep1 = Geno(rep1,condlist,control,modified)
            self.rep2 = Geno(rep2,condlist,control,modified)
        
class Replicate:
    def __init__(self,pep,prot):
        self.pep = pep
        self.prot = prot
        #separate peptides based on modifications
        self.pepacepho,self.peppho,self.pepace = omin.modSel(pep,"Phospho","Acetyl")
        #Mitocarta2.0 Calls this part could probably be done better
        dex = self.pepacepho["Master Protein Accessions"].dropna()
        self.mpa = pd.DataFrame([i.split(";")[0] for i in dex],index=dex.index)
        self.mpa.columns = ["Accession"]
        self.mpa = self.mpa.merge(prot[["Accession","Gene ID"]],on="Accession",how="left")
        self.mpa.index = dex.index
        self.mitoprot = omin.mitoCartaCall.mitoProt(prot)
        self.mitopep = self.mpa.merge(self.mitoprot,on="Accession",how='left')
        self.mitopep.index = self.mpa.index
        self.bdex = self.mitopep[self.mitopep.MitoCarta2_List==1.0]
        self.wdex = self.mitopep[self.mitopep.MitoCarta2_List!=1.0]
        #sepate all peptide abundance columns
        #FIXME phoacemod is probably redundant
        phoacemod = pep.Modifications.str.contains("Acetyl")^pep.Modifications.str.contains("Phospho")
        mod = pep[phoacemod].copy()
        nonmod = pep[~phoacemod].copy()
        self.modpepabun = omin.sep(mod,"Abundance:")
        self.nonmodpepabun = omin.sep(nonmod,"Abundance:")
        
    def diffGeno(self,condlist,control):
        self.mod = Rep(self.modpepabun,condlist,control,self.mitopep)
        self.nonmod = Rep(self.nonmodpepabun,condlist,control,modified=False)
