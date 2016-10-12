# -*- coding: utf-8 -*-
import pandas as pd
import omin
class Ligarep:
    def __init__(self,condlist,control):
        testlist = condlist.copy()
        testlist.remove(control)

        self.lfc = pd.concat([omin.tooLogFC(self.ko,self.wt,i) for i in condlist],axis=1)
        self.pval = pd.concat([omin.tooPvalr(self.ko,self.wt,i) for i in condlist],axis=1)
        
        #self.epval_ko = pd.concat([omin.pvalr(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
        #self.elfc_ko = pd.concat([omin.logFC(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
        #self.epval_wt = pd.concat([omin.pvalr(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)
        #self.elfc_wt = pd.concat([omin.logFC(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)

#class Ligarep:
#    def __init__(self,ko,wt,condlist,control):
#        testlist = condlist.copy()
#        testlist.remove(control)
#        
#        self.ko = ko
#        self.wt = wt
#        
#        self.lfc = pd.concat([omin.tooLogFC(self.ko,self.wt,i) for i in condlist],axis=1)
#        self.pval = pd.concat([omin.tooPvalr(self.ko,self.wt,i) for i in condlist],axis=1)
#        
#        self.epval_ko = pd.concat([omin.pvalr(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
#        self.elfc_ko = pd.concat([omin.logFC(omin.overPooler(self.ko),i,control) for i in testlist],axis=1)
#        self.epval_wt = pd.concat([omin.pvalr(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)
#        self.elfc_wt = pd.concat([omin.logFC(omin.overPooler(self.wt),i,control) for i in testlist],axis=1)
