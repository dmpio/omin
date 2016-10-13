# -*- coding: utf-8 -*-
import pandas as pd
import omin

class Ligarep:
    """This class is deprecated.
    """
    def __init__(self,condlist,control):
        testlist = condlist.copy()
        testlist.remove(control)
        self.lfc = pd.concat([omin.tooLogFC(self.ko,self.wt,i) for i in condlist],axis=1)
        self.pval = pd.concat([omin.tooPvalr(self.ko,self.wt,i) for i in condlist],axis=1)
