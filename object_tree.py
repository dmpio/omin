# 2016-10-20 23:34:21
# LOAD BOILER-PLATE CODE
import pandas as pd
import numpy as np

from modulocator import modulocator
modulocator("Notebook")
import omin

# CrAT Normalization to Pool with omin's Experiment class.
cratpepfile = "ExampleData/_E749_4154_010716_PeptideGroups.txt"
cratprotfile = "ExampleData/_E749_4154_010716_Proteins.txt"
modifications = ["Aceyl","Phospho"]
treatments = ["NonEx","Immediate Post","60 min post"]
crat_raw = omin.RawData(cratpepfile,cratprotfile)
# crat_exp = omin.Experiment(crat_raw,modifications)
crat_exp = omin.Experiment(crat_raw,modifications,treatments=treatments)

def myTypeCheck(mysterious_object):
    """Returns type of object as string.

    Parameters
    ----------
    mysterious_object : (:obj)
        Literally any type of object.

    Returns
    -------
    result_formatted : str
        The type of the mysterious_object as a string.
    """
    import re
    type_as_string = str(type(mysterious_object))
    #Define the search pattern. Here the the pattern is 'match all characters between single quotes'.
    search_pattern = "\\'.+\\'"
    result = re.search(search_pattern,type_as_string)
    result_formatted = result.group().strip("\\'")
    return result_formatted

def attForm(my_object):
    object_list = []
    if "omin" in myTypeCheck(my_object):
        attr_list = my_object.__dict__.keys()
        for attr in attr_list:
            if "omin" in myTypeCheck(my_object.__dict__[attr]):
                attr = [attr,attForm(my_object.__dict__[attr])]
                object_list.append(attr)
            else:
                object_list.append([attr,[]])
    return object_list

import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class Window(QWidget):

    def __init__(self):

        QWidget.__init__(self)
        self.treeView = QTreeView()
        self.model = QStandardItemModel()
        data = attForm(crat_exp)
        self.addItems(self.model, data)
        self.treeView.setModel(self.model)
        self.model.setHorizontalHeaderLabels([self.tr("Object")])
        layout = QVBoxLayout()
        layout.addWidget(self.treeView)
        self.setLayout(layout)

    def addItems(self, parent, elements):
        for text, children in elements:
            item = QStandardItem(text)
            parent.appendRow(item)
            if children:
                self.addItems(item, children)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
