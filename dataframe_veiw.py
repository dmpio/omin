import numpy as np
import pandas as pd
from PyQt4 import QtCore, QtGui

class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = np.array(data.values)
        self._cols = data.columns
        self.r, self.c = np.shape(self._data)

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return self._data[index.row(),index.column()]
        return None

    def headerData(self, p_int, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return self._cols[p_int]
            elif orientation == QtCore.Qt.Vertical:
                return p_int
        return None

class TableView(QtGui.QTableView):
    """
    A simple table to demonstrate the QComboBox delegate.
    """
    def __init__(self, *args, **kwargs):
        QtGui.QTableView.__init__(self, *args, **kwargs)


# if __name__=="__main__":
#     from sys import argv, exit
#
#     class Widget(QtGui.QWidget):
#         """
#         A simple test widget to contain and own the model and table.
#         """
#         def __init__(self, parent=None):
#             QtGui.QWidget.__init__(self, parent)
#
#             l=QtGui.QVBoxLayout(self)
#             cdf = self.get_data_frame()
#             self._tm=PandasModel(cdf)
#             # self._tm.update(cdf)
#             self._tv=TableView(self)
#             self._tv.setModel(self._tm)
#             l.addWidget(self._tv)
#
#         def get_data_frame(self):
#             df = pd.DataFrame({'Name':['a','b','c','d'],
#             'First':[2.3,5.4,3.1,7.7], 'Last':[23.4,11.2,65.3,88.8], 'Class':[1,1,2,1], 'Valid':[True, True, True, False]})
#             return df
#
#     a=QtGui.QApplication(argv)
#     w=Widget()
#     w.show()
#     w.raise_()
# exit(a.exec_())
