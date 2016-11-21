import numpy as np
import pandas as pd
from PyQt4.QtGui import *
from PyQt4.QtCore import *

class PandasModel(QAbstractTableModel):
    """
    Populates a TableView with a pandas DataFrame.
    """
    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = np.array(data.values)
        self._cols = data.columns
        self.r, self.c = np.shape(self._data)

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                row_column = tuple([index.row(),index.column()])
                return self._data.item(row_column)
        return None

    def headerData(self, p_int, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._cols[p_int]
            elif orientation == Qt.Vertical:
                return p_int
        return None

class TableView(QTableView):
    """
    A simple implementation of QTableView.
    """
    def __init__(self, *args, **kwargs):
        QTableView.__init__(self, *args, **kwargs)


if __name__=="__main__":
    from sys import argv, exit
    class TestWidget(QWidget):
        """
        A widget to test the classes above.
        """
        def __init__(self, parent=None):
            QWidget.__init__(self, parent)

            layout = QVBoxLayout(self)
            cdf = self.get_data_frame()
            self._tm=PandasModel(cdf)
            # self._tm.update(cdf)
            self._tv=TableView(self)
            self._tv.setModel(self._tm)
            layout.addWidget(self._tv)

        def get_data_frame(self):
            file_path = "ExampleData\_E749_4154_010716_PeptideGroups.txt"
            dataframe = pd.read_csv(file_path,delimiter="\t",low_memory=False)
            return dataframe

    app=QApplication(argv)
    w=TestWidget()
    w.show()
    app.exec_()
