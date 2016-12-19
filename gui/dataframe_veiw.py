import numpy as np
import pandas as pd
from PyQt4.QtGui import *
from PyQt4.QtCore import *

from spyder.utils.qthelpers import (add_actions, create_action,
                                    keybinding, qapplication)


class PandasModel(QAbstractTableModel):
    """
    Populates a TableView with a pandas DataFrame.
    """
    def __init__(self, dataframe, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.df = dataframe
        self.data = np.array(self.df.values)
        self.cols = self.df.columns
        self.r, self.c = np.shape(self.data)

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                row_column = tuple([index.row(), index.column()])
                return self.data.item(row_column)
        return None

    def headerData(self, p_int, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.cols[p_int]
            elif orientation == Qt.Vertical:
                return p_int
        return None


class TableView(QTableView):
    """
    A simple implementation of QTableView.

    You can find all of the relevent documentation for this at the following
    link;
    http://pyqt.sourceforge.net/Docs/PyQt4/qtableview.html

    """
    def __init__(self, parent, model):
        QTableView.__init__(self, parent)
        self.setModel(model)
        # RIGHT-CLICK CONTEXT MENU
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.openMenu)

        # Create defaults for column headers.
        # Documentation for qheaderview can be found at the following link;
        # http://pyqt.sourceforge.net/Docs/PyQt4/qheaderview.html
        self.header = self.horizontalHeader()
        self.header.setDefaultSectionSize(100)
        self.header.setDefaultAlignment(Qt.AlignLeft)
        #HEADER CLICK SORT FUNCTION CONNECTION
        # self.header.setContextMenuPolicy(Qt.CustomContextMenu)
        # self.header.customContextMenuRequested.connect(self.openMenu)

    def openMenu(self, position):
        """Right-click context menu.
        """
        print(self.__dict__)
        menu = QMenu()
        menu.addAction("Edit")#self.tr("Edit"))
        menu.exec_(self.viewport().mapToGlobal(position))

    # def headerClicked(self, logicalIndex):
    #     """When a header is clicked corresponding column is sorted.
    #
    #     Parameters
    #     ----------
    #     logicalIndex : int
    #
    #     """
    #     menu = QMenu()
    #     menu.addAction("Edit")#self.tr("Edit"))
    #     menu.exec_(self.viewport().mapToGlobal(logicalIndex))
        # print(logicalIndex)
        # self.order = self.header.sortIndicatorOrder()
        # self.pdata = self.get_data_frame()
        # self.pdata.sort_values(by=self.pdata.columns[logicalIndex],
        #                        ascending=self.order,
        #                        inplace=True,
        #                        kind='mergesort')
        # self._tm = PandasModel(self.pdata)
        # self._tv.setModel(self._tm)
        # self._tv.update()



if __name__ == "__main__":
    from sys import argv, exit

    class TestWidget(QWidget):
        """
        A widget to test the classes above.
        """
        def __init__(self, parent=None):
            super(TestWidget, self).__init__(parent)
            # QWidget.__init__(self, parent)

            layout = QVBoxLayout(self)
            cdf = self.get_data_frame()
            self._tm = PandasModel(cdf)
            # self._tm.update(cdf)
            self.tv = TableView(self,PandasModel(cdf))
            # self.tv.setModel(self._tm)
            layout.addWidget(self.tv)

        def get_data_frame(self):
            file_path = "ExampleData\_E749_4154_010716_PeptideGroups.txt"
            dataframe = pd.read_csv(file_path,
                                    delimiter="\t",
                                    low_memory=False)
            return dataframe

    app = QApplication(argv)
    w = TestWidget()
    w.show()
    app.exec_()
