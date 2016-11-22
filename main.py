import os
import platform
import stat
import sys
import pandas as pd
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import file_select
from dataframe_veiw import *

import modulocator
from modulocator import modulocator
modulocator("Notebook")
import omin


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
    # Define Pattern; 'match all characters between single quotes'.
    search_pattern = "\\'.+\\'"
    result = re.search(search_pattern, type_as_string)
    result_formatted = result.group().strip("\\'")
    return result_formatted


def attForm(my_object):
    object_list = []
    if "omin" in myTypeCheck(my_object):
        attr_list = my_object.__dict__.keys()
        for attr in attr_list:
            if "omin" in myTypeCheck(my_object.__dict__[attr]):
                attr = [attr, attForm(my_object.__dict__[attr])]
                object_list.append(attr)
            else:
                object_list.append([attr, [my_object.__dict__[attr]]])
    return object_list


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.filename = None
        self.peptides_path = None
        self.proteins_path = None
        self.omin_object = None
        # Set window TITLE
        self.setWindowTitle("omin-omics analytics")
        # TREE WIDGET
        self.treeView = QTreeView()
        self.model = QStandardItemModel()
        self.setCentralWidget(self.treeView)
        self.treeView.setModel(self.model)
        self.model.setHorizontalHeaderLabels([self.tr("Object")])
        # TABLE WIDGET
        # AKA DATAFRAME WIDGET
        tableDockWidget = QDockWidget("Table", self)
        tableDockWidget.setObjectName("TableDockWidget")
        tableDockWidget.setAllowedAreas(
            Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.tv = TableView()
        tableDockWidget.setWidget(self.tv)
        self.addDockWidget(Qt.RightDockWidgetArea, tableDockWidget)
        # LOG WIDGET
        logDockWidget = QDockWidget("Log", self)
        logDockWidget.setObjectName("LogDockWidget")
        # Set allowed areas for LogDockWidget
        logDockWidget.setAllowedAreas(
            Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea |
            Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
        # Create listWidget for logDockWidget
        self.listWidget = QListWidget()
        logDockWidget.setWidget(self.listWidget)
        self.addDockWidget(Qt.BottomDockWidgetArea, logDockWidget)
        # FILE MENU
        fileNewAction = self.createAction(
            "&New...", self.fileNew, QKeySequence.New,
            "filenew", "Create an image file")
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenuActions = (None, fileNewAction)
        self.connect(self.fileMenu, SIGNAL(
            "aboutToShow()"), self.updateFileMenu)
        # LOAD DATA
        self.populateTree()
        self.treeView.clicked.connect(self.on_tree_item_clicked)
        # self.populateTable()

    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=False, signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/{}.png".format(icon)))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def fileNew(self):
        dialog = file_select.FileSelect(self)
        dialog.show()
        if dialog.exec_():
            self.peptides_path = dialog.peptides_path.text()
            self.proteins_path = dialog.proteins_path.text()
            # FIXME : Make another dialog box that allows the user to select
            # the modifications
            modifications = ["Aceyl", "Phospho"]
            treatments = ["NonEx", "Immediate Post", "60 min post"]
            omin_raw = omin.RawData(self.peptides_path, self.proteins_path)
            self.omin_object = omin.Experiment(
                omin_raw, modifications, treatments=treatments)
            self.updateSatus("omin project loaded.")
            self.populateTree()
            # self.imageLabel.setText(self.peptides_path)

    def addActions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def updateSatus(self, message):
        self.listWidget.addItem(message)
        # self.populateTree()

    def updateFileMenu(self):
        self.fileMenu.clear()
        self.addActions(self.fileMenu, self.fileMenuActions)

    def populateTree(self):
        self.treeView.reset()
        # self.treeView.setHeaderLabels(["project"])
        if self.omin_object is not None:
            # self.treeView.
            self.omin_object_data = attForm(self.omin_object)
            self.addItems(self.model, self.omin_object_data)
        else:
            pass

    def populateTable(self, dataframe=None):
        if dataframe is not None:
            self.tv.reset()
            self.pm = PandasModel(dataframe)
            self.tv.setModel(self.pm)
            self.tv.update()
        else:
            pass

    @pyqtSlot(QModelIndex)
    def on_tree_item_clicked(self, index):
        selection = self.model.itemFromIndex(index)
        self.updateSatus(selection.text())
        # print(selection.data())
        self.populateTable(selection.data())

    def addItems(self, parent, elements):
        for text, children in elements:
            item = QStandardItem()
            item.setText(text)
            # item.setData("Hello!")
            parent.appendRow(item)
            if children:
                if "pandas" not in myTypeCheck(children[0]):
                    self.addItems(item, children)
                else:
                    item.setData(children[0])


def main():
    app = QApplication(sys.argv)
    # app.setWindowIcon(QIcon(":/icon.png"))
    main = MainWindow()
    main.show()
    app.exec_()

main()
