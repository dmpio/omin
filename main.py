import os
import platform
import stat
import sys

import modulocator
from modulocator import modulocator
modulocator("Notebook")
import omin

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import file_select


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.filename = None
        self.peptides_path = None
        self.proteins_path = None
        self.omin_object = None
        #Set window TITLE
        self.setWindowTitle("omin-omics analytics")

        treeLabel = QLabel("Tre&e")
        self.treeWidget = QTreeWidget()
        treeLabel.setBuddy(self.treeWidget)
        self.setCentralWidget(self.treeWidget)
        #TABLE WIDGET
        tableDockWidget = QDockWidget("Table",self)
        tableDockWidget.setObjectName("TableDockWidget")
        tableDockWidget.setAllowedAreas(Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)
        self.tableWidget = QTableWidget()
        tableDockWidget.setWidget(self.tableWidget)
        self.addDockWidget(Qt.RightDockWidgetArea,tableDockWidget)
        #LOG WIDGET
        logDockWidget = QDockWidget("Log", self)
        logDockWidget.setObjectName("LogDockWidget")
        logDockWidget.setAllowedAreas(Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea|Qt.TopDockWidgetArea|Qt.BottomDockWidgetArea)
        self.listWidget = QListWidget()
        logDockWidget.setWidget(self.listWidget)
        self.addDockWidget(Qt.BottomDockWidgetArea, logDockWidget)

        fileNewAction = self.createAction("&New...",self.fileNew,QKeySequence.New, "filenew", "Create an image file")

        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenuActions = (None, fileNewAction)
        self.connect(self.fileMenu, SIGNAL("aboutToShow()"), self.updateFileMenu)
        self.populateTree()
        self.populateTable()
        self.connect(self.treeWidget,SIGNAL('selectionChanged(QItemSelection, QItemSelection)'), self.test)
        # self.connect(self.treeWidget,SIGNAL('clicked()'),self.updateSatus("hey"))
        # self.connect(self.treeWidget.currentItem(),SIGNAL('selectionChanged()'),self.updateSatus("hey"))
        # self.fileMenuActions.addActions(self.fileMenu,self.fileMenuActions)

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
            self.omin_object = omin.RawData(self.peptides_path,self.proteins_path)
            self.updateSatus("omin project loaded.")
            # self.imageLabel.setText(self.peptides_path)

    def addActions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def updateSatus(self,message):
        self.listWidget.addItem(message)
        self.populateTree()

    def updateFileMenu(self):
        self.fileMenu.clear()
        self.addActions(self.fileMenu, self.fileMenuActions)

    def populateTree(self):
        self.treeWidget.clear()
        self.treeWidget.setHeaderLabels(["project"])
        if self.omin_object is not None:
            # return
            att_list = self.omin_object.__dict__.keys()
            for att in att_list:
                QTreeWidgetItem(self.treeWidget, [att])
        else:
            att_list = ["None Open"]
            for att in att_list:
                QTreeWidgetItem(self.treeWidget, [att])

    def populateTable(self):
        selected = None
        self.tableWidget.clear()

    def treeSelect(self):
        term = self.treeWidget.currentItem()
        self.updateSatus(term)

    @pyqtSlot("QItemSelection, QItemSelection")
    def test(self, selected, deselected):
        print("hello!")
        print(selected)
        print(deselected)



def main():
    app = QApplication(sys.argv)
    # app.setWindowIcon(QIcon(":/icon.png"))
    main = MainWindow()
    main.show()
    app.exec_()

main()
