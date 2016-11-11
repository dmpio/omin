import sys
from math import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class FileSelect(QDialog):

    def __init__(self, parent=None):
        super(FileSelect, self).__init__(parent)
        #Create QLineEdit boxes for peptides_path and protiens_path
        self.peptides_path = QLineEdit("Peptides file")
        self.protiens_path = QLineEdit("Protiens file")

        self.peptide_browse_button = QPushButton("Browse")
        self.protien_browse_button = QPushButton("Browse")
        self.next_button = QPushButton("Next")
        self.back_button = QPushButton("Back")
        self.testlabel = QLabel("")
        # self.peptides_path.selectAll()
        grid = QGridLayout()
        grid.addWidget(self.peptides_path,0,0)
        grid.addWidget(self.peptide_browse_button,0,1)
        grid.addWidget(self.protiens_path,1,0)
        grid.addWidget(self.protien_browse_button,1,1)
        grid.addWidget(self.next_button,2,1)
        grid.addWidget(self.back_button,2,0)
        grid.addWidget(self.testlabel,3,0)

        self.setLayout(grid)
        self.peptides_path.setFocus()
        # self.connect(self.peptide_browse_button, SIGNAL("returnPressed()"),self.updateUi)
        self.connect(self.peptide_browse_button, SIGNAL("clicked()"), lambda: self.setPath("peptides"))
        self.connect(self.protien_browse_button,SIGNAL("clicked()"),lambda: self.setPath("protiens"))
        self.connect(self.next_button,SIGNAL("clicked()"),self.testLabelChange)
        self.setWindowTitle("Select Files")

    def setPath(self,tool):
        path = QFileDialog.getExistingDirectory(self, "Select File")# self.__dict__[tool+"_path"])
        if path:
            self.__dict__[tool+"_path"].setText(QDir.toNativeSeparators(path))

    def testLabelChange(self):
        self.testlabel.setText(self.peptides_path.text())

app = QApplication(sys.argv)
form = FileSelect()
form.show()
app.exec_()
