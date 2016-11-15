import sys
from math import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import modulocator
from modulocator import modulocator
modulocator("Notebook")
import omin

class FileSelect(QDialog):

    def __init__(self, parent=None):

        super(FileSelect, self).__init__(parent)
        #BUTTONS
        self.peptides_path = QLineEdit("Peptides file")
        self.proteins_path = QLineEdit("Proteins file")
        self.peptide_browse_button = QPushButton("Browse")
        self.protein_browse_button = QPushButton("Browse")
        # self.next_button = QPushButton("Next")
        # self.back_button = QPushButton("Back")
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        #LABELS
        self.testlabel = QLabel("")
        # self.peptides_path.selectAll()

        #LAYOUT
        grid = QGridLayout()
        grid.addWidget(self.peptides_path,0,0)
        grid.addWidget(self.peptide_browse_button,0,1)
        grid.addWidget(self.proteins_path,1,0)
        grid.addWidget(self.protein_browse_button,1,1)
        # grid.addWidget(self.next_button,2,1)
        # grid.addWidget(self.back_button,2,0)
        grid.addWidget(self.testlabel,3,0)
        grid.addWidget(self.buttonBox,4,1)
        self.setLayout(grid)
        self.peptides_path.setFocus()
        #CONNECTIONS
        self.connect(self.peptide_browse_button, SIGNAL("clicked()"), lambda: self.setPath("peptides"))
        self.connect(self.protein_browse_button, SIGNAL("clicked()"), lambda: self.setPath("proteins"))
        self.connect(self.buttonBox,SIGNAL("accepted()"),self.accept)
        self.connect(self.buttonBox,SIGNAL("rejected()"),self.reject)
        #WINDOW TITLE
        self.setWindowTitle("Select Files")

    def setPath(self,tool):
        # path = QFileDialog.getExistingDirectory(self, "Select File")
        path = QFileDialog.getOpenFileName(self, "Select File")
        if path:
            self.__dict__[tool+"_path"].setText(QDir.toNativeSeparators(path))

    def testLabelChange(self):
        data = omin.RawData(self.peptides_path.text(),self.proteins_path.text()).peptides.shape
        self.testlabel.setText("Peptide rows:"+str(data[0]))

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    form = FileSelect()
    form.show()
    app.exec_()
