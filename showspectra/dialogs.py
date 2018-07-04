import os
from PyQt5.QtCore import pyqtSignal,QObject
from PyQt5.QtWidgets import QDialog, QPushButton, QGroupBox, QHBoxLayout, QVBoxLayout, QGridLayout, QRadioButton, QButtonGroup, QFileDialog


class selectTelescope(QDialog):
    """ Select telescope """

    def __init__(self, k=0, parent=None):
        super().__init__(parent)

        self.k = k
        self.setupUI()

    def setupUI(self):

        self.telescope = self.createGroup('Telescope',['WIYN','VIMOS','SDSS'])

        hgroup = QGroupBox()
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)
        
        grid = QGridLayout()
        grid.addWidget(self.telescope,0,0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        self.setWindowTitle('Select a telescope to start')
        self.resize(400,300)

    def createGroup(self, title, items, default=0):    
        """ creates a group of radio buttons  """
        group = QGroupBox(title)
        group.buttons = QButtonGroup()
        vbox = QVBoxLayout()
        buttons = []
        i = 0
        for item in items:
            buttons.append(QRadioButton(item))
            group.buttons.addButton(buttons[-1], i)
            vbox.addWidget(buttons[-1])
            i += 1
        vbox.addStretch(1)
        # Set 1st option as default
        buttons[default].setChecked(True)
        group.setLayout(vbox)
        return group


    def OK(self):
        self.done(1)

    def save(self):
        telescope  = self.telescope.buttons.checkedButton().text()
        return telescope
            
    def Cancel(self):
        self.done(0)

class selectFiles(QDialog):
    """ Selection of files (Flux, Error, and Sky) """
    def __init__(self, parent=None):
        super().__init__(parent)

        self.flux = self.err = self.sky = None
        
        self.setupUI()

    def setupUI(self):


        hgroup = QGroupBox()
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)

        vgroup = QGroupBox()
        vbox = QVBoxLayout()
        self.b1 = QPushButton("Flux")
        #self.b1.setCheckable(True)
        #self.b1.toggle()
        self.b1.clicked.connect(self.selectFlux)
        self.b2 = QPushButton("Err")
        #self.b2.setCheckable(True)
        #self.b2.toggle()
        self.b2.clicked.connect(self.selectErr)
        self.b3 = QPushButton("Sky")
        #self.b3.setCheckable(True)
        #self.b3.toggle()
        self.b3.clicked.connect(self.selectSky)
        vbox.addWidget(self.b1)
        vbox.addWidget(self.b2)
        vbox.addWidget(self.b3)
        vgroup.setLayout(vbox)

        grid = QGridLayout()
        grid.addWidget(vgroup,0,0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        
        self.setWindowTitle("Select spectral files")
        self.resize(400,300)


    def selectFlux(self):
        self.flux = self.selectFile('Flux')
        if self.flux is not None:
            head, tail = os.path.split(self.flux)
            self.b1.setText('Flux: '+tail)

    def selectErr(self):
        self.err = self.selectFile('Error')
        if self.err is not None:
            head, tail = os.path.split(self.flux)
            self.b2.setText('Err: '+tail)
            
    def selectSky(self):
        self.sky = self.selectFile('Sky')
        if self.sky is not None:
            head, tail = os.path.split(self.flux)
            self.b3.setText('Sky: '+tail)
        
    def selectFile(self,label):
        fd = QFileDialog()
        fd.setWindowTitle('Open '+label+' File')
        fd.setLabelText(QFileDialog.Accept, "Select "+label)
        fd.setNameFilters(["Fits Files (*.fits, *.fit)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            fileNames = fd.selectedFiles()
            return fileNames[0]
        else:
            return None
            print("No file selected")

    def OK(self):
        if (self.flux is not None) & (self.err is not None ) & (self.sky is not None):
            self.done(1)
        else:
            print('Please, complete all the fields')

    def save(self):
        return self.flux, self.err, self.sky
            
    def Cancel(self):
        self.done(0)