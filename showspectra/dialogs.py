import os
from PyQt5.QtWidgets import (QDialog, QPushButton, QGroupBox,
                             QHBoxLayout, QVBoxLayout, QGridLayout,
                             QRadioButton, QButtonGroup, QFileDialog,
                             QListWidget, QListWidgetItem)
from PyQt5.QtCore import QSize
import numpy as np


class selectTelescope(QDialog):
    """Select telescope."""

    def __init__(self, k=0, parent=None):
        super().__init__(parent)

        self.k = k
        self.setupUI()

    def setupUI(self):

        self.telescope = self.createGroup('Telescope', ['WIYN', 'VIMOS', 'SDSS','MUSE'])

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
        grid.addWidget(self.telescope, 0, 0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        self.setWindowTitle('Select a telescope to start')
        self.resize(400, 300)

    def createGroup(self, title, items, default=0):
        """Creates a group of radio buttons."""
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
        telescope = self.telescope.buttons.checkedButton().text()
        return telescope

    def Cancel(self):
        self.done(0)


class selectFiles(QDialog):
    """ Selection of files (Flux, Error, and Sky) """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.flux = self.err = self.sky = None
        self.setupUI()
        # self.selectFlux()

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
        self.b1.clicked.connect(self.selectFlux)
        self.b2 = QPushButton("Err")
        self.b2.clicked.connect(self.selectErr)
        self.b3 = QPushButton("Sky")
        self.b3.clicked.connect(self.selectSky)
        vbox.addWidget(self.b1)
        vbox.addWidget(self.b2)
        vbox.addWidget(self.b3)
        vgroup.setLayout(vbox)
        grid = QGridLayout()
        grid.addWidget(vgroup, 0, 0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        self.setWindowTitle("Select spectral files")
        self.resize(400, 300)

    def selectFlux(self):
        self.flux = self.selectFile('Flux')
        if self.flux is not None:
            path, file = os.path.split(self.flux)
            self.b1.setText('Flux: ' + file)
            # Try to guess error and sky files from the same directory
            from glob import glob as gb
            import re
            
            match1 = re.search('SDSS', file)
            match2 = re.search('MUSE', file)
            if match1 or match2:
                self.b2.setText('Err: '+file)
                self.b3.setText('Sky: '+file)
                self.err = file
                self.sky = file
            else:
                # check if pattern is found to guess the two other files             
                print('path is ', path)
                files = gb(path + '/*.fit*')
                print('Files are ', files)
                r = re.compile(".*err.*\.fit.*", flags=re.IGNORECASE)
                rErr = list(filter(r.match, files))
                if rErr is not None:
                    self.err = rErr[0]
                    head, tail = os.path.split(self.err)
                    self.b2.setText('Err: ' + tail)
                r = re.compile(".*sky.*\.fit.*", flags=re.IGNORECASE)
                rSky = list(filter(r.match, files))
                if rSky is not None:
                    self.sky = rSky[0]
                    head, tail = os.path.split(self.sky)
                    self.b3.setText('Sky: ' + tail)

    def selectErr(self):
        self.err = self.selectFile('Error')
        if self.err is not None:
            head, tail = os.path.split(self.err)
            self.b2.setText('Err: ' + tail)

    def selectSky(self):
        self.sky = self.selectFile('Sky')
        if self.sky is not None:
            head, tail = os.path.split(self.sky)
            self.b3.setText('Sky: ' + tail)

    def selectFile(self, label):
        fd = QFileDialog()
        fd.setWindowTitle('Open ' + label + ' File')
        fd.setLabelText(QFileDialog.Accept, "Select " + label)
        fd.setNameFilters(["Fits Files (*.fit* *.FIT*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            fileNames = fd.selectedFiles()
            return fileNames[0]
        else:
            return None
            print("No file selected")

    def OK(self):
        if (self.flux is not None) & (self.err is not None) & (self.sky is not None):
            self.done(1)
        else:
            print('Please, complete all the fields')

    def save(self):
        return self.flux, self.err, self.sky

    def Cancel(self):
        self.done(0)


class selectRedshift(QDialog):
    """Selection of optimal redshift/template from a list of cross-correlations with templates."""

    def __init__(self, z, sz, snr, temps, parent=None):
        super().__init__()
        self.z = z
        self.sz = sz
        self.snr = snr
        self.temps = temps
        self.zlist = ["z = %.5f +/- %.5f (snr %.1f) %s" % (z, sz, snr, t) for (z, sz, snr, t)
                      in zip(self.z, self.sz, self.snr, self.temps)]
        self.setupUI()

    def setupUI(self):
        self.setWindowTitle("Select best template redshift")
        layout = QVBoxLayout()
        hbox = QHBoxLayout()
        hgroup = QGroupBox()
        self.list = QListWidget(self)
        for z in self.zlist:
            item = QListWidgetItem(self.list)
            item.setText(z)
        self.list.setMaximumSize(QSize(400, 100))
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1)
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)
        # Layout
        layout.addWidget(self.list)
        layout.addWidget(hgroup)
        self.setLayout(layout)

    def OK(self):
        self.done(1)

    def Cancel(self):
        self.done(0)


class selectLine(QDialog):
    """Visual selection of line in the spectrum."""

    def __init__(self, w, lines, parent=None):
        super().__init__()
        self.names = []
        self.waves = []
        for line in lines.copy():
            lin = lines[line]
            self.names.append(lin[0])
            self.waves.append(lin[1])
        self.names = np.array(self.names)
        self.waves = np.array(self.waves)
        self.z = (w - self.waves) / self.waves
        self.zlist = ["%s %.5f %.4f" % (n, w, z) for (n, w, z)
                      in zip(self.names, self.waves, self.z)]
        self.setupUI()

    def setupUI(self):
        self.setWindowTitle("Select best line")
        layout = QVBoxLayout()
        hbox = QHBoxLayout()
        hgroup = QGroupBox()
        self.list = QListWidget(self)
        for z in self.zlist:
            item = QListWidgetItem(self.list)
            item.setText(z)
        self.list.setMaximumSize(QSize(400, 100))
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1)
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)
        # Layout
        layout.addWidget(self.list)
        layout.addWidget(hgroup)
        self.setLayout(layout)

    def OK(self):
        self.done(1)

    def Cancel(self):
        self.done(0)


class guessParams(QDialog):
    """ Dialog window to define guess parameters of the continuum and lines fit """

    def __init__(self, parent=None):
        super().__init__()
        self.setupUI()

    def setupUI(self):

        self.continuum = self.createGroup('Continuum', ['Constant', 'Slope'], default=0)
        self.emission = self.createGroup('Emission lines', ['0', '1', '2', '3'], default=1)
        self.absorption = self.createGroup('Absorption lines', ['0', '1', '2', '3'], default=0)

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
        grid.addWidget(self.continuum, 0, 0)
        grid.addWidget(self.emission, 1, 0)
        grid.addWidget(self.absorption, 2, 0)
        grid.addWidget(hgroup, 3, 0)
        self.setLayout(grid)
        self.setWindowTitle('Guess parameters')
        self.resize(400, 300)

    def createGroup(self, title, items, default=0):
        """Creates a group of radio buttons."""
        group = QGroupBox(title)
        group.buttons = QButtonGroup()
        hbox = QHBoxLayout()
        buttons = []
        i = 0
        for item in items:
            buttons.append(QRadioButton(item))
            group.buttons.addButton(buttons[-1], i)
            hbox.addWidget(buttons[-1])
            i += 1
        hbox.addStretch(1)
        # Set 1st option as default
        buttons[default].setChecked(True)
        group.setLayout(hbox)
        return group

    def OK(self):
        self.done(1)

    def save(self):
        continuum = self.continuum.buttons.checkedButton().text()
        emission = self.emission.buttons.checkedButton().text()
        absorption = self.absorption.buttons.checkedButton().text()
        return continuum, emission, absorption

    def Cancel(self):
        self.done(0)
