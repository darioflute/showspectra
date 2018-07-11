#!/usr/bin/env python

import sys,os
import numpy as np

from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QTabBar,QHBoxLayout,
                             QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog,  QTableView, QComboBox, QAbstractItemView,
                             QMessageBox, QInputDialog, QDialog, QLabel)
from PyQt5.QtGui import QIcon, QImage, QStandardItem, QStandardItemModel, QPixmap, QMovie
from PyQt5.QtCore import Qt, QSize, QTimer, QThread, QObject, pyqtSignal

from showspectra.graphics import SpectrumCanvas, NavigationToolbar
from showspectra.dialogs  import selectTelescope, selectFiles
from showspectra.inout    import saveAnalysis, recoverAnalysis


class GUI (QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Show Spectra'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480

        # initialize galaxies
        self.galaxies = None

        # Ordinate closing
        self.setAttribute(Qt.WA_DeleteOnClose)

        # Path of the package
        self.path0, file0 = os.path.split(__file__)
        self.setWindowIcon(QIcon(QPixmap(self.path0+'/icons/showspectra.png')))

        # Style
        # Define style
        with open(self.path0+'/stylesheet.css',"r") as fh:
            self.setStyleSheet(fh.read())

        self.initUI()

        

    def initUI(self):
        """ User interface """

        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top,self.width,self.height)

        # Main widget
        wid = QWidget()
        self.setCentralWidget(wid)

        # Main layout
        mainLayout = QHBoxLayout()
        
        # main panel
        self.createSpectralPanel()

        # Add panel to main layout
        mainLayout.addWidget(self.spectralPanel)
        wid.setLayout(mainLayout)
        self.show()

        # Menu
        self.createMenu()

    def selTelescope(self):
        """ Select telescope """
        TD = selectTelescope()
        if TD.exec_() == QDialog.Accepted:
            self.telescope = TD.save()
            print('Selected telescope: ',self.telescope)
            #self.importTelescope()
            
    def importTelescope(self):
        """ Import library of the telescope """

        if self.telescope == 'WIYN':
            from showspectra.wiyn import getSky, getGalaxies, getErrors
        elif self.telescope == 'VIMOS':
            from showspectra.vimos import getSky, getGalaxies, getErrors
        else:
            print("Telescope "+self.telescope+" is not yet supported.")
        
    def createMenu(self):
        """ Menu """

        # File import/save
        bar = self.menuBar()
        file = bar.addMenu("File")
        file.addAction(QAction("Quit",self,shortcut='Ctrl+q',triggered=self.fileQuit))

        # Help 
        help = bar.addMenu("Help")

        bar.setNativeMenuBar(False)


    def createSpectralPanel(self):
        """ Panel to plot the spectrum """

        self.spectralPanel = QWidget()
        layout = QVBoxLayout(self.spectralPanel)

        # Plotting panel
        self.sp = SpectrumCanvas(self.spectralPanel,width=4, height=2.5,dpi=100)
        self.sp.toolbar = NavigationToolbar(self.sp, self)
        
        # Toolbar
        self.createToolbar()

        # Status bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to Show Spectra !", 10000)

        # Footer for the spectrum
        footer = QWidget()
        footer.layout = QHBoxLayout(footer)
        footer.layout.addWidget(self.sp.toolbar)
        footer.layout.addWidget(self.tb)
        footer.layout.addWidget(self.sb)
        
        # Add widgets to the panel
        layout.addWidget(self.sp)
        layout.addWidget(footer)
        
    def createToolbar(self):
        """ Toolbar with main commands """

        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')

        # Actions
        self.openAction = self.createAction(self.path0+'/icons/open.png','Open files','Ctrl+o',self.fileOpen)
        self.teleAction = self.createAction(self.path0+'/icons/telescope.png','Select telescope','Ctrl+T',self.selTelescope)
        self.quitAction = self.createAction(self.path0+'/icons/exit.png','Quit program','Ctrl+q',self.fileQuit)
        #self.helpAction = self.createAction(self.path0+'/icons/help.png','Help','Ctrl+q',self.onHelp)

        # Add actions
        self.tb.addAction(self.teleAction)
        self.tb.addAction(self.openAction)
        self.tb.addAction(self.quitAction)


    def createAction(self,icon,text,shortcut,action):
        act = QAction(QIcon(icon),text, self)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act

    # Actions
    
    def fileQuit(self):
        """ Quitting the program """
        self.close()

    def fileOpen(self):
        """ Opening spectral files """

        # Import telescope library
        if self.telescope == 'WIYN':
            from showspectra.wiyn import getSky, getGalaxies, getErrors
        elif self.telescope == 'VIMOS':
            from showspectra.vimos import getSky, getGalaxies, getErrors
        else:
            print("Telescope "+self.telescope+" is not yet supported.")

        FD = selectFiles()
        if FD.exec_() == QDialog.Accepted:
            flux,err,sky = FD.save()
            print (FD.save())

            # Opening files
            print('Opening ',flux)
            self.galaxies = getGalaxies(flux)
            self.dirname, file = os.path.split(flux)
            print('Opening ',err)
            getErrors(self, err)
            print('Opening ',sky)
            self.sky = getSky(sky)
            # Recover previous analysis
            if (os.path.exists(self.dirname+'/showspectra.fits')):
                print("Recovering previous analysis ...")
                recoverAnalysis(self)
            else:
                self.ngal = 0
                self.ngalaxies = len(self.galaxies)

            self.sp.computeInitialFigure(self)
        

#if __name__ == '__main__':
def main():
    app = QApplication(sys.argv)
    gui = GUI()
        
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(width*0.025, 0, width*0.95, width*0.45)
    # Add an icon for the application
    app.setWindowIcon(QIcon(QPixmap(gui.path0+'/icons/showspectra.png')))
    app.setApplicationName('SHOWSPECTRA')
    app.setApplicationVersion('0.01-beta')

    # Select telescope
    TD = selectTelescope()
    if TD.exec_() == QDialog.Accepted:
        gui.telescope = TD.save()
        print('Selected telescope: ',gui.telescope)

    sys.exit(app.exec_())
