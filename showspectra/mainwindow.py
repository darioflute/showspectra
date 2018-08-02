#!/usr/bin/env python

import sys
import os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QHBoxLayout,
                             QVBoxLayout, QStatusBar, QToolBar, QAction, QDialog)
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import Qt

from showspectra.graphics import SpectrumCanvas, NavigationToolbar
from showspectra.dialogs import selectTelescope, selectFiles, guessParams
from showspectra.inout import recoverAnalysis
from showspectra.templates import readTemplates
from showspectra.xcorr import cross_correlation
from showspectra.interactors import SegmentsSelector, SegmentsInteractor, LineInteractor

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
        self.setWindowIcon(QIcon(QPixmap(self.path0 + '/icons/showspectra.png')))

        # Style
        # Define style
        with open(self.path0 + '/stylesheet.css', "r") as fh:
            self.setStyleSheet(fh.read())

        self.initUI()
                

    def initUI(self):
        """User interface."""

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

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
        """Select telescope."""
        TD = selectTelescope()
        if TD.exec_() == QDialog.Accepted:
            self.telescope = TD.save()
            print('Selected telescope: ', self.telescope)

    # def importTelescope(self):
    #     """Import library of the telescope."""

    #     if self.telescope == 'WIYN':
    #         from showspectra.wiyn import getSky, getGalaxies, getErrors
    #     elif self.telescope == 'VIMOS':
    #         from showspectra.vimos import getSky, getGalaxies, getErrors
    #     else:
    #         print("Telescope " + self.telescope + " is not yet supported.")

    def createMenu(self):
        """Menu."""

        # File import/save
        bar = self.menuBar()
        file = bar.addMenu("File")
        file.addAction(QAction("Select telescope", self, shortcut='Ctrl+T', triggered=self.selTelescope))
        file.addAction(QAction("Open fits files", self, shortcut='Ctrl+o', triggered=self.fileOpen))
        file.addAction(QAction("Quit", self, shortcut='Ctrl+q', triggered=self.fileQuit))
        
        # Tools 
        tools = bar.addMenu("Tools")
        tools.addAction(QAction("Mask/unmask", self, shortcut='Ctrl+m', triggered=self.maskSpectrum))
        tools.addAction(QAction("Cross-correlate", self, shortcut='Ctrl+x', triggered=self.xcorrSpectrum))

        # Help
        # help = bar.addMenu("Help")
        bar.setNativeMenuBar(False)

    def createSpectralPanel(self):
        """Panel to plot the spectrum."""

        self.spectralPanel = QWidget()
        layout = QVBoxLayout(self.spectralPanel)

        # Plotting panel
        self.sp = SpectrumCanvas(self.spectralPanel, width=4, height=2.5, dpi=100)
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
        """Toolbar with main commands."""

        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')

        # Actions
        self.maskAction = self.createAction(self.path0 + '/icons/mask.png',
                                            'Mask/unmask spectrum', 'Ctrl+m', self.maskSpectrum,checkable=True)
        self.xcorrAction = self.createAction(self.path0 + '/icons/galaxy.png',
                                            'Cross-correlate with templates', 'Ctrl+x', self.xcorrSpectrum)
        self.guessAction = self.createAction(self.path0 + '/icons/gauss.png',
                                            'Guess continuum and lines', 'Ctrl+g', self.guessSpectrum)
        self.fitAction = self.createAction(self.path0 + '/icons/compute.png',
                                            'Fit continuum and lines', 'Ctrl+f', self.fitSpectrum)
        self.openAction = self.createAction(self.path0 + '/icons/open.png',
                                            'Open files', 'Ctrl+o', self.fileOpen)
        self.teleAction = self.createAction(self.path0 + '/icons/telescope.png',
                                            'Select telescope', 'Ctrl+T', self.selTelescope)
        self.quitAction = self.createAction(self.path0 + '/icons/exit.png',
                                            'Quit program', 'Ctrl+q', self.fileQuit)
        # self.helpAction = self.createAction(self.path0+'/icons/help.png',
        # 'Help','Ctrl+q',self.onHelp)

        # Add actions
        self.tb.addAction(self.maskAction)
        self.tb.addAction(self.xcorrAction)
        self.tb.addAction(self.guessAction)
        self.tb.addAction(self.fitAction)
        self.tb.addAction(self.teleAction)
        self.tb.addAction(self.openAction)
        self.tb.addAction(self.quitAction)

    def createAction(self, icon, text, shortcut, action, checkable=False):
        act = QAction(QIcon(icon), text, self,checkable=checkable)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act

    # Actions
    def fileQuit(self):
        """Quitting the program."""
        self.close()

    def maskSpectrum(self):
        """Mask spectrum mode."""
        try:
            self.sp.span.active ^= True
            if self.sp.span.active:
                self.sp.showMask = True
        except BaseException:
            print('No spectrum defined')
            pass
        
    def xcorrSpectrum(self):
        """Cross-correlate a spectrum with SDSS templates."""   
        self.sp.showTemplate = True
        cross_correlation(self)

    def guessSpectrum(self):
        """Create a guess of continuum and lines."""
        self.GP = guessParams()
        if self.GP.exec_() == QDialog.Accepted:
            cont, em, ab = self.GP.save()
            if cont == 'Constant':
                self.zeroDeg = True
            else:
                self.zeroDeg = False
            self.nem = em
            self.nab = ab
            try:
                self.onRemoveContinuum('segments deleted')
                self.onRemoveContinuum('line deleted')
            except BaseException:
                pass
            self.CS = SegmentsSelector(self.sp.axes,self.sp.fig, self.onContinuumSelect,zD=self.zeroDeg)
        else:
            return        
        
    def onContinuumSelect(self, verts):
        # Order the x coordinates of the verts
        x, y = zip(*verts)
        x = np.array(x); y = np.array(y)
        # Order increasing if wavelength, decreasing if frequency
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]
        verts = [(i,j) for (i,j) in zip(x,y)]

        SI = SegmentsInteractor(self.sp.axes, verts, self.zeroDeg)
        SI.modSignal.connect(self.onModifiedGuess)
        SI.mySignal.connect(self.onRemoveContinuum)
        self.sp.guess = SI
        # Once the continuum is selected draw guesses for the lines required
        # Check with one component only
        x0 = (x[1] + x[2]) * 0.5
        fwhm = (x[2] - x[1]) / 2.
        # Find the max of the spectrum in the interval
        idx = (self.sp.wave > x[1]) & (self.sp.wave < x[2])
        A = np.nanmax(self.sp.flux[idx]) - self.sp.guess.intcpt - self.sp.guess.slope * x0
        LI = LineInteractor(self.sp.axes, self.sp.guess.intcpt, self.sp.guess.slope, x0, A, fwhm)
        LI.modSignal.connect(self.onModifiedGuess)
        LI.mySignal.connect(self.onRemoveContinuum)
        # This step will have a list of lines 
        self.sp.line = LI

    def onModifiedGuess(self, event):
        """Reacts to modifications of the continuum guess."""
        if event == 'continuum guess modified':
            if self.sp.line is not None:
                # Change continuum data in the line
                oldc = self.sp.line.c0 + self.sp.line.cs * self.sp.line.x0
                self.sp.line.c0 = self.sp.guess.intcpt
                self.sp.line.cs = self.sp.guess.slope
                newc = self.sp.line.c0 + self.sp.line.cs * self.sp.line.x0
                if np.abs(oldc - newc) > 0.:
                    if self.sp.line.A >= 0:
                        self.sp.line.A += oldc - newc
                    else:
                        self.sp.line.A -= oldc - newc
                self.sp.line.grab_background()
                self.sp.line.updateCurves()
        else:
            pass

    def onRemoveContinuum(self, event):
        if event == 'segments deleted':
            self.sp.guess.disconnect()  
            self.sp.guess = None
            self.sp.fig.canvas.draw_idle()
        elif event == 'line deleted':
            self.sp.line.disconnect()
            self.sp.line = None
            self.sp.fig.canvas.draw_idle()

    def fitSpectrum(self):
        """Fit defined guess."""
        print('Fit the defined guess')

    def fileOpen(self):
        """Opening spectral files."""

        # Import telescope library
        if self.telescope == 'WIYN':
            from showspectra.wiyn import getSky, getGalaxies, getErrors
        elif self.telescope == 'VIMOS':
            from showspectra.vimos import getSky, getGalaxies, getErrors
        else:
            print("Telescope " + self.telescope + " is not yet supported.")

        FD = selectFiles()
        if FD.exec_() == QDialog.Accepted:
            flux, err, sky = FD.save()
            print(FD.save())

            # Opening files
            print('Opening ', flux)
            self.galaxies = getGalaxies(flux)
            self.dirname, file = os.path.split(flux)
            print('Opening ', err)
            getErrors(self, err)
            print('Opening ', sky)
            self.sky = getSky(sky)
            # Recover previous analysis
            if (os.path.exists(self.dirname + '/showspectra.fits')):
                print("Recovering previous analysis ...")
                recoverAnalysis(self)
            else:
                self.ngal = 0
                self.ngalaxies = len(self.galaxies)
            self.sp.computeInitialFigure(self)


def main():
    app = QApplication(sys.argv)
    gui = GUI()

    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(width * 0.025, 0, width * 0.95, width * 0.45)
    # Add an icon for the application
    app.setWindowIcon(QIcon(QPixmap(gui.path0 + '/icons/showspectra.png')))
    app.setApplicationName('SHOWSPECTRA')
    app.setApplicationVersion('0.01-beta')

    # Select telescope and open first files
    gui.selTelescope()
    gui.fileOpen()
    # Read spectral templates into the gui
    readTemplates(gui)

    sys.exit(app.exec_())
