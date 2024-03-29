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
from showspectra.inout import importAnalysis, exportAnalysis
from showspectra.templates import readTemplates
from showspectra.xcorr import cross_correlation
from showspectra.interactors import (SegmentsSelector, SegmentsInteractor, 
                                     LineInteractor, LineManager)
from showspectra.lines import fitContinuum, fitLines
from showspectra.spectra import Line


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
        with open(self.path0 + '/stylesheet.css', "r") as fh:
            self.setStyleSheet(fh.read())
        # Start UI
        self.initUI()

    def initUI(self):
        """User interface."""
        # Read spectral templates
        # readTemplates(self)
        # Start user interface
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # Main widget
        wid = QWidget()
        self.setCentralWidget(wid)
        # Main layout
        mainLayout = QHBoxLayout()
        # Main panel
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
            self.fileOpen()

    def createMenu(self):
        """Menu."""
        # File import/save
        bar = self.menuBar()
        file = bar.addMenu("File")
        file.addAction(QAction("Select telescope", self, shortcut='Ctrl+T',
                               triggered=self.selTelescope))
        file.addAction(QAction("Open fits files", self, shortcut='Ctrl+o', triggered=self.fileOpen))
        file.addAction(QAction("Quit", self, shortcut='Ctrl+q', triggered=self.fileQuit))
        # Tools
        tools = bar.addMenu("Tools")
        tools.addAction(QAction("Mask/unmask sky", self, shortcut='Ctrl+m',
                                triggered=self.maskSpectrum))
        tools.addAction(QAction("Mask/unmask glitch", self, shortcut='Ctrl+g',
                                triggered=self.maskGlitch))
        tools.addAction(QAction("Cross-correlate", self, shortcut='Ctrl+x',
                                triggered=self.xcorrSpectrum))
        # Help
        # help = bar.addMenu("Help")
        bar.setNativeMenuBar(False)

    def createSpectralPanel(self):
        """Panel to plot the spectrum."""
        self.spectralPanel = QWidget()
        layout = QVBoxLayout(self.spectralPanel)
        # Plotting panel
        self.sp = SpectrumCanvas(self.spectralPanel, width=4, height=2.5, dpi=100)
        self.sp.maskSignal.connect(self.maskOtherSpectra)
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
        iconpath = self.path0 + '/icons/'
        self.maskAction = self.createAction(iconpath + 'masksky.png', 'Mask/unmask sky features', 'Ctrl+m',
                                            self.maskSpectrum, checkable=True)
        self.glitchAction = self.createAction(iconpath + 'maskglitch.png', 'Mask/unmask glitch', 'Ctrl+g',
                                            self.maskGlitch, checkable=True)
        self.xresizeAction = self.createAction(iconpath + 'wresize.png', 'Resize wavelength',
                                               'Ctrl+r', self.resizeWavelength)
        self.zresetAction = self.createAction(iconpath + 'zreset.png', 'Reset redshift',
                                              'Ctrl+r', self.resetRedshift)
        self.xcorrAction = self.createAction(iconpath + 'galaxy.png',
                                             'Cross-correlate with templates',
                                             'Ctrl+x', self.xcorrSpectrum)
        self.idlineAction = self.createAction(iconpath + 'identifyline.png', 'Identify line',
                                              'Ctrl+l', self.identifyLine)
        self.guessAction = self.createAction(iconpath + 'guess.png', 'Guess continuum and lines',
                                             'Ctrl+g', self.guessSpectrum)
        self.fitAction = self.createAction(iconpath + 'fitline.png', 'Fit continuum and lines',
                                           'Ctrl+f', self.fitSpectrum)
        self.removeAction = self.createAction(iconpath + 'remove.png', 'Remove fitted line',
                                              'Ctrl+r', self.removeFittedLine)
        self.openAction = self.createAction(iconpath + 'open.png', 'Open files', 'Ctrl+o',
                                            self.fileOpen)
        self.teleAction = self.createAction(iconpath + 'telescope.png', 'Select telescope',
                                            'Ctrl+T', self.selTelescope)
        self.quitAction = self.createAction(iconpath + 'exit.png',
                                            'Quit program', 'Ctrl+q', self.fileQuit)
        # self.helpAction = self.createAction(iconpath+'help.png', 'Help', 'Ctrl+q', self.onHelp)
        # Add actions
        self.tb.addAction(self.maskAction)
        self.tb.addAction(self.glitchAction)
        self.tb.addAction(self.xresizeAction)
        self.tb.addAction(self.zresetAction)
        self.tb.addAction(self.xcorrAction)
        self.tb.addAction(self.idlineAction)
        self.tb.addAction(self.guessAction)
        self.tb.addAction(self.fitAction)
        self.tb.addAction(self.removeAction)
        self.tb.addAction(self.teleAction)
        self.tb.addAction(self.openAction)
        self.tb.addAction(self.quitAction)

    def createAction(self, icon, text, shortcut, action, checkable=False):
        act = QAction(QIcon(icon), text, self, checkable=checkable)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act

    def resizeWavelength(self):
        """Resize to full wavelength."""
        self.sp.gal.limits()
        self.sp.drawSpectrum()

    def resetRedshift(self):
        """Reset redshift."""
        self.sp.gal.z = 0.0
        self.sp.gal.dz = 0.0
        self.sp.zTemplate = None
        self.sp.gal.quality = '?'
        if self.sp.gal.spectype != 'sky':
            self.sp.gal.spectype = '?'
        self.resizeWavelength()
        # Save new analysis
        exportAnalysis(self.galaxies, self.ngal, self.dirname)

    def removeFittedLine(self):
        """Program awaits for deleting a line with a click of the mouse."""
        self.sp.removeFittedLine = True

    def identifyLine(self):
        """Identify a line by eye."""
        self.sp.identifyLine = True

    def fileQuit(self):
        """Quitting the program."""
        if self.ngal is not None:
            exportAnalysis(self.galaxies, self.ngal, self.dirname)
        self.close()


    def maskSpectrum(self):
        self.maskSpectra(False)

    def maskSpectra(self, glitch):
        """Mask spectrum mode."""
        # print('glitch is ', glitch)
        # print('glitch masking ', self.glitchAction.isChecked())
        # print('sky masking ', self.maskAction.isChecked())
        try:
            if self.sp.maskGlitch and glitch:
                self.sp.span.active ^= True
            elif self.sp.maskSky and not glitch :
                self.sp.span.active ^= True
            elif self.sp.maskGlitch == False and self.sp.maskSky == False:
                self.sp.span.active ^= True
            if self.sp.span.active:
                if glitch == True:
                    self.sp.maskGlitch = True
                    self.sp.maskSky = False
                    self.glitchAction.setChecked(True)
                    self.maskAction.setChecked(False)
                else:
                    self.sp.maskSky = True
                    self.sp.maskGlitch = False
                    self.maskAction.setChecked(True)
                    self.glitchAction.setChecked(False)
                if self.sp.showLines:
                    self.sp.changeVisibility('Lines')
                self.sp.showMask = True
                self.sp.setMaskVisibility(self.sp.showMask)
                self.sp.leg.get_texts()[5].set_alpha(1.0)
                self.sp.leg.get_lines()[5].set_alpha(1.0)
                self.sp.leg.get_texts()[3].set_alpha(0.3)
                self.sp.leg.get_lines()[3].set_alpha(0.3)
            else:
                if glitch == True:
                    self.sp.maskGlitch = False
                    self.glitchAction.setChecked(False)
                else:
                    self.sp.maskSky = False
                    self.maskAction.setChecked(False)
                if not self.sp.showLines:
                    self.sp.changeVisibility('Lines')
                self.sp.leg.get_texts()[3].set_alpha(1.0)
                self.sp.leg.get_lines()[3].set_alpha(1.0)
                pass
            self.sp.draw_idle()
        except BaseException:
            print('No spectrum defined')
            pass
        
    def maskGlitch(self):
        self.maskSpectra(True)

    def maskOtherSpectra(self, message):
        """Spread to the other spectra the action on one spectrum."""
        indmin = self.sp.masklimits[0]
        indmax = self.sp.masklimits[1]
        if message == 'mask':
            for galaxy in self.galaxies:
                galaxy.c[indmin:indmax] = 0
                galaxy.clip()
        elif message == 'unmask':
            for galaxy in self.galaxies:
                galaxy.c[indmin:indmax] = 1
                galaxy.clip()
        # Export the new analysis - to save stuff in case of crash
        exportAnalysis(self.galaxies, self.ngal, self.dirname)

    def xcorrSpectrum(self):
        """Cross-correlate a spectrum with SDSS templates."""
        self.sp.showTemplate = True
        xc = cross_correlation(self)
        if xc == 1:
            self.sp.gal.quality = 'XCorr'
            self.sp.gal.spectype = 'Galaxy'
            self.sp.updateQualityAnnotation()
            # Export the new analysis - to save stuff in case of crash
            exportAnalysis(self.galaxies, self.ngal, self.dirname)
        else:
            self.sp.removeTemplate()

    def guessSpectrum(self):
        """Create a guess of continuum and lines."""
        # First close any masking 
        if self.sp.maskGlitch:
            self.maskGlitch()
        elif self.sp.span.active:
            self.maskSpectra(False)
        # Then ask for type of fit
        if self.sp.gal.quality == "?":
            message = 'Please, find first the redshift of the spectrum !'
            self.sb.showMessage(message, 10000)
            print(message)
            return

        self.GP = guessParams()
        if self.GP.exec_() == QDialog.Accepted:
            if self.sp.showLines:
                self.sp.changeVisibility('Lines')
                self.sp.leg.get_texts()[3].set_alpha(0.3)
                self.sp.leg.get_lines()[3].set_alpha(0.3)
            cont, em, ab = self.GP.save()
            if cont == 'Constant':
                self.zeroDeg = True
            else:
                self.zeroDeg = False
            self.nem = int(em)
            self.nab = int(ab)
            try:
                self.onRemoveContinuum('all')
                self.sp.fig.canvas.draw_idle()
            except BaseException:
                pass
            self.sp.guessContinuum = True
            self.sp.span.active = False
            self.maskAction.setChecked(False)
            self.CS = SegmentsSelector(self.sp.axes, self.sp.fig, self.onContinuumSelect,
                                       zD=self.zeroDeg)
        else:
            return

    def onContinuumSelect(self, verts):
        # Order the x coordinates of the verts
        x, y = zip(*verts)
        x = np.array(x)
        y = np.array(y)
        # Order increasing if wavelength, decreasing if frequency
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]
        verts = [(i, j) for (i, j) in zip(x, y)]
        self.sp.guess = SegmentsInteractor(self.sp.axes, verts, self.zeroDeg)
        self.sp.guess.modSignal.connect(self.onModifiedGuess)
        self.sp.guess.mySignal.connect(self.onRemoveContinuum)
        interactors = [self.sp.guess]
        # Once the continuum is selected draw guesses for the lines required
        if self.nem > 0:
            self.sp.emlines = self.addLines(self.nem, x, 'emission')
            interactors.extend(self.sp.emlines)
        else:
            self.sp.emlines = []
        if self.nab > 0:
            self.sp.ablines = self.addLines(self.nab, x, 'absorption')
            interactors.extend(self.sp.ablines)
        else:
            self.sp.ablines = []
        self.sp.lineManager = LineManager(self.sp.axes, interactors)
        self.sp.guessContinuum = False
        if not self.sp.showLines:
            self.sp.changeVisibility('Lines')
            self.sp.leg.get_texts()[3].set_alpha(1.0)
            self.sp.leg.get_lines()[3].set_alpha(1.0)

    def addLines(self, n, x, type):
        lines = []
        for i in range(n):
            dx = (x[2] - x[1]) / (2 * n)
            x0 = x[1] + dx + i * 2 * dx
            fwhm = dx
            idx = (self.sp.wave > (x0 - dx)) & (self.sp.wave < (x0 + dx))
            y0 = self.sp.guess.intcpt + self.sp.guess.slope * x0
            if type == 'emission':
                A = np.nanmax(self.sp.flux[idx]) - y0
            else:
                A = np.nanmin(self.sp.flux[idx]) - y0
            LI = LineInteractor(self.sp.axes, y0,
                                self.sp.guess.slope, x0, A, fwhm)
            LI.modSignal.connect(self.onModifiedGuess)
            LI.mySignal.connect(self.onRemoveContinuum)
            lines.append(LI)
        return lines

    def onModifiedGuess(self, event):
        """Reacts to modifications of the guess."""
        if event == 'continuum guess modified':
            if self.sp.emlines is not None:
                self.sp.modifyGuess = True
                # Change continuum data in the line
                for line in self.sp.emlines + self.sp.ablines:
                    oldc = line.c0
                    line.c0 = self.sp.guess.intcpt + line.x0 * self.sp.guess.slope
                    line.cs = self.sp.guess.slope
                    newc = line.c0
                    if np.abs(oldc - newc) > 0.:
                        if line.A >= 0:
                            line.A += oldc - newc
                        else:
                            line.A -= oldc - newc
                    line.updateCurves()
        elif event == 'line guess modified':
            # Redraw lines after modifying one line
            print('line modified')
        else:
            pass

    def onRemoveContinuum(self, event):
        if event == 'segments deleted':
            self.sp.guess.disconnect()
            self.sp.guess = None
            self.sp.fig.canvas.draw_idle()
        elif event == 'line deleted':
            for line in self.sp.emlines + self.sp.ablines:
                line.disconnect()
                line = None
            self.sp.emlines = []
            self.sp.ablines = []
            self.sp.fig.canvas.draw_idle()
        elif event == 'all':
            self.sp.guess.disconnect()
            self.sp.guess = None
            for line in self.sp.emlines + self.sp.ablines:
                line.disconnect()
                line = None
            self.sp.lineManager.disconnect()
            self.sp.emlines = []
            self.sp.ablines = []
            #self.sp.fig.canvas.draw_idle()

    def fitSpectrum(self):
        """Fit defined guess."""
        # Continuum  self.sp.guess
        # Lines:     self.sp.emlines+self.sp.ablines
        intercept, interceptErr, slope, slopeErr = fitContinuum(self.sp)
        linefitpars = fitLines(self.sp, intercept, slope)
        z = self.sp.gal.z
        xg, yg = zip(*self.sp.guess.xy)
        xg = np.array(xg) * (1. + self.sp.gal.z)
        LinesNames = list(self.sp.Lines.keys())
        LinesCenters = [list(self.sp.Lines.values())[i][1] for i in range(len(LinesNames))]
        # print('Cont. ', intercept, slope)
        for pars in linefitpars:
            loc, locErr, amp, ampErr, scale, scaleErr = pars
            # print('loc, scale, amp ', loc, scale, amp)
            # print('redshift is: ',z)
            dl = np.abs(LinesCenters - loc / (1. + z))
            idx, = np.where(dl == min(dl))
            if len(idx) > 1:
                if amp > 0:
                    linename = LinesNames[idx[0]]  # Emission
                else:
                    linename = LinesNames[idx[1]]  # Absorption
            else:
                linename = LinesNames[idx[0]]
            # print('line is: ', linename)
            center = LinesCenters[idx[0]]
            line = Line(xg[0], xg[3], center, intercept, interceptErr, 
                        slope, slopeErr, loc, locErr, scale, scaleErr, amp, ampErr)
            # Add to dictionary
            self.sp.gal.lines[linename] = line
        # Take out guess
        self.onRemoveContinuum('all')
        # Draw fitted lines and continuum
        self.sp.drawSpectrum()
        # Export the new analysis - to save stuff in case of crash
        exportAnalysis(self.galaxies, self.ngal, self.dirname)

    def fileOpen(self):
        """Opening spectral files."""
        # Import telescope library
        if self.telescope == 'WIYN':
            from showspectra.wiyn import getSky, getGalaxies, getErrors
        elif self.telescope == 'VIMOS':
            from showspectra.vimos import getSky, getGalaxies, getErrors
        elif self.telescope == 'SDSS':
            from showspectra.sdss import getGalaxies
        elif self.telescope == 'MUSE':
            from showspectra.muse import getGalaxies
        else:
            print("Telescope " + self.telescope + " is not yet supported.")
        FD = selectFiles()
        if FD.exec_() == QDialog.Accepted:
            flux, err, sky = FD.save()
            print(FD.save())
            # Opening files
            self.ngal = None
            print('Opening ', flux)
            try:
                self.dirname, file = os.path.split(flux)
                if self.telescope in ['VIMOS','WIYN']:
                    self.galaxies = getGalaxies(flux)
                    print('Opening ', err)
                    getErrors(self, err)
                    print('Opening ', sky)
                    self.sky = getSky(sky)
                elif self.telescope in ['SDSS','MUSE']:
                    self.galaxies, self.sky = getGalaxies(flux)
            except BaseException:
                print('Files are not correct for telescope ', self. telescope)
                return
            # Recover previous analysis
            # if (os.path.exists(self.dirname + '/showspectra.fits')):
            #    print("Recovering previous analysis ...")
            #    recoverAnalysis(self)
            if (os.path.exists(self.dirname + '/showspectra.json')):
                print("Recovering previous analysis ...")
                analysisFile = self.dirname + '/showspectra.json'
                self.ngal, self.ngalaxies, self.galaxies = importAnalysis(analysisFile,
                                                                          self.galaxies)
            else:
                self.ngal = 0
                self.ngalaxies = len(self.galaxies)
            self.sp.computeInitialFigure(self)
            self.resizeWavelength()


def main():
    from showspectra import __version__
    app = QApplication(sys.argv)
    gui = GUI()
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(width * 0.025, 0, width * 0.95, width * 0.45)
    # Add an icon for the application
    app.setWindowIcon(QIcon(QPixmap(gui.path0 + '/icons/showspectra.png')))
    app.setApplicationName('SHOWSPECTRA')
    app.setApplicationVersion(__version__)
    # Read templates
    readTemplates(gui)
    # Select telescope and open first files
    gui.selTelescope()
    # gui.fileOpen()
    # Exec app
    sys.exit(app.exec_())
