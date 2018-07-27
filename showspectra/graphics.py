import numpy as np
import os
from PyQt5.QtWidgets import (QSizePolicy, QInputDialog)
from PyQt5.QtCore import Qt, QSize
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.font_manager import FontProperties
from matplotlib.widgets import Button
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.widgets import SpanSelector
# Matplotlib parameters
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
rcParams['font.size'] = 13
rcParams['mathtext.fontset'] = 'stix'
rcParams['legend.numpoints'] = 1
from scipy.signal import savgol_filter

class MplCanvas(FigureCanvas):
    """Basic matplotlib canvas class."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        FigureCanvas.updateGeometry(self)
        self.computeInitialFigure()

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(5, 5)

    def computeInitialFigure(self):
        pass


class SpectrumCanvas(MplCanvas):
    """ Canvas to plot spectra """

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

        # Import lines
        from showspectra.lines import define_lines
        self.Lines = define_lines()
        # Display defaults
        self.displayFlux = True
        self.displayErr = True
        self.displaySky = True
        self.displayLines = True
        self.displayMask = True
        self.displayTemplate = True
        self.xlimits = None
        self.ylimits = None

        # Activate focus
        self.setFocusPolicy(Qt.ClickFocus)
        self.setFocus()


    def computeInitialFigure(self, parent=None):

        if parent is None:
            pass
        else:
            try:
                self.fig.delaxes(self.axes)
                self.axes = None
            except BaseException:
                pass

            # Define figure
            self.fig.set_edgecolor('none')
            self.axes = self.fig.add_axes([0.04, 0.10, .85, .78], label='mainax')
            self.axes.format_coord = lambda x, y: "{:6.2f} \u212B  {:10.4e} W/m\u00B2/Hz".format(x, y)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.grid(True, which='both')

            # Next/Previous
            self.axprev = self.fig.add_axes([0.44, 0.85, 0.02, 0.04], label='axprev')
            self.axnext = self.fig.add_axes([0.47, 0.85, 0.02, 0.04], label='axnext')
            for a in ['top','bottom','right','left']:
                self.axprev.spines[a].set_visible(False)
                self.axnext.spines[a].set_visible(False)
            # self.bnext = Button(self.axnext, '>>', color='powderblue', hovercolor='lightskyblue')
            self.bnext = Button(self.axnext, '>>', color='white', hovercolor='white')
            self.bnext.on_clicked(self.nextspec)
            self.bprev = Button(self.axprev, '<<', color='white', hovercolor='white')
            self.bprev.on_clicked(self.prevspec)
            # Filter
            self.filter = False
            # self.axfilt = self.fig.add_axes([0.03, 0.80, 0.04, 0.04], label='axnext')
            # self.bfilt = Button(self.axfilt, 'Filter off', color='powderblue', 
            #                    hovercolor='lightskyblue')
            # self.bfilt.on_clicked(self.applyFilter)
            # Initial set for visibility
            self.showFlux = True
            self.showSky = True
            self.showErr = True
            self.showLines = True
            self.showMask = True

            # Plot spectrum
            self.drawSpectrum(parent)

            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()

    def drawSpectrum(self, parent):

        # self.axes.clear()  # leaves stuff around and the code crushes
        self.axes.cla()
        self.axes.set_xlabel('Wavelength [$\AA$]')
        self.axes.set_ylabel('Flux [$W\,m^{-2}\,Hz^{-1}$]')
        self.parent = parent
        self.ngal = parent.ngal
        self.gal = parent.galaxies[self.ngal]
        gal = parent.galaxies[self.ngal]
        self.wave = gal.wc / (1. + gal.z)
        if self.filter:
            flux = savgol_filter(gal.fc, 7, 3)
        else:
            flux = gal.fc
        good = self.gal.c == True
        self.fluxLine = self.axes.plot(self.wave[good], flux[good], label='Flux')
        self.galspec, = self.fluxLine
        self.axes.set_xlim([gal.xlim1, gal.xlim2])
        self.axes.set_ylim([gal.ylim1, gal.ylim2])
        # dx = gal.xlim2 - gal.xlim1
        dy = gal.ylim2 - gal.ylim1
        self.galspec.set_visible(self.showFlux)
        # Sky spectrum
        f = parent.sky.f
        fgal = parent.galaxies[parent.ngal].f
        f = (f - np.median(f)) / max(f) * gal.ylim2 + np.median(fgal) - gal.ylim1 / 10.
        self.skyLine = self.axes.plot(parent.sky.w / (1. + gal.z), f, color='r', label='Sky')
        self.skyspec, = self.skyLine
        self.skyspec.set_visible(self.showSky)
        # Spectrum error
        self.errLine = self.axes.plot(self.wave, gal.ec, color='g', label='Error')
        self.errspec, = self.errLine
        self.errspec.set_visible(self.showErr)
        # Fake line to have the lines in the legend
        self.linesLine = self.axes.plot([0, 0.1], [0, 0], color='purple',
                                        alpha=1.0, label='Lines', zorder=11)
        self.linesLayer, = self.linesLine
        # Fake line to have the lines in the legend
        if self.filter:
            filtlab = 'F on'
            alphaf = 1.0
        else:
            filtlab = 'F off'
            alphaf = 0.3
        self.filterLine = self.axes.plot([0, 0.1], [0, 0], ':', color='skyblue',
                                        alpha=alphaf, label=filtlab, zorder=11)
        self.filterLayer, = self.filterLine
        # Fake line to have masks in the legend
        if self.showMask:
            alpham = 1.0
        else:
            alpham = 0.3
        self.maskLine = self.axes.plot([0, 0.1], [0, 0], ':', color='yellow',
                                       alpha=alpham, label='Mask', zorder=11)
        self.maskLayer, = self.maskLine
        # Annotations
        # Number of galaxy spectrum
        self.sannotation = self.axes.annotate("{:3d}".format(self.ngal), xytext=(0.50, 0.98),
                                              xy=(0.50,0.98), ha='center', picker=5, 
                                              textcoords='axes fraction', xycoords='axes fraction')
        # Redshift
        self.zannotation = self.axes.annotate(" z = {:.4f}".format(self.gal.z), xy=(1.04, 0.75),
                                              picker=5, xycoords='axes fraction',ha='center')
        # Line names
        self.annotations = []
        font = FontProperties(family='DejaVu Sans', size=12)
        xlim0, xlim1 = self.axes.get_xlim()
        ylim0, ylim1 = self.axes.get_ylim()
        dy = ylim1 - ylim0
        for line in self.Lines.keys():
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]
            if (wline > xlim0 and wline < xlim1):
                wdiff = abs(self.wave - wline)
                imin = np.argmin(wdiff)
                y = flux[imin]
                y1 = y
                if nline[0:2] == 'A:':
                    y2 = ylim0 + 0.05 * dy
                    nline0 = nline[2:]
                else:
                    nline0 = nline
                    if (ylim1 - (y + 0.2 * dy)) > ((y - 0.2 * dy) - ylim0):
                        y2 = y + 0.2 * dy
                    else:
                        y2 = y - 0.2 * dy
                annotation = self.axes.annotate(nline0, xy=(wline, y1), xytext=(wline, y2),
                                                color='purple', alpha=0.4, visible=self.showLines,
                                                ha='center',
                                                arrowprops=dict(edgecolor='purple', facecolor='y',
                                                                arrowstyle='-', alpha=0.4,
                                                                connectionstyle='angle, angleA=0,' +
                                                                ' angleB=90, rad=10'),
                                                rotation=90, fontstyle='italic',
                                                fontproperties=font)
                annotation.draggable()
                self.annotations.append(annotation)
        # Prepare legend
        lns = self.fluxLine + self.errLine + self.skyLine + self.linesLine + \
            self.filterLine + self.maskLine
        self.lines = [self.galspec, self.errspec, self.skyspec, self.linesLayer,\
                      self.filterLayer,self.maskLayer]
        self.labs = [l.get_label() for l in lns]
        self.leg = self.axes.legend(lns, self.labs, loc='center right', bbox_to_anchor=(1.08, 0.2),
                                    fancybox=True, shadow=True, ncol=1)
        # Check alpha of filter
        self.leg.get_texts()[4].set_alpha(alphaf)
        self.leg.get_texts()[5].set_alpha(alpham)
        self.leg.draggable()
        self.lined = dict()
        self.labed = dict()
        for legline, origline, txt in zip(self.leg.get_lines(), self.lines, self.leg.texts):
            legline.set_picker(5)  # 5pts tolerance
            self.lined[legline] = origline
            self.labed[legline] = txt
        # Show mask
        if self.showMask:
            # Found the extremes of the masked patches
            nmask = self.gal.c == False
            if np.sum(nmask) > 0:
                c = self.gal.c.astype(int)
                c = np.append(c,1)
                dc = c[1:]-c[:-1]
                dc1 = dc == -1
                dc2 = dc == 1
                istart = np.where(dc1)
                if c[0] == 0:
                    istart.insert(0, 0)
                iend = np.where(dc2)
                # Mark rectangles
                istart = np.concatenate(istart); iend = np.concatenate(iend)
                print('Starts ',istart)
                print('Ends ',iend)
                print('There are ',np.size(istart),' masks')
                self.axrectangles = []
                for ist,ien in zip(istart,iend):
                    axrectangle = self.axes.axvspan(self.wave[ist],self.wave[ien],\
                                                    facecolor='LightYellow',\
                                                    alpha=1,linewidth=0,zorder=1)              
                    self.axrectangles.append(axrectangle)
            else:
                self.axrectangles = None
                print('No masked elements')
            
        # Connect canvas to events
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.dragged = None
        self.pick_pos = None
        self.draw_idle()
        # Start the span selector 
        self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                   rectprops=dict(alpha=0.5, facecolor='LightSalmon'))
        self.span.active = False

    def nextspec(self, event):
        if self.parent.ngal < (self.parent.ngalaxies - 1):
            self.parent.ngal += 1
            self.drawSpectrum(self.parent)
        else:
            print('There are no spectra left')

    def prevspec(self, event):
        if self.parent.ngal > 0:
            self.parent.ngal -= 1
            self.drawSpectrum(self.parent)
        else:
            print('First spectrum reached')
            
    def applyFilter(self, event):
        self.filter ^= True
        self.drawSpectrum(self.parent)
        #if self.filter:
        #    # self.bfilt.label.set_text('Filter on')
        #    self.leg.get_texts()[4].set_text('F on')
        #else:
        #    # self.bfilt.label.set_text('Filter off')
        #    self.leg.get_texts()[4].set_text('F off')

    def changeVisibility(self, label):
        print('Called change visibility')
        if label == 'Flux':
            self.showFlux ^= True
            self.galspec.set_visible(self.showFlux)
        elif label == 'Err':
            self.showErr ^= True
            self.errspec.set_visible(self.showErr)
        elif label == 'Sky':
            self.showSky ^= True
            self.skyspec.set_visible(self.showSky)
        elif label == 'Lines':
            self.showLines ^= True
            for annotation in self.annotations:
                annotation.set_visible(self.displayLines)
                
    def onSelect(self, xmin, xmax):
        """Define a masked region."""
        # Add region to clip mask
        indmin, indmax = np.searchsorted(self.wave, (xmin, xmax))
        print('indices ', indmin, indmax)
        self.gal.c[indmin:indmax]=0
        # Plot rectangle
        self.axes.axvspan(xmin,xmax,facecolor='LightYellow',alpha=1,linewidth=0,zorder=1)

    def onpick(self, event):
        if isinstance(event.artist, Line2D):
            legline = event.artist
            label = legline.get_label()
            if label == 'Lines':
                self.showLines ^= True
                for annotation in self.annotations:
                    annotation.set_visible(self.showLines)
                vis = self.showLines
                i=3
            elif label == 'Flux':
                self.showFlux ^= True
                self.galspec.set_visible(self.showFlux)
                vis = self.showFlux
                i=0
            elif label == 'Error':
                self.showErr ^= True
                self.errspec.set_visible(self.showErr)
                vis = self.showErr
                i=1
            elif label == 'Sky':
                self.showSky ^= True
                self.skyspec.set_visible(self.showSky)
                vis = self.showSky
                i=2
            elif (label == 'F on') | (label == 'F off'):
                self.applyFilter(event)
                vis = self.filter
                i=4
            elif label == 'Mask':
                self.showMask ^= True
                vis = self.showMask
                if self.axrectangles is not None:
                    for rect in self.axrectangles:
                        rect.set_visible(self.showMask)
                i=5
            else:
                print('Unknown label')
            # Transparency of legend
            if vis:
                alpha = 1.0
            else:
                alpha = 0.3
            legline.set_alpha(alpha)
            texts = self.leg.get_texts()
            texts[i].set_alpha(alpha)
            self.draw_idle()
        elif isinstance(event.artist, Text):
            # text = event.artist.get_text()
            if event.artist == self.zannotation:
                # c = 299792.458  # km/s
                znew = self.getDouble(self.gal.z)
                if znew is not None:
                    if znew != self.gal.z:
                        self.gal.z = znew
                        self.removeAnnotations()
                        self.drawSpectrum(self.parent)
            elif event.artist == self.sannotation:
                nnew = self.getInt(self.ngal)
                if nnew is not None:
                    if nnew != self.parent.ngal:
                        self.parent.ngal = nnew
                        self.removeAnnotations()
                        self.drawSpectrum(self.parent)
            else:
                self.dragged = event.artist
                self.pick_pos = event.mouseevent.xdata
        else:
            if event.artist == 'Legend':
                pass
            else:
                pass
                #print('Other pick event ', event.artist)
            pass
        return True

    def removeAnnotations(self):
        for annotation in self.annotations:
            annotation.remove()
        self.sannotation.remove()
        self.zannotation.remove()

    def getDouble(self, z):
        znew, okPressed = QInputDialog.getDouble(self, "Redshift", "z", z, -10000., 50000., 2)
        if okPressed:
            return znew
        else:
            return None

    def getInt(self, n):
        nnew, okPressed = QInputDialog.getInt(self, "Spectrum no", "n", n, 0,
                                              self.parent.ngalaxies, 2)
        if okPressed:
            return nnew
        else:
            return None

    def onrelease(self, event):
        # print('release event is ', event)
        if self.dragged is not None and self.pick_pos is not None:
            x1 = event.xdata
            x0 = self.pick_pos
            # w0 = np.array([self.Lines[line][1]  for  line in self.Lines.keys()])
            # wz = w0 * (1. + self.gal.z)
            z = (x1 - x0) / x0
            self.gal.z = (1. + self.gal.z) * (1 + z) - 1.
            for annotation in self.annotations:
                annotation.remove()
            self.zannotation.remove()
            self.drawSpectrum(self.parent)
            self.dragged = None
        return True


class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self, canvas, parent):
        self.iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")
        self.toolitems = [
            ('Home', 'Go back to original limits', 'home', 'home'),
            ('Pan', 'Pan figure', 'move', 'pan'),
            ('Zoom', 'Zoom in', 'zoom_to_rect', 'zoom'),
        ]
        self.parent = parent
        super().__init__(canvas, parent)
