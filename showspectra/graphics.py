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
# from matplotlib.widgets import SpanSelector
# Matplotlib parameters
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
rcParams['font.size'] = 13
rcParams['mathtext.fontset'] = 'stix'
rcParams['legend.numpoints'] = 1


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
            self.axes = self.fig.add_axes([0.10, 0.10, .85, .78], label='mainax')
            self.axes.format_coord = lambda x, y: "{:8.4f} A  {:10.4f} W/m2/Hz".format(x, y)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.grid(True, which='both')

            # Next/Previous
            self.axprev = self.fig.add_axes([0.03, 0.3, 0.04, 0.04], label='axprev')
            self.axnext = self.fig.add_axes([0.03, 0.35, 0.04, 0.04], label='axnext')
            self.bnext = Button(self.axnext, '>>', color='powderblue', hovercolor='lightskyblue')
            self.bnext.on_clicked(self.next)
            self.bprev = Button(self.axprev, '<<', color='powderblue', hovercolor='lightskyblue')
            self.bprev.on_clicked(self.prev)

            # Initial set for visibility
            self.showFlux = True
            self.showSky = True
            self.showErr = True
            self.showLines = True

            # Plot spectrum
            self.drawSpectrum(parent)

            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()

    def drawSpectrum(self, parent):

        self.axes.clear()
        self.parent = parent
        self.ngal = parent.ngal
        self.gal = parent.galaxies[self.ngal]
        gal = parent.galaxies[self.ngal]
        # print('z is ', gal.z)
        wave = gal.wc / (1. + gal.z)
        flux = gal.fc
        # if self.filter:
        #    flux = savgol_filter(flux, 7, 3)
        self.fluxLine = self.axes.plot(wave, flux, label='Flux')
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
        self.errLine = self.axes.plot(wave, gal.ec, color='g', label='Err')
        self.errspec, = self.errLine
        self.errspec.set_visible(self.showErr)

        # Fake line to have the lines in the legend
        self.linesLine = self.axes.plot([0, 0.1], [0, 0], color='purple',
                                        alpha=0.4, label='Lines', zorder=11)
        self.linesLayer, = self.linesLine
        self.zannotation = self.axes.annotate(" z = {:.4f}".format(self.gal.z), xy=(-0.1, 0.9),
                                              picker=5, xycoords='axes fraction')
        # Lines
        self.annotations = []
        font = FontProperties(family='DejaVu Sans', size=12)
        xlim0, xlim1 = self.axes.get_xlim()
        ylim0, ylim1 = self.axes.get_ylim()
        dy = ylim1 - ylim0
        for line in self.Lines.keys():
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]
            if (wline > xlim0 and wline < xlim1):
                wdiff = abs(wave - wline)
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
                                                arrowprops=dict(edgecolor='purple', facecolor='y',
                                                                arrowstyle='-', alpha=0.4,
                                                                connectionstyle='angle, angleA=0,' +
                                                                ' angleB=90, rad=10'),
                                                rotation=90, fontstyle='italic',
                                                fontproperties=font)
                annotation.draggable()
                self.annotations.append(annotation)
        # Prepare legend
        lns = self.fluxLine + self.errLine + self.skyLine + self.linesLine
        lines = [self.galspec, self.errspec, self.skyspec, self.linesLayer]
        self.labs = [l.get_label() for l in lns]
        leg = self.axes.legend(lns, self.labs, loc='center left', bbox_to_anchor=(-0.13, 0.6),
                               fancybox=True, shadow=True, ncol=1)
        leg.draggable()
        self.lined = dict()
        self.labed = dict()
        for legline, origline, txt in zip(leg.get_lines(), lines, leg.texts):
            legline.set_picker(5)  # 5pts tolerance
            self.lined[legline] = origline
            self.labed[legline] = txt

        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.dragged = None
        self.pick_pos = None
        self.draw_idle()

    def next(self, event):
        self.parent.ngal += 1
        self.drawSpectrum(self.parent)

    def prev(self, event):
        if self.parent.ngal > 0:
            self.parent.ngal -= 1
            self.drawSpectrum(self.parent)

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

    def onpick(self, event):
        if isinstance(event.artist, Line2D):
            legline = event.artist
            label = legline.get_label()
            if label == 'Lines':
                self.showLines ^= True
                for annotation in self.annotations:
                    annotation.set_visible(self.showLines)
            elif label == 'Flux':
                self.showFlux ^= True
                self.galspec.set_visible(self.showFlux)
            elif label == 'Err':
                self.showErr ^= True
                self.errspec.set_visible(self.showErr)
            elif label == 'Sky':
                self.showSky ^= True
                self.skyspec.set_visible(self.showSky)
            else:
                origline = self.lined[legline]
                # txt = self.labed[legline]
                vis = not origline.get_visible()
                # print("label ",label," origline ",origline,'txt ',txt)
                origline.set_visible(vis)
            self.draw_idle()
        elif isinstance(event.artist, Text):
            # text = event.artist.get_text()
            if event.artist == self.zannotation:
                # c = 299792.458  # km/s
                znew = self.getDouble(self.gal.z)
                if znew is not None:
                    if znew != self.gal.z:
                        self.gal.z = znew
                        for annotation in self.annotations:
                            annotation.remove()
                        self.zannotation.remove()
                        self.drawSpectrum(self.parent)
            else:
                self.dragged = event.artist
                self.pick_pos = event.mouseevent.xdata
        else:
            print('Other pick event ', event.artist)
            pass
        return True

    def getDouble(self, z):
        znew, okPressed = QInputDialog.getDouble(self, "Redshift", "z", z, -10000., 50000., 2)
        if okPressed:
            return znew
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
