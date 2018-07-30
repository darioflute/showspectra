import numpy as np
import os
from PyQt5.QtWidgets import (QSizePolicy, QInputDialog, QApplication)
from PyQt5.QtCore import Qt, QSize
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.font_manager import FontProperties
from matplotlib.widgets import Button
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.widgets import SpanSelector
from scipy.signal import savgol_filter
from scipy.ndimage.filters import median_filter
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

        # Activate focus
        # self.setFocusPolicy(Qt.ClickFocus)
        # self.setFocus()
        
    def computeInitialFigure(self, parent=None):

        if parent is None:
            pass
        else:
            try:
                self.fig.delaxes(self.axes)
                self.axes = None
            except BaseException:
                pass

            self.parent = parent
            # Define figure
            self.fig.set_edgecolor('none')
            self.axes = self.fig.add_axes([0.04, 0.10, .85, .78], label='mainax')
            self.axes.format_coord = lambda x, y: \
                "{:6.2f} \u212B {:10.4e} W/m\u00B2/Hz".format(x, y)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.grid(True, which='both')

            # Next/Previous
            self.axprev = self.fig.add_axes([0.44, 0.85, 0.02, 0.04], label='axprev')
            self.axnext = self.fig.add_axes([0.48, 0.85, 0.02, 0.04], label='axnext')
            for a in ['top', 'bottom', 'right', 'left']:
                self.axprev.spines[a].set_visible(False)
                self.axnext.spines[a].set_visible(False)
            self.bnext = Button(self.axnext, '>>', color='white', hovercolor='white')
            self.bnext.on_clicked(self.nextspec)
            self.bprev = Button(self.axprev, '<<', color='white', hovercolor='white')
            self.bprev.on_clicked(self.prevspec)
            # Filter
            self.filter = False
            # Initial set for visibility
            self.showFlux = True
            self.showSky = True
            self.showErr = True
            self.showLines = True
            self.showMask = True
            self.showTemplate = False

            # Plot spectrum
            self.drawSpectrum()

            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
            
            # Activate mouse wheel for zooming
            self.idw = self.mpl_connect('scroll_event', self.onWheel)
            # Activate pressing button event (for panning with middle button)
            self.mpl_connect('button_press_event', self.onPan)
            self.mpl_connect('button_release_event', self.onRelease)
            self.mpl_connect('motion_notify_event', self.onPan)
            self._event = None  # Event for panning

            # Start the span selector
            self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                     rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
            self.span.active = False


    def drawSpectrum(self):

        self.axes.cla()
        self.axes.set_xlabel('Wavelength [$\AA$]')
        self.axes.set_ylabel('Flux [$W\,m^{-2}\,Hz^{-1}$]')
        self.ngal = self.parent.ngal
        self.gal = self.parent.galaxies[self.ngal]
        gal = self.parent.galaxies[self.ngal]
        self.wave = gal.wc / (1. + gal.z)
        if self.filter:
            flux = savgol_filter(gal.fc, 7, 3)
        else:
            flux = gal.fc
        self.fluxLine = self.axes.plot(self.wave[self.gal.c], flux[self.gal.c], label='Flux')
        self.galspec, = self.fluxLine
        self.axes.set_xlim([gal.xlim1, gal.xlim2])
        self.axes.set_ylim([gal.ylim1, gal.ylim2])
        # dx = gal.xlim2 - gal.xlim1
        dy = gal.ylim2 - gal.ylim1
        self.galspec.set_visible(self.showFlux)
        # Sky spectrum
        f = self.parent.sky.f
        fgal = self.parent.galaxies[self.parent.ngal].f
        f = (f - np.median(f)) / max(f) * gal.ylim2 + np.median(fgal) - gal.ylim1 / 10.
        self.skyLine = self.axes.plot(self.parent.sky.w / (1. + gal.z), f, color='r', label='Sky', alpha=0.5)
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
                                              xy=(0.50, 0.98), ha='center', picker=5,
                                              textcoords='axes fraction', xycoords='axes fraction')
        # Redshift
        self.zannotation = self.axes.annotate(" z = {:.4f}".format(self.gal.z), xy=(1.04, 0.75),
                                              picker=5, xycoords='axes fraction', ha='center')
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
        self.lines = [self.galspec, self.errspec, self.skyspec, self.linesLayer,
                      self.filterLayer, self.maskLayer]
        # Redshift template
        if self.gal.zTemplate is not None:
            self.drawTemplate()
            lns += self.templateLine
            self.lines.append(self.templateLayer)

        self.labs = [l.get_label() for l in lns]
        self.leg = self.axes.legend(lns, self.labs, loc='center right', bbox_to_anchor=(1.10, 0.2),
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
            nmask = ~self.gal.c
            if np.sum(nmask) > 0:
                c = self.gal.c.astype(int)
                c = np.append(c, 1)
                dc = c[1:] - c[:-1]
                # dc1 = dc == -1
                # dc2 = dc == 1
                istart = np.where(dc == -1)
                if c[0] == 0:
                    istart = np.append(0, istart)
                iend = np.where(dc == 1)
                istart = np.ravel(istart)
                iend = np.ravel(iend)
                print('istart ',istart)
                print('iend ', iend)
                # Mark rectangles
                #istart = np.array(istart)
                #iend = np.array(iend)
                self.axrectangles = []
                for ist, ien in zip(istart, iend):
                    print(ist,ien)
                    axrectangle = self.axes.axvspan(self.wave[ist], self.wave[ien],
                                                    facecolor='LightYellow',
                                                    alpha=1, linewidth=0, zorder=1)
                    self.axrectangles.append(axrectangle)
            else:
                self.axrectangles = None
                print('No masked elements')

        # Connect canvas to events
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.dragged = None
        self.pick_pos = None
        self.draw_idle()
        
    def drawTemplate(self):
        """ Overplot the template """
        wg = self.gal.wc.copy()
        fg = self.gal.fc.copy()
        template = self.parent.templates[self.gal.zTemplate]
        z = self.gal.z
        bgr = median_filter(fg,128,mode='mirror')
        wt = template.w.copy()
        ft = template.f.copy()
        wg /= 1.+z
        minwg  = np.min(wg)
        maxwg  = np.max(wg)
        
        mask = (wt >= minwg) & (wt <= maxwg)
        wt = wt[mask]
        ft = ft[mask]
        fi = np.interp(wt, wg, fg)
        bgr = median_filter(fi,128,mode='mirror')
        fi -= bgr

        alpha = np.sum(fi*ft)/np.sum(ft*ft)
        ft *= alpha

        self.templateLine = self.axes.plot(wt,ft+bgr,color='Cyan',label='Template', alpha=0.4)
        self.templateLayer, = self.templateLine
        self.templateLayer.set_visible(self.showTemplate)

        
    def updateTemplate(self, newRow):
        """Manual selection of template."""
        self.gal.z = self.parent.zxcorr[newRow]
        self.gal.dz = self.parent.szxcorr[newRow]
        self.gal.zTemplate = self.parent.txcorr[newRow]
        self.drawSpectrum()


    def nextspec(self, event):
        if self.parent.ngal < (self.parent.ngalaxies - 1):
            self.parent.ngal += 1
            self.parent.galaxies[self.parent.ngal].limits()  ## recompute original limits
            self.drawSpectrum()
        else:
            print('There are no spectra left')

    def prevspec(self, event):
        if self.parent.ngal > 0:
            self.parent.ngal -= 1
            self.parent.galaxies[self.parent.ngal].limits()  ## recompute original limits
            self.drawSpectrum()
        else:
            print('First spectrum reached')

    def applyFilter(self, event):
        self.filter ^= True
        self.drawSpectrum()

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
        """
            Define a masked region.
            If SHIFT or CTRL (Command on Apple) is pressed, previous masks are deleted.
        """
        # Add region to clip mask
        indmin, indmax = np.searchsorted(self.wave, (xmin, xmax))
        # Check if SHIFT is on
        modifiers = QApplication.keyboardModifiers()
        if (modifiers == Qt.ShiftModifier) or (modifiers == Qt.ControlModifier):
            # Unmask
            self.gal.c[indmin:indmax] = 1
        else:
            # Mask
            self.gal.c[indmin:indmax] = 0
            # Plot rectangle
            # self.axes.axvspan(xmin, xmax, facecolor='LightYellow', alpha=1, linewidth=0, zorder=1)
        self.drawSpectrum()

    def onpick(self, event):
        if event.mouseevent.button != 1: return 
        if isinstance(event.artist, Line2D):
            legline = event.artist
            label = legline.get_label()
            if label == 'Lines':
                self.showLines ^= True
                for annotation in self.annotations:
                    annotation.set_visible(self.showLines)
                vis = self.showLines
                i = 3
            elif label == 'Flux':
                self.showFlux ^= True
                self.galspec.set_visible(self.showFlux)
                vis = self.showFlux
                i = 0
            elif label == 'Error':
                self.showErr ^= True
                self.errspec.set_visible(self.showErr)
                vis = self.showErr
                i = 1
            elif label == 'Sky':
                self.showSky ^= True
                self.skyspec.set_visible(self.showSky)
                vis = self.showSky
                i = 2
            elif (label == 'F on') | (label == 'F off'):
                self.applyFilter(event)
                vis = self.filter
                i = 4
            elif label == 'Mask':
                self.showMask ^= True
                vis = self.showMask
                if self.axrectangles is not None:
                    for rect in self.axrectangles:
                        rect.set_visible(self.showMask)
                i = 5
            elif label == 'Template':
                self.showTemplate ^= True
                vis = self.showTemplate
                i = 6
                self.templateLayer.set_visible(self.showTemplate)
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
            if event.artist == self.zannotation:
                znew = self.getDouble(self.gal.z)
                if znew is not None:
                    if znew != self.gal.z:
                        self.gal.z = znew
                        self.removeAnnotations()
                        self.drawSpectrum()
            elif event.artist == self.sannotation:
                nnew = self.getInt(self.ngal)
                if nnew is not None:
                    if nnew != self.parent.ngal:
                        self.parent.ngal = nnew
                        self.removeAnnotations()
                        self.drawSpectrum()
            else:
                self.dragged = event.artist
                self.pick_pos = event.mouseevent.xdata
        else:
            if event.artist == 'Legend':
                pass
            else:
                pass
            pass
        return True
    
    def onPan(self, event):
        """Routine to pan using the middle button of the mouse."""
        if event.button != 2: return
        if event.name == 'button_press_event':
            self._event = event
        elif event.name == 'button_release_event':
            self._event = None
        elif event.name == 'motion_notify_event':
            if self._event is None:
                return
            if event.x != self._event.x or event.y != self._event.y:
                pixel_to_data = self.axes.transData.inverted()
                data = pixel_to_data.transform_point((event.x, event.y))
                data_ = pixel_to_data.transform_point((self._event.x, self._event.y))
                dx = data[0]-data_[0]
                dy = data[1]-data_[1]
                if event.x != self._event.x:
                    xlim = self.axes.get_xlim()
                    self.axes.set_xlim(xlim[0]-dx, xlim[1]-dx)
                    self.gal.xlim1 = xlim[0]-dx
                    self.gal.xlim2 = xlim[1]-dx
                if event.ydata != self._event.ydata:
                    ylim = self.axes.get_ylim()
                    self.axes.set_ylim(ylim[0]-dy, ylim[1]-dy)
                    self.gal.ylim1 = ylim[0]-dy
                    self.gal.ylim2 = ylim[1]-dy
                self.draw_idle()
            self._event = event
    
    def onWheel(self, event):
        """Zoom-in and out by rolling the wheel of the mouse."""
        eb = event.button
        x0, y0 = event.xdata, event.ydata
        if eb == 'up':
            factor = 0.9
        elif eb == 'down':
            factor = 1.1
        else:
            factor = 1.0
        if event.key == 'control':
            xlim = self.axes.get_xlim()
            wl = (x0 - xlim[0]) * factor
            wr = (xlim[1] - x0) * factor
            self.axes.set_xlim(x0 - wl, x0 + wr)
            self.gal.xlim1 = x0 - wl
            self.gal.xlim2 = x0 + wr
        else:
            ylim = self.axes.get_ylim()
            hd = (y0 - ylim[0]) * factor
            hu = (ylim[1] - y0) * factor
            self.axes.set_ylim(y0 - hd, y0 + hu)
            self.gal.ylim1 = y0 - hd
            self.gal.ylim2 = y0 + hu
        self.draw_idle()

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

    def onRelease(self, event):
        if event.button == 1:
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
                self.drawSpectrum()
                self.dragged = None
                return True
            
            # Deselect pan & zoom options on mouse release
            if self.toolbar._active == "PAN":
                self.toolbar.pan()
            if self.toolbar._active == "ZOOM":
                self.toolbar.zoom()
        elif event.button == 2:
            self.onPan(event)

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
