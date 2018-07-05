import numpy as np
import os
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT 
from matplotlib.ticker import ScalarFormatter
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon,FancyArrowPatch
#from matplotlib.colorbar import Colorbar

# Matplotlib parameters
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

from matplotlib.lines import Line2D
from matplotlib.text import Text
#from matplotlib.colors import Normalize
#import matplotlib.cbook as cbook

from matplotlib.widgets import SpanSelector

from PyQt5.QtWidgets import (QVBoxLayout, QSizePolicy, QInputDialog, QDialog, QListWidget,
                             QListWidgetItem,QPushButton,QLabel,QMessageBox,QScrollArea,QWidget)
from PyQt5.QtGui import QIcon,QFont
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtTest import QTest




class MplCanvas(FigureCanvas):
    """ Basic matplotlib canvas class """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,QSizePolicy.MinimumExpanding,QSizePolicy.MinimumExpanding)
        FigureCanvas.updateGeometry(self)
        self.computeInitialFigure()

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(5,5)
    
    def computeInitialFigure(self):
        pass


class SpectrumCanvas(MplCanvas):
    """ Canvas to plot spectra """

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self,*args, **kwargs)

        # Import lines
        from lines import define_lines
        self.Lines = define_lines()


        # Display defaults
        self.displayFlux     = True
        self.displayErr      = True
        self.displaySky      = True
        self.displayLines    = True
        self.displayMask     = True
        self.displayTemplate = True
        self.xlimits = None
        self.ylimits = None

    def computeInitialFigure(self, parent=None):

        # Refresh axes
        #self.axes.clear()
        
        if parent is None:
            pass
        else:
            try:
                self.fig.delaxes(self.axes)
                self.axes = None
            except:
                pass
            # Focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
            
            # Define figure
            self.fig.set_edgecolor('none')
            self.axes = self.fig.add_axes([0.12,0.15,.8,.78])
            self.axes.format_coord = lambda x,y: "{:8.4f} A  {:10.4f} W/m2/Hz".format(x,y)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            
            # Plot spectrum
            gal = parent.galaxies[parent.ngal]
            wave = gal.wc/(1.+gal.z)
            flux = gal.fc
            #if self.filter:
            #    flux = savgol_filter(flux, 7, 3)
            self.galspec, = self.axes.plot(wave,flux)
            self.axes.set_xlim([gal.xlim1,gal.xlim2])
            self.axes.set_ylim([gal.ylim1,gal.ylim2])
            #drawLines(self,wave,flux)
            # Sky spectrum
            f = parent.sky.f
            fgal = parent.galaxies[parent.ngal].f
            f=(f-np.median(f))/max(f)*gal.ylim2+np.median(fgal)-gal.ylim1/10.
            self.skyspec, = self.axes.plot(parent.sky.w/(1.+gal.z),f,color='r')
            self.skyspec.set_visible(parent.showSky)
            # Spectrum error
            self.errspec, = self.axes.plot(wave,gal.ec,color='g')
            self.errspec.set_visible(parent.showErr)

            # Lines
            self.annotations = []
            font = FontProperties(family='DejaVu Sans', size=12)
            xlim0,xlim1 = self.axes.get_xlim()
            ylim0,ylim1 = self.axes.get_ylim()
            dy = ylim1-ylim0

            for line in self.Lines.keys():
                nline = self.Lines[line][0]
                wline = self.Lines[line][1]
                if (wline > xlim0 and wline < xlim1):
                    wdiff = abs(wave - wline)
                    imin = np.argmin(wdiff)
                    y = flux[imin]
                    y1 = y
                    if (ylim1 - (y+0.2*dy)) > ((y-0.2*dy)-ylim0):
                        y2 = y + 0.2*dy
                    else:
                        y2 = y - 0.2*dy
                        annotation = self.axes.annotate(nline, xy=(wline,y1), xytext=(wline,y2),
                                                        color='purple', alpha = 0.4, arrowsprops = dict(edgecolor='purple',
                                                                                                        facecolor = 'y',arrowstyle='-',
                                                                                                        alpha=0.4,
                                                                                                        connectionstyle = 'angle, angleA=0,angleB=90,rad=10'),
                                                        rotation = 90, fontstyle='italic',fontproperties=font, visible = True)
                        annotation.draggable()
                        self.annotations.append(annotation)
            
            self.draw_idle()
        
        

        
class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self,canvas,parent):
        # Select only a few buttons
        #self.toolitems = [t for t in NavigationToolbar2QT.toolitems if
        #                  t[0] in ('Home', 'Pan', 'Zoom', 'Save')]
        self.iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"icons")

        #icon = QIcon(self.iconDir+'exit.png')
        #print('path of icons is ', icon)
        self.toolitems = [
            ('Home','Go back to original limits','home','home'),
            ('Pan','Pan figure','move','pan'),
            ('Zoom','Zoom in','zoom_to_rect','zoom'),
        ]
        self.parent = parent
        super().__init__(canvas,parent)

