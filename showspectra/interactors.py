import numpy as np
from PyQt5.QtCore import pyqtSignal, QObject
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.patches import Polygon


class SegmentsSelector(QObject):
    """This class allows one to mark two segments to define the continuum around a line."""
    
    def __init__(self, ax, fig, callback, color='#7ec0ee', zD = True):
        super().__init__()

        self.x = []
        self.y = []
        self.line1 = None
        self.line2 = None
        self.color = color
        self.fig = fig
        self.ax = ax
        self.callback = callback
        self.zeroDeg = zD

        self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
        self.__ID2 = self.fig.canvas.mpl_connect('button_press_event', self.__button_press_callback)
        self.__ID3 = self.fig.canvas.mpl_connect('button_release_event', self.__button_release_callback)

    def __motion_notify_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            if (event.button == None or event.button == 1):
                if self.line1 != None: # Move line around
                    if self.line2 == None:
                        if self.zeroDeg: self.y[0]=y
                        self.line1.set_data([self.x[0], x],
                                            [self.y[0], y])
                    else:
                        if self.zeroDeg:
                            self.y=[y,y,y,y]
                            self.line1.set_data([self.x[0], self.x[1]],
                                                [self.y[0], self.y[1]])
                        self.line2.set_data([self.x[2], x],
                                            [self.y[2], y])
                    self.fig.canvas.draw_idle()


    def __button_release_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            if event.button == 1:
                if self.line2 == None:  # Segment 1 completed
                    self.x.append(x)
                    self.y.append(y)
                    if self.zeroDeg:
                        self.y[-2]=self.y[-1]
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    self.fig.canvas.draw_idle()
                    self.fig.canvas.mpl_disconnect(self.__ID1)
                    self.line2 = 'start'
                else:
                    self.x.append(x)
                    self.y.append(y)
                    if self.zeroDeg:
                        self.y[-1]=self.y[-2]
                        m = 0.
                    else:
                        # Adjust to the same slope between first and last point
                        m = (self.y[3]-self.y[0])/(self.x[3]-self.x[0])
                    for i in range(4):
                        self.y[i] = self.y[0]+m*(self.x[i]-self.x[0])
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    self.line2.set_data([self.x[2], self.x[3]],
                                        [self.y[2], self.y[3]])
                    self.fig.canvas.draw_idle()
                    self.xy = [(i,j) for (i,j) in zip(self.x,self.y)]
                    # Disconnect
                    self.fig.canvas.mpl_disconnect(self.__ID1) 
                    self.fig.canvas.mpl_disconnect(self.__ID2) 
                    self.fig.canvas.mpl_disconnect(self.__ID3) 
                    # Callback function, pass the vertices
                    self.callback(self.xy)
                    # Remove lines
                    self.remove()

    def __button_press_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1:
                if self.line2 == None:
                    if self.line1 == None:  # If you press the left button, single click
                        self.line1 = Line2D([x, x],
                                            [y, y],
                                            marker='o',
                                            color=self.color)
                        self.start_point = [x,y]
                        self.previous_point =  self.start_point
                        self.x=[x]
                        self.y=[y]
                        self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
                        ax.add_line(self.line1)
                        # add a segment
                        self.fig.canvas.draw_idle()
                else:
                    if self.zeroDeg:
                        self.y = [y,y]
                        self.line1.set_data([self.x[0], self.x[1]],
                                            [self.y[0], self.y[1]])
                    self.line2 = Line2D([x, x],
                                        [y, y],
                                        marker='o',
                                        color=self.color)
                    self.start_point = [x,y]
                    self.previous_point =  self.start_point
                    self.x.append(x)
                    self.y.append(y)
                    self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
                        
                    ax.add_line(self.line2)
                    self.fig.canvas.draw()

    def remove(self):
        """ Remove lines from plot """
        try:
            self.line1.remove()
            self.line2.remove()
        except:
            print('no lines to remove')


class SegmentsInteractor(QObject):
    """
    A segmented continuum interactor.
    """
    
    showverts = True
    epsilon = 10  # max pixel distance to count as a vertex hit
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, verts, zeroDeg=True):
        super().__init__()


        self.ax = ax
        self.type = 'Continuum'
        #        color = 'skyblue'
        color = '#7ec0ee'

        self.zeroDeg = zeroDeg
        x, y = zip(*verts)
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.computeSlope()
        self.line1 = Line2D(x[:2],y[:2],color=color,linewidth=2, animated = True)
        self.line2 = Line2D(x[2:],y[2:],color=color,linewidth=2, animated = True)

        self.canvas = ax.figure.canvas
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., 
                           markerfacecolor=color, animated=True)                
        self.ax.add_line(self.line1)
        self.ax.add_line(self.line2)
        self.ax.add_line(self.line)

        self.cid = self.line1.add_callback(self.si_changed)
        self._ind = None  # the active vert
        self.connect()

    def computeSlope(self):

        xg,yg = zip(*self.xy)
        xg = np.array(xg); yg = np.array(yg)
        if self.zeroDeg:
            self.slope = 0
        else:
            self.slope = (yg[3]-yg[0])/(xg[3]-xg[0])
        self.intcpt = yg[0]-self.slope*xg[0]


    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', 
                                                   self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', 
                                                  self.motion_notify_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.canvas.mpl_disconnect(self.cid_key)
        try:
            self.line1.remove()
        except:
            print('no line 1')
        try:
            self.line2.remove()
        except:
            print('no line 2')
        try:
            self.line.remove()
        except:
            print('no markers')
        self.canvas.draw_idle()
        # self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.line1)
        self.ax.draw_artist(self.line2)
        self.ax.draw_artist(self.line)

    def si_changed(self, line1):
        'this method is called whenever the line1 object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, line1)
        self.line.set_visible(vis)  

    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'

        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform(self.xy)
        x, y = zip(*xy)
        x = np.array(x); y = np.array(y)
        d = np.hypot(x - event.x, y - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        # print('distance is ',d[ind])
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return

        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            self.mySignal.emit('line deleted')
        self.canvas.draw_idle()

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None


    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x_, y_ = event.xdata, event.ydata
        # Rebuild line collection
        x,y = zip(*self.xy)
        x = np.asarray(x)
        y = np.asarray(y)
        # Update point
        y[self._ind] = y_
        if self._ind > 0:
            if x_ < x[self._ind-1]:
                x[self._ind] = x[self._ind-1]
            else:
                x[self._ind] = x_
        if self._ind < 3:
            if x_ > x[self._ind+1]:
                x[self._ind] = x[self._ind+1]
            else:
                x[self._ind] = x_
        if self.zeroDeg:
            m = 0
        else:
            if self._ind < 2:
                m = (y[3]-y[self._ind])/(x[3]-x[self._ind])
            else:
                m = (y[self._ind]-y[0])/(x[self._ind]-x[0])
        for i in range(4):
            y[i] = y[self._ind]+m*(x[i]-x[self._ind])
            self.xy[i] = (x[i],y[i])
        # Update segments and markers
        self.updateLinesMarkers()
        #self.canvas.restore_region(self.background)
        #self.ax.draw_artist(self.line1)
        #self.ax.draw_artist(self.line2)
        #self.ax.draw_artist(self.line)
        #self.canvas.update()
        #self.canvas.flush_events()
        self.canvas.draw_idle()
        # Notify callback
        self.modSignal.emit('continuum guess modified')

    def updateLinesMarkers(self):
        self.line1.set_data(zip(*self.xy[:2]))
        self.line2.set_data(zip(*self.xy[2:]))
        self.line.set_data(zip(*self.xy))
        # Update parameters of the continuum
        self.computeSlope()


class LineInteractor(QObject):
    """
    A Gaussian line interactor.
    """
    
    showverts = True
    epsilon = 10  # max pixel distance to count as a vertex hit
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, c0, cs, x0, A, fwhm):
        super().__init__()

        self.ax = ax
        self.canvas = ax.figure.canvas
        self.c0 = c0  # Value of continuum at origin
        self.cs = cs  # Slope of the continuum
        self.type = 'Line'
        color = '#7ec0ee'
        self.x0 = x0
        self.A = A
        self.fwhm = fwhm
        self.computeMarkers()
        self.computeGaussian()
        self.gauss = Polygon(list(self.verts), animated=True, fill=False, closed=False, color=color)
        self.ax.add_patch(self.gauss)
        x, y = zip(*self.xy)
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor=color, animated=True)                
        self.ax.add_line(self.line)
        # Callback for changes
        self.cid = self.gauss.add_callback(self.si_changed)
        self._ind = None  # the active vert
        self.connect()

    def computeMarkers(self):
        'Compute position of markers.'
        x = self.x0 + 0.5 * self.fwhm * np.array([-1, 0, 1])
        y = self.c0 + x * self.cs + self.A * np.array([0.5, 1., 0.5])
        self.xy = [(i,j) for (i,j) in zip(x,y)]

    def computeGaussian(self):
        'Compute the Gaussian polygon from the position of the markers.'
        self.sigma = self.fwhm/(2 * np.sqrt(2 * np.log(2)))
        # Create an array of x values and compute the value of the Gaussian on it
        x = np.linspace(self.x0 - self.fwhm, self.x0 + self.fwhm, 30)
        dx = (x - self.x0) / self.sigma / np.sqrt(2.)
        y = self.c0 + x * self.cs + self.A * np.exp(-dx * dx)
        self.verts = [(x_,y_) for x_,y_ in zip(x,y)]
        return 

    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.canvas.mpl_disconnect(self.cid_key)
        try:
            self.line.remove()
        except:
            print('no markers')
        try:
            self.gauss.remove()
        except:
            print('no line')
        self.canvas.draw_idle()
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.gauss)
        self.ax.draw_artist(self.line)

    def si_changed(self, artist):
        'this method is called whenever the artist object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, artist)
        self.line.set_visible(vis)  

    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'
        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform(self.xy)
        x, y = zip(*xy)
        x = np.array(x); y = np.array(y)
        d = np.hypot(x - event.x, y - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
        # print('distance is ',d[ind])
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            self.mySignal.emit('line deleted')
        self.canvas.draw_idle()

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None


    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x_, y_ = event.xdata, event.ydata
        # Rebuild line collection
        x, y = zip(*self.xy)
        x = np.asarray(x)
        y = np.asarray(y)
        # Update markers and Gaussian parameters
        if self._ind == 0:
            if x_ > x[1]:
                pass
            else:
                #x[self._ind] = x_
                #x[2] = x[1] + (x[1]-x_)  # Symmetrical change
                #dx = x[1] - x_
                self.fwhm = 2 * (x[1] - x_)
        elif self._ind == 1:
            dx = x_ - x[1]
            #x[0] += dx
            #x[1] += dx
            #x[2] += dx
            self.x0 += dx
            dy = y_ - y[1]
            if (self.A > 0) & (dy < -self.A):  # Emission line 
                pass
            elif (self.A < 0) & (dy > -self.A):  # Absorption line
                pass
            else:
                #y[1] += dy/2.
                #y[0] += dy
                #y[2] += dy/2.
                self.A += dy
        elif self._ind == 2:
            if x_ < x[1]:
                pass
            else:
                #x[2] = x_
                #x[0] = x[1] - (x[2]-x[1])  # Symmetrical change
                self.fwhm = 2 * (x_ - x[1])
        # self.xy = [(i,j) for (i,j) in zip(x,y)]
        # self.fwhm = x[2] - x[0]
        self.updateCurves()
        # Update markers and Gaussian
        # self.line.set_data(zip(*self.xy))
        # self.computeGaussian()
        # self.gauss.xy = self.verts
        #self.canvas.restore_region(self.background)
        #self.ax.draw_artist(self.gauss)
        #self.ax.draw_artist(self.line)
        #self.canvas.update()
        #self.canvas.flush_events()
        # self.canvas.draw_idle()
        # Notify callback
        self.modSignal.emit('line guess modified')
        
    def updateCurves(self):
        self.computeGaussian()  
        self.computeMarkers()
        self.ax.draw_artist(self.gauss)
        self.line.set_data(zip(*self.xy))
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.line)
        self.gauss.xy = self.verts
        self.ax.draw_artist(self.gauss)
        self.canvas.update()
        self.canvas.flush_events()
        #self.canvas.draw_idle()