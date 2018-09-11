"""Classes defining spectra."""

import numpy as np


class Galaxy(object):
    """Galaxy with spectrum."""

    def __init__(self, wave, flux):
        self.w = wave
        self.f = flux
        self.e = np.ones(len(wave))  # initialize to same error
        self.c = np.ones(len(wave), dtype=bool)  # inititialization of clipping mask
        self.wc = wave  # initialize clipped wavelenght
        self.fc = flux  # initialize clipped spectrum
        self.ec = np.ones(len(wave))  # initialize to same error
        self.z = 0.  # initialize redshift
        self.dz = 0.  # initialize redshift error
        self.zTemplate = None  # template with best cross-correlation
        self.spectype = None  # type (?,star,galaxy,broadAGN,sky)
        self.quality = '?'  # quality ( OK,guess,?)
        # Set limits
        self.limits()
        self.lines = {}
        # Identifications
        self.ra = None
        self.dec = None
        self.fiber = None
        self.aperture = None

    # Redefine limits in wavelenght and flux
    def limits(self):
        w = self.wc / (1 + self.z)
        f = self.fc
        minf = np.nanmin(f[self.c])
        maxf = np.nanmax(f[self.c])
        df = maxf - minf
        mf = 0.5 * (maxf + minf)
        dw = np.nanmax(w)-np.nanmin(w)
        self.xlim1 = np.nanmin(w)-dw/30.
        self.xlim2 = np.nanmax(w)+dw/30.
        self.ylim1 = mf - df * 0.55
        self.ylim2 = mf + df * 0.55

    # Update clipped values
    def clip(self):
        self.wc = self.w[self.c]
        self.fc = self.f[self.c]
        self.ec = self.e[self.c]


class Sky(object):
    """Sky spectrum."""
    def __init__(self, wave, flux):
        self.w = wave
        self.f = flux
