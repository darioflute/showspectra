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
        self.quality = None  # quality ( OK,guess,?)
        self.xlim1 = min(wave)
        self.xlim2 = max(wave)
        self.ylim1 = min(flux)
        self.ylim2 = max(flux)
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
        minf = min(f)
        maxf = max(f)
        df = maxf - minf
        mf = 0.5 * (maxf + minf)
        self.xlim1 = min(w)
        self.xlim2 = max(w)
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
