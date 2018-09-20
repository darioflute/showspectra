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
        self.z = 0.  # initialize redshift [from xcorr]
        self.dz = 0.  # initialize redshift error
        self.ze = None  # z from emission lines
        self.za = None  # z from absorption lines
        self.zTemplate = None  # template with best cross-correlation
        self.spectype = '?'  # type (?,star,galaxy,broadAGN,sky)
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
        dw = np.nanmax(w) - np.nanmin(w)
        self.xlim1 = np.nanmin(w) - dw / 30.
        self.xlim2 = np.nanmax(w) + dw / 30.
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


class Line(object):
    """Line definition."""
    def __init__(self, w1, w2, center, intcpt, intcptErr, slope, slopeErr,
                 loc, locErr, scale, scaleErr, amp, ampErr):
        self.w1 = w1
        self.w2 = w2
        self.center = center
        self.location = loc
        self.eLocation = locErr
        self.scale = scale
        self.eScale = scaleErr
        self.amplitude = amp
        self.eAmplitude = ampErr
        self.intercept = intcpt
        self.eIntercept = intcptErr
        self.slope = slope
        self.eSlope = slopeErr
        self.computeAll()

    def computeAll(self):
        self.z = (self.location - self.center) / self.center  # z from the line
        self.eZ = self.eLocation / self.center
        # print('z is ', self.z)
        w0 = (self.intercept + self.slope * self.location)  # observed
        # eW0 = self.eIntercept + self.location * self.eSlope
        # We consider only error in intercept since the continuum is evaluated around the line
        eW0 = self.eIntercept
        self.flux = np.sqrt(2 * np.pi) * np.abs(self.amplitude) * self.scale
        self.eFlux = self.flux * (self.eAmplitude/self.amplitude + self.eScale/self.scale) 
        self.EW = self.flux / w0 / (1. + self.z)  # At rest
        self.eEW = self.EW * (self.eFlux / self.flux + eW0 / w0)
        self.FWHM = 2.35483 * self.scale / (1 + self.z)  # At rest
        self.eFWHM = self.FWHM * self.eScale / self.scale

class Spectrum(object):
    """Analysis of spectrum."""

    def __init__(self):
        self.z = 0.  # initialize redshift [from xcorr]
        self.dz = 0.  # initialize redshift error
        self.ze = None  # z from emission lines
        self.za = None  # z from absorption lines
        self.zTemplate = None  # template with best cross-correlation
        self.spectype = '?'  # type (?,star,galaxy,broadAGN,sky)
        self.quality = '?'  # quality ( OK,guess,?)
        self.lines = {}
        # Identifications
        self.ra = None
        self.dec = None
        self.fiber = None
        self.aperture = None

