from astropy.io.fits import getdata
import numpy as np
from showspectra.spectra import Sky, Galaxy
import re


def getSky(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    data, hdr = getdata(file, 0, header=True)
    crpix = hdr['crpix1'] - 1
    cdelt = hdr['cd1_1']
    crval = hdr['crval1']
    n = hdr['naxis1']
    wave = crval + (np.arange(n) - crpix) * cdelt
    return Sky(wave, data)


def getGalaxies(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    data, hdr = getdata(file, 0, header=True)
    crpix = hdr['crpix1'] - 1
    cdelt = hdr['cd1_1']
    crval = hdr['crval1']
    wave = crval + (np.arange(len(data[0])) - crpix) * cdelt
    galaxies = [Galaxy(wave, data[i]) for i in range(len(data))]

    for i in range(len(galaxies)):
        apstring = hdr['apid' + str(i + 1)]
        galaxies[i].aperture = i + 1
        pattern = '\(([:\d\.]+) ([\+-]*[:\d\.]+)\)'
        m = re.search(pattern, apstring)
        ra = m.group(1)
        dec = m.group(2)
        galaxies[i].ra = ra
        galaxies[i].dec = dec
        m = re.search('\((\d+)\)', apstring)
        galaxies[i].fiber = m.group(1)
        s = re.search('sky', apstring)
        if s:
            galaxies[i].spectype = 'sky'
    return galaxies


def getErrors(self, file):
    # Read FITS file and define the list of objects (galaxy spectra)
    data, hdr = getdata(file, 0, header=True)
    for i in range(len(data)):
        self.galaxies[i].e = data[i]
        # Update clipped errors only if unclipped data (from recovery)
        if len(self.galaxies[i].wc) == len(self.galaxies[i].ec):
            self.galaxies[i].ec = data[i]
