from astropy.io.fits import getdata
import numpy as np
from showspectra.spectra import Sky, Galaxy


def getSky(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    data, hdr = getdata(file, 0, header=True)
    crpix = hdr['crpix1'] - 1
    cdelt = hdr['cd1_1']
    crval = hdr['crval1']
    n = hdr['naxis1']
    wave = crval + (np.arange(n) - crpix) * cdelt
    skyflux = data[0]
    return Sky(wave, skyflux)


def getGalaxies(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    data, hdr = getdata(file, 0, header=True)
    crpix = hdr['crpix1'] - 1
    cdelt = hdr['cd1_1']
    crval = hdr['crval1']
    wave = crval + (np.arange(len(data[0])) - crpix) * cdelt
    galaxies = [Galaxy(wave, data[i]) for i in range(len(data))]

    # Source identification added to the header
    # HIERARCH SLIT_??_ID, RA, DEC
    for i in range(len(galaxies)):
        galaxies[i].fiber = hdr['HIERARCH SLIT_{:02d}_ID'.format(i)]
        galaxies[i].ra = hdr['HIERARCH SLIT_{:02d}_RA'.format(i)]
        galaxies[i].dec = hdr['HIERARCH SLIT_{:02d}_DEC'.format(i)]

    return galaxies


def getErrors(self, file):
    # Read FITS file and define the list of objects (galaxy spectra)
    data, hdr = getdata(file, 0, header=True)
    for i in range(len(data)):
        self.galaxies[i].e = data[i]
        # Update clipped errors only if unclipped data (from recovery)
        if len(self.galaxies[i].wc) == len(self.galaxies[i].ec):
            self.galaxies[i].ec = data[i]
