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

    # Incorrect: Rebecca says they are all galaxies (two slits have more than one galaxy)
    # HIERARCH ESO INS SLIT18 ID = 22 / ID of SLIT i
    # HIERARCH ESO INS SLIT18 OBJ RA = 10.366204 / RA (J2000) of object
    # HIERARCH ESO INS SLIT18 OBJ DEC = -9.160961 / DEC (J2000) of object
    eso = 'HIERARCH ESO INS'
    n = len(galaxies)
    fiber = np.zeros(n)
    ra = np.zeros(n)
    dec = np.zeros(n)
    x = np.zeros(n)
    galsky = np.zeros(n)
    i = 0
    while True:
        fiber[i] = hdr[eso + ' SLIT' + str(i + 1) + ' ID']
        ra[i] = hdr[eso + ' SLIT' + str(i + 1) + ' OBJ RA']
        dec[i] = hdr[eso + ' SLIT' + str(i + 1) + ' OBJ DEC']
        x[i] = hdr[eso + ' SLIT' + str(i + 1) + ' X']
        try:
            hdr[eso + ' SLIT' + str(i + 2) + ' X']
        except BaseException:
            break
        i += 1
    j = 0
    i += 1
    while True:
        fiber[i] = hdr[eso + ' REF' + str(j + 1) + ' ID']
        ra[i] = hdr[eso + ' REF' + str(j + 1) + ' OBJ RA']
        dec[i] = hdr[eso + ' REF' + str(j + 1) + ' OBJ DEC']
        x[i] = hdr[eso + ' REF' + str(j + 1) + ' X']
        galsky[i] = 1
        try:
            hdr[eso + ' REF' + str(j + 2) + ' X']
        except BaseException:
            break
        j += 1
        i += 1

    # order according to x
    idx = np.argsort(x)
    fiber = fiber[idx]
    ra = ra[idx]
    dec = dec[idx]
    x = x[idx]
    galsky = galsky[idx]

    for i in range(len(galaxies)):
        galaxies[i].fiber = str(fiber[i])
        galaxies[i].ra = str(ra[i])
        galaxies[i].dec = str(dec[i])
        if galsky[i]:
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
