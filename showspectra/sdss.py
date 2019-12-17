from showspectra.spectra import Sky, Galaxy

def getSky(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    from astropy.io import fits
    with fits.open(file) as hdu:
        data  = hdu['SKY'].data
        logw = hdu['LOGWAVE'].data

    wave = 10**logw
    return Sky(wave, data)


def getGalaxies(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    from astropy.io import fits

    with fits.open(file) as hdu:
        header = hdu['PRIMARY'].header
        flux = hdu['FLUX'].data
        err  = hdu['ERR'].data
        logw = hdu['LOGWAVE'].data
        
    wave = 10**logw
    nf = len(flux[0,:])
    galaxies = [Galaxy(wave, flux[:,i]) for i in range(nf)]

    for i in range(len(galaxies)):
        galaxies[i].aperture = i + 1
        ra = header['HIERARCH FIBER_{:03d}_RA'.format(i)]
        dec = header['HIERARCH FIBER_{:03d}_DEC'.format(i)]
        galaxies[i].ra = ra
        galaxies[i].dec = dec
        galaxies[i].fiber = header['HIERARCH FIBER_{:03d}_ID'.format(i)]
        galaxies[i].spectype = '?'
        galaxies[i].e = err[:,i]
        
    return galaxies

