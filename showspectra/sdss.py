from showspectra.spectra import Sky, Galaxy

def getSky(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    from astropy.io import fits
    import numpy as np
    print('Reading sky for ', file)
    with fits.open(file) as hdu:
        data  = hdu['SKY'].data
        print('shape ', np.shape(data))
        logw = hdu['LOGWAVE'].data

    wave = 10**logw
    print('sky defined')
    return Sky(wave, data)


def getGalaxies(file):
    """Read FITS file and define the list of objects (galaxy spectra)."""
    from astropy.io import fits

    with fits.open(file) as hdu:
        header = hdu['PRIMARY'].header
        flux = hdu['FLUX'].data * 1.e-17
        err  = hdu['ERR'].data * 1.e-17
        logw = hdu['LOGWAVE'].data
        sky = hdu['SKY'].data * 1.e-17
        
    wave = 10**logw
    nf = len(flux[0,:])
    print('there are ', nf, ' galaxies')
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
        
    
    return galaxies, Sky(wave, sky)
