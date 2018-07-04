from astropy.io import fits
import numpy as np

def saveAnalysis(self):
    """ save wc,fc,ec, redshift, and other parameters for each galaxy in a multi-fits.
    Each extension contains a 2D array:  [3,n] with the w,f,e clipped
    and a header with: redshift, lines tuples (w,intensity,..), xlimits, ylimits
    
    """
    print ("Saving analysis ...")
    hdr = fits.Header()
    hdr['ngal']=self.ngal
    hdu = fits.PrimaryHDU(header=hdr)
    hdulist = fits.HDUList([hdu])
    
    for gal in self.galaxies:
        data = np.array([gal.w,gal.f,gal.e,gal.c])
        hdr = fits.Header()
        hdr['Redshift'] = gal.z
        hdr['zError'] = gal.dz
        hdr['Template'] = gal.zTemplate
        hdr['SpecType'] = gal.spectype
        hdr['Quality'] = gal.quality
        hdr['xlim1'] = gal.xlim1
        hdr['xlim2'] = gal.xlim2
        hdr['ylim1'] = gal.ylim1
        hdr['ylim2'] = gal.ylim2
        # Add lines
        for i in range(len(gal.lines)):
            line = gal.lines.keys()[i]
            gl = gal.lines[line]
            strline = "%.2e %.1f %.2e %.2f %.2f %.2e %.2f" % (gl[0],gl[1],gl[2],gl[3],gl[4],gl[5],gl[6])
            hdr['hierarch '+line] = strline
            
        hdu = fits.ImageHDU(data=data, header=hdr,name='SCI')
        hdulist.append(hdu)

    hdulist.writeto('showspectra.fits',clobber=True)  # overwrite


def recoverAnalysis(self):
    """ Recover previously saved data """
    hdulist = fits.open(self.dirname+"/showspectra.fits")
    igal = -1
    LinesNames = self.Lines.keys()
    for hdu in hdulist:
        hdr  = hdu.header
        # Read through Science data
        if igal >= 0:
            data = hdu.data
            gal = self.galaxies[igal]
            gal.w = data[0]
            gal.f = data[1]
            gal.e = data[2]
            c = data[3]
            gal.c = c.astype(bool) # casting to boolean
            gal.z  = hdr['Redshift']
            gal.dz  = hdr['zError']
            gal.zTemplate = hdr['Template']
            gal.quality = hdr['Quality']
            if gal.spectype != 'sky':
                gal.spectype = hdr['SpecType']
                gal.xlim1  = hdr['xlim1']
                gal.xlim2  = hdr['xlim2']
                gal.ylim1  = hdr['ylim1']
                gal.ylim2  = hdr['ylim2']
                for key in hdr.keys():
                    if key in LinesNames:
                        values = hdr[key]                        
                        #  arguments are [cont,center,amplitude,sigma,fwhm,flux,EW]
                        gal.lines[key] = [float(x) for x in values.split(" ")]
            gal.clip()
        # Read Primary header
        else:
            self.ngal = hdr['ngal']
        igal += 1
            
    hdulist.close()

