from astropy.io import fits
import numpy as np


def saveAnalysis(self):
    """Save wc,fc,ec, redshift, and other parameters for each galaxy in a multi-fits.
    Each extension contains a 2D array:  [3,n] with the w,f,e clipped
    and a header with: redshift, lines tuples (w,intensity,..), xlimits, ylimits.
    """
    print("Saving analysis ...")
    hdr = fits.Header()
    hdr['ngal'] = self.ngal
    hdu = fits.PrimaryHDU(header=hdr)
    hdulist = fits.HDUList([hdu])

    for gal in self.galaxies:
        data = np.array([gal.w, gal.f, gal.e, gal.c])
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
            strline = "%.2e %.1f %.2e %.2f %.2f %.2e %.2f" %\
                      (gl[0], gl[1], gl[2], gl[3], gl[4], gl[5], gl[6])
            hdr['hierarch ' + line] = strline
        hdu = fits.ImageHDU(data=data, header=hdr, name='SCI')
        hdulist.append(hdu)

    hdulist.writeto('showspectra.fits', overwrite=True)


def recoverAnalysis(self):
    """Recover previously saved data."""
    hdulist = fits.open(self.dirname + "/showspectra.fits")
    igal = -1
    LinesNames = self.Lines.keys()
    for hdu in hdulist:
        hdr = hdu.header
        # Read through Science data
        if igal >= 0:
            data = hdu.data
            gal = self.galaxies[igal]
            gal.w = data[0]
            gal.f = data[1]
            gal.e = data[2]
            c = data[3]
            gal.c = c.astype(bool)  # casting to boolean
            gal.z = hdr['Redshift']
            gal.dz = hdr['zError']
            gal.zTemplate = hdr['Template']
            gal.quality = hdr['Quality']
            if gal.spectype != 'sky':
                gal.spectype = hdr['SpecType']
                gal.xlim1 = hdr['xlim1']
                gal.xlim2 = hdr['xlim2']
                gal.ylim1 = hdr['ylim1']
                gal.ylim2 = hdr['ylim2']
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


import json
class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

def exportAnalysis(galaxies, ngal, dirname, name=None):
    """Export results of analysis."""
    import json, io  
    from collections import OrderedDict
    data = OrderedDict([
            ('ngal', ngal),
            ('ngalaxies', len(galaxies)),
            ('waveUnit', 'micrometer'),
            ('fluxUnit', 'Jy')
            ])
    for i, galaxy in enumerate(galaxies):
        # Define list of lines
        lines = {}
        for line in galaxy.lines.copy():
            l = galaxy.lines[line]
            l.computeAll()
            lines[line] = {
                    'w1': l.w1,
                    'w2': l.w2,
                    'z': l.z,
                    'location': l.location,
                    'scale': l.scale,
                    'amplitude': l.amplitude,
                    'intercept': l.intercept,
                    'slope': l.slope,
                    'FWHM': l.FWHM,
                    'EW': l.EW,
                    'flux': l.flux
                    }
        # Find masked regions
        nmask = ~galaxy.c
        if np.sum(nmask) > 0:
            c = galaxy.c.astype(int)
            c = np.append(c, 1)
            dc = c[1:] - c[:-1]
            istart = np.where(dc == -1)
            if c[0] == 0:
                istart = np.append(0, istart)
            iend = np.where(dc == 1)
            istart = np.ravel(istart)
            iend = np.ravel(iend)
            masked = [(i,j) for (i,j) in zip(istart,iend)]
        else:
            masked = []
        # Define data for the single galaxy
        data[i] = {
                'z': galaxy.z,
                'dz': galaxy.dz,
                'ra': galaxy.ra,
                'dec': galaxy.dec,
                'quality': galaxy.quality,
                'spectype': galaxy.spectype,
                'template': galaxy.zTemplate,
                'aperture': galaxy.aperture,
                'fiber': galaxy.fiber,
                'lines': lines,
                'masked': masked
                }
        data.move_to_end(i, last=True) # Move element to the end
    with io.open(dirname+'/showspectra.json', mode='w') as f:
            str_= json.dumps(data, indent=2, separators=(',',': '),
                             ensure_ascii=False, cls=MyEncoder)
            f.write(str_)

def importAnalysis(file, galaxies):
    """Import results from previous analysis."""
    import json
    from collections import OrderedDict
    from showspectra.spectra import Line
    
    with open(file) as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    
    ngal = data['ngal']
    ngalaxies = data['ngalaxies']
    print('total galaxies: ', ngalaxies)
    for key in range(ngalaxies):
        # print('galaxy no: ', key)
        g = galaxies[key]
        d = data[str(key)]
        g.z = d['z']
        g.dz = d['dz']
        g.quality = d['quality']
        g.spectype = d['spectype']
        try:
            g.zTemplate = d['template']
        except:
            g.zTemplate = None
        masks = d['masked']
        for mask in masks:
            g.c[mask[0]:mask[1]] = 0
        lines = d['lines']
        for line in lines.copy():
            l = lines[line]
            g.lines[line] = Line(l['w1'], l['w2'], l['intercept'], l['slope'],
                    l['location'], l['scale'], l['amplitude'], l['z'])
            
    return ngal, ngalaxies, galaxies