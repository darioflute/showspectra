from astropy.io import fits
import numpy as np
import json


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
    import io
    from collections import OrderedDict
    data = OrderedDict([('ngal', ngal),
                        ('ngalaxies', len(galaxies)),
                        ('waveUnit', 'Angstrom'),
                        ('fluxUnit', 'W/m2/Hz')
                        ])
    for i, galaxy in enumerate(galaxies):
        # Define list of lines
        lines = {}
        # Redshifts
        zem = []
        zabs = []
        for line in galaxy.lines.copy():
            li = galaxy.lines[line]
            li.computeAll()
            lines[line] = {'w1': li.w1,
                           'w2': li.w2,
                           'z': [li.z, li.eZ],
                           'center': li.center,
                           'location': [li.location, li.eLocation],
                           'scale': [li.scale, li.eScale],
                           'amplitude': [li.amplitude, li.eAmplitude],
                           'intercept': [li.intercept, li.eIntercept],
                           'slope': [li.slope, li.eSlope],
                           'FWHM': [li.FWHM, li.eFWHM],
                           'EW': [li.EW, li.eEW],
                           'flux': [li.flux, li.eFlux]
                           }
            if li.amplitude > 0:
                zem.append(li.z)
            else:
                zabs.append(li.z)
        # Find masked regions
        nmask = ~galaxy.c
        if np.sum(nmask) > 0:
            c = galaxy.c.astype(int)
            c = np.append(c, 1)
            dc = c[1:] - c[:-1]
            istart = np.where(dc == -1)
            if c[0] == 0:
                istart = np.append(-1, istart)
            iend = np.where(dc == 1)
            istart = np.ravel(istart) + 1
            iend = np.ravel(iend) + 1
            masked = [(i, j) for (i, j) in zip(istart, iend)]
        else:
            masked = []
        # Compute emission and absorption spectra
        if len(zem) > 0:
            zem = np.array(zem)
            ze = np.nanmean(zem)
        else:
            ze = 0.
        if len(zabs) > 0:
            zabs = np.array(zabs)
            za = np.nanmean(zabs)
        else:
            za = 0.
        # Define data for the single galaxy
        data[i] = {'z': galaxy.z,
                   'dz': galaxy.dz,
                   'ze': ze,
                   'za': za,
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
        data.move_to_end(i, last=True)  # Move element to the end
    with io.open(dirname + '/showspectra.json', mode='w') as f:
            str_ = json.dumps(data, indent=2, separators=(',', ':'),
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
    # print('total galaxies: ', ngalaxies)
    for key in range(ngalaxies):
        g = galaxies[key]
        # print('galaxy no: ', key, np.size(g.c))
        d = data[str(key)]
        g.z = d['z']
        g.dz = d['dz']
        try:
            g.za = d['za']
        except BaseException:
            g.za = None
        try:
            g.ze = d['ze']
        except BaseException:
            g.za = None
        g.quality = d['quality']
        g.spectype = d['spectype']
        try:
            g.zTemplate = d['template']
        except BaseException:
            g.zTemplate = None
        masks = d['masked']
        for mask in masks:
            if mask[1] == np.size(g.c):
                g.c[mask[0]:] = 0
            else:
                g.c[mask[0]:mask[1]] = 0
        lines = d['lines']
        for line in lines.copy():
            li = lines[line]
            g.lines[line] = Line(li['w1'], li['w2'], li['center'], 
                                 li['intercept'][0], li['intercept'][1],
                                 li['slope'][0], li['slope'][1],
                                 li['location'][0], li['location'][1],
                                 li['scale'][0], li['scale'][1],
                                 li['amplitude'][0], li['amplitude'][1])
    return ngal, ngalaxies, galaxies

def getSpectra(file):
    """Import results from previous analysis."""
    import json
    from collections import OrderedDict
    from showspectra.spectra import Line, Spectrum
    with open(file) as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
        
    n = data['ngalaxies']
    spectra = []
    for i in range(n):
        g = Spectrum()
        ii = str(i)
        d = data[ii]
        g.z = d['z']
        g.dz = d['dz']
        g.ra = d['ra']
        g.dec = d['dec']
        try:
            g.za = d['za']
        except BaseException:
            g.za = None
        try:
            g.ze = d['ze']
        except BaseException:
            g.za = None
        try:
            g.aperture = d['aperture']
        except BaseException:
            g.aperture = None
        try:
            g.fiber = d['fiber']
        except BaseException:
            g.fiber = None
        g.quality = d['quality']
        g.spectype = d['spectype']
        try:
            g.zTemplate = d['template']
        except BaseException:
            g.zTemplate = None
        lines = d['lines']
        for line in lines.copy():
            li = lines[line]
            g.lines[line] = Line(li['w1'], li['w2'], li['center'], 
                                 li['intercept'][0], li['intercept'][1],
                                 li['slope'][0], li['slope'][1],
                                 li['location'][0], li['location'][1],
                                 li['scale'][0], li['scale'][1],
                                 li['amplitude'][0], li['amplitude'][1])
        spectra.append(g)
    return spectra