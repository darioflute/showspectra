import numpy as np
from lmfit import Parameters, minimize
from scipy.ndimage.filters import median_filter
from showspectra.dialogs import selectRedshift

def MAD(data, axis=None):
    """Median of absolute deviations."""
    return np.nanmedian(np.absolute(data - np.nanmedian(data, axis)), axis)

def gaussResiduals(p, x, data=None):
    v = p.valuesdict()
    x0 = v['center']
    s = v['sigma']
    A = v['amplitude']
    model = A / (np.sqrt(2 * np.pi) * s) * np.exp(-0.5 * ((x - x0) / s)**2)
    if data is None:
        return model
    else:
        return (model - data)


def computeXcorr(wg, fg, eg, template):
    # B) compare with templates
    # interpolate the galaxy spectrum on the same logarithmic wavelength scale as the template
    # select only the common wavelength space
    wt = template.w.copy()
    ft = template.f.copy()
    ft /= np.std(ft)
    wt0 = np.log10(wt[0])
    wtd = 0.0001
    minwg = np.min(wg)
    maxwg = np.max(wg)
    # If the spectrum extends at wavelengths longer than those of the template ...
    if maxwg > np.max(wt):
        # Extend wt
        nmax = (np.log10(maxwg) - wt0) / wtd
        nmax = np.rint(nmax)
        wave = wt0 + np.arange(nmax) * wtd
        wt = np.exp(np.log(10.) * wave)
        nzero = len(wt) - len(ft)
        ft = np.append(ft, np.zeros(nzero))
    # Zero wg for wavelength lower than min(wt) or cut it ?
    # if minwg < np.min(wt):
    #    minwg = np.min(wt)
    fi = np.interp(wt, wg, fg)
    mask = (wt <= minwg) | (wt >= maxwg)
    idx = np.where(mask)
    fi[idx] = 0
    cc = np.correlate(fi, ft, mode='full')
    # Consider only the positive part
    cc = cc[len(fi) - 1:]
    # Compute the number of points used
    ncc = np.zeros(len(fi))
    for i in range(len(fi)):
        ncc[i] = (fi[i:] != 0).sum()
    # Look for the most significant maximum
    maxc = []
    maxs = []
    ccsigma = MAD(cc)
    cc1 = cc[2:-2] - cc[1:-3]
    cc2 = cc[2:-2] - cc[3:-1]
    idmax, = np.where((cc1 > 0) & (cc2 > 0))
    for im in idmax:
        i = im + 2
        n1 = i - 50
        n2 = i + 50
        if n1 < 0:
            n1 = 0
        if n2 > len(cc):
            n2 = len(cc)
        ind = list(range(n1, i - 2)) + list(range(i + 2, n2))
        # val = (cc[i]-np.median(cc[n1:n2]))/np.std(cc)  # this works but misses some cases
        # val = (cc[i]-np.median(cc[ind]))/self.mad(cc[ind])  # this works pretty well
        # val = np.median(cc[i]-cc[n1:n2])/self.mad(cc[n1:n2])  # this works pretty well
        # val = (cc[i]-np.median(cc[n1:n2]))/np.std(cc[n1:n2])  # this does not work !
        # val = (cc[i]-np.median(cc[n1:n2]))/ccsigma  # this works pretty well
        # val = cc[i]/self.mad(cc[n1:n2])                         # Does not work
        val = (cc[i] - np.median(cc[ind])) / ccsigma
        if (val > 3):
            # print " z {0:.5f}".format(np.power(10,0.0001*i))," SNR {0:.3f} ".format(val)
            maxc.append(i)
            maxs.append(val)
    maxc = np.array(maxc)
    maxs = np.array(maxs)
    lminwt = np.log10(np.min(wt))
    lmaxwt = np.log10(np.max(wt))
    if len(maxs) > 0:
        idx = np.argsort(maxs)
        nm = maxc[idx[-2:]]
        snr = maxs[idx[-2:]]
    else:
        return ([0], [0], [0])
    dl = 0.0001  # cd1_1 of SDSS templates
    x = np.arange(len(cc))
    ztot = []
    sztot = []
    for i in range(len(nm)):
        z = np.power(10, dl * nm[i]) - 1.0
        # Number of overlapping bins at maximum correlation
        lminwg = np.log10(minwg) - nm[i] * dl
        lmaxwg = np.log10(maxwg) - nm[i] * dl
        lminw = np.max([lminwt, lminwg])
        lmaxw = np.min([lmaxwt, lmaxwg])
        N = (lmaxw - lminw) / dl
        # Fitting the peak with a Gaussian
        y = cc / N  # Normalize cross-correlation to number of points used
        fit_params = Parameters()
        fit_params.add('center', nm[i])
        fit_params.add('sigma', 10.)
        fit_params.add('amplitude', 1000.)
        out = minimize(gaussResiduals, fit_params, args=(x, ), kws={'data': y}, method='Nelder')
        S = out.params['sigma'].value            # sigma of the fitting Gaussian
        A = y[nm[i]]
        sigmaz = dl * (1 + z) * np.log(10.) * S / np.sqrt(N) * np.sqrt(1 / A**2 - 1.)
        # print " z {0:.5f}".format(z), "{0:.5f}".format(sigmaz),
        # " SNR {0:.3f} ".format(snr[i]), " N ",N
        ztot.append(z)
        sztot.append(sigmaz)
    return (ztot, sztot, snr)


def cross_correlation(self):
    """Cross correlation between galaxy spectrum and templates."""
    # A) subtract continuum from galaxy spectrum
    self.gal = self.galaxies[self.ngal]
    wg = self.gal.wc.copy()
    fg = self.gal.fc.copy()
    eg = self.gal.ec.copy()
    # Does not work with VIMOS data, but it works better with WIYN ...
    if self.telescope == 'wiyn':
        fg /= eg
    bgr = median_filter(fg, 128, mode='mirror')
    fg -= bgr
    stdfg = np.std(fg)
    fg /= stdfg
    # B) compute the max of cross-correlation
    z = []
    sz = []
    snr = []
    temps = []
    for template in self.templates.keys():
        # print template
        a, b, c = computeXcorr(wg, fg, eg, self.templates[template])
        if a != [0]:
            z.extend(a)
            sz.extend(b)
            snr.extend(c)
            if len(a) > 0:
                temps.append(template)
            if len(a) > 1:
                temps.append(template)

    z = np.array(z)
    sz = np.array(sz)
    snr = np.array(snr)
    temps = np.array(temps)

    # Print the top five x-correlations
    idx = np.argsort(snr)
    idxs = idx[-7:]
    idxs = idxs[::-1]  # reverse list
    self.zxcorr = z[idxs]
    self.szxcorr = sz[idxs]
    self.txcorr = temps[idxs]
    self.snrxcorr = snr[idxs]
    for i in range(len(idxs)):
        ii = idxs[i]
        print("Best xcorr: z {0:.5f}".format(z[ii]),
              " sz {0:.5f}".format(sz[ii]),
              " snr {0:.3f}".format(snr[ii]),
              temps[ii])
    idx = (z > 0) & (snr > 3)
    if idx.sum() > 0:
        snr0 = snr[idx]
        idx, = np.where(snr == max(snr0))
        self.gal.z = z[idx[0]]
        self.gal.dz = sz[idx[0]]
        self.gal.zTemplate = temps[idx[0]]
        self.gal.limits()
        self.showTemplate = True
        self.sp.drawSpectrum()
        # Prompt a dialogue window to show the best 5 x-corr and
        # the possibility to choose a different solution
        self.selectZ = selectRedshift(self.zxcorr, self.szxcorr, self.snrxcorr, self.txcorr)
        self.selectZ.list.currentRowChanged.connect(self.sp.updateTemplate)
        self.selectZ.exec_()
    else:
        print("No template fits the data")
