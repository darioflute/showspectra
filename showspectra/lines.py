import numpy as np


def define_lines(telescope):
    import collections
    alpha = u'\u03B1'
    beta = u'\u03B2'
    gamma = u'\u03B3'
    delta = u'\u03B4'
    eps = u'\u03B5'
    #        alpha.encode('utf8')
    if telescope == 'WIYN' or telescope == 'VIMOS' or telescope == 'SDSS':
        return collections.OrderedDict([
                ('OVI 1033', ['OVI', 1033.82]),
                ('Ly-alpha 1215', ['Ly' + alpha, 1215.24]),
                ('N-V 1240', ['N-V', 1240.81]),
                ('OI 1304', ['OI', 1304.53]),
                ('Si-IV 1397', ['Si-IV', 1397.61]),
                ('C-IV 1549', ['C-IV', 1549.48]),
                ('He-II 1640', ['He-II', 1640.40]),
                ('O-III 1666', ['O-III', 1665.85]),
                ('Al-III 1857', ['Al-III', 1857.40]),
                ('C-III 1908', ['C-III', 1908.734]),
                ('C-II 2326', ['C-II', 2326.00]),
                ('Ne-IV 2439', ['Ne-IV', 2439.50]),
                ('Mg-II 2799', ['Mg-II', 2799.117]),  # SDSS
                ('Ne-V 3346', ['Ne-V', 3346.79]),
                ('Ne-VI 3426', ['Ne-VI', 3426.85]),
                ('[OII] 3728', ['[OII]', 3728.48]),
                ('[NeIII] 3869', ['[NeIII]', 3869.85]),
                ('H-8 3890', ['H-8', 3890.15]),
                ('H-eps 3971', ['H' + eps, 3971.20]),
                ('H-delta 4102', ['H' + delta, 4102.89]),
                ('H-gamma 4341', ['H' + gamma, 4341.68]),
                ('[OIII] 4363', ['[OIII]', 4363.2]),
                ('H-beta 4862', ['H' + beta, 4862.68]),
                ('[OIII] 4960', ['[OIII]', 4960.30]),
                ('[OIII] 5008', ['[OIII]', 5008.24]),
                ('[NI] 5199', ['[NI]', 5199]),
                ('[NII] 5755', ['[NII]', 5754.6]),
                ('HeI 5877', ['HeI', 5877.29]),
                ('[OI] 6302', ['[OI]', 6302.05]),
                ('[NII] 6549', ['[NII]', 6549.85]),
                ('H-alpha 6564', ['H' + alpha, 6564.61]),
                ('[NII] 6585', ['[NII]', 6585.28]),
                ('SII 6718', ['SII', 6718.29]),
                ('SII 6732', ['SII', 6732.67]),
                ('A:Ca(H) 3934', ['A:Ca(H)', 3934.78]),
                ('A:Ca(K) 3969', ['A:Ca(K)', 3969.59]),
                ('A:G-band 4300', ['A:G-band', 4300.4]),
                ('A:Mg-1 5167', ['A:Mg-1', 5167.3222]),
                ('A:Mg-2 5172', ['A:Mg-2', 5172.6847]),
                ('A:Mg-3 5183', ['A:Mg-3', 5183.6046]),
                ('A:Na 5894', ['A:Na', 5894.57]),
                ('A:H-delta 4102', ['A:H' + delta, 4102.89]),
                ('A:H-gamma 4341', ['A:H' + gamma, 4341.68]),
                ('A:H-beta 4862', ['A:H' + beta, 4862.68]),
                ('A:H-alpha 6564', ['A:H' + alpha, 6564.61])
                ])
    else:  # Higher resolution can see the [OII] doublet
        return collections.OrderedDict([
                ('OVI 1033', ['OVI', 1033.82]),
                ('Ly-alpha 1215', ['Ly' + alpha, 1215.24]),
                ('N-V 1240', ['N-V', 1240.81]),
                ('OI 1304', ['OI', 1304.53]),
                ('Si-IV 1397', ['Si-IV', 1397.61]),
                ('C-IV 1549', ['C-IV', 1549.48]),
                ('He-II 1640', ['He-II', 1640.40]),
                ('O-III 1666', ['O-III', 1665.85]),
                ('Al-III 1857', ['Al-III', 1857.40]),
                ('C-III 1908', ['C-III', 1908.734]),
                ('C-II 2326', ['C-II', 2326.00]),
                ('Ne-IV 2439', ['Ne-IV', 2439.50]),
                ('Mg-II 2799', ['Mg-II', 2799.117]),  # SDSS
                ('Ne-V 3346', ['Ne-V', 3346.79]),
                ('Ne-VI 3426', ['Ne-VI', 3426.85]),
                ('[OII] 3727', ['[OII]', 3727.092]),  # [OII] doublet
                ('[OII] 3730', ['[OII]', 3729.875]),
                ('[NeIII] 3869', ['[NeIII]', 3869.85]),
                ('H-8 3890', ['H-8', 3890.15]),
                ('H-eps 3971', ['H' + eps, 3971.20]),
                ('H-delta 4102', ['H' + delta, 4102.89]),
                ('H-gamma 4341', ['H' + gamma, 4341.68]),
                ('H-beta 4862', ['H' + beta, 4862.683]),
                ('[OIII] 4960', ['[OIII]', 4960.295]),
                ('[OIII] 5008', ['[OIII]', 5008.240]),
                ('[NI] 5199', ['[NI]', 5199]),
                ('HeI 5877', ['HeI', 5877.29]),
                ('[OI] 6302', ['[OI]', 6302.05]),
                ('[NII] 6549', ['[NII]', 6549.85]),
                ('H-alpha 6564', ['H' + alpha, 6564.61]),
                ('[NII] 6585', ['[NII]', 6585.28]),
                ('SII 6718', ['SII', 6718.29]),
                ('SII 6732', ['SII', 6732.67]),
                ('A:Ca(H) 3934', ['A:Ca(H)', 3934.78]),
                ('A:Ca(K) 3969', ['A:Ca(K)', 3969.59]),
                ('A:G-band 4300', ['A:G-band', 4300.4]),
                ('A:Mg-1 5167', ['A:Mg-1', 5167.3222]),
                ('A:Mg-2 5172', ['A:Mg-2', 5172.6847]),
                ('A:Mg-3 5183', ['A:Mg-3', 5183.6046]),
                ('A:Na 5894', ['A:Na', 5894.57]),
                ('A:H-delta 4102', ['A:H' + delta, 4102.89]),
                ('A:H-gamma 4341', ['A:H' + gamma, 4341.68]),
                ('A:H-beta 4862', ['A:H' + beta, 4862.68]),
                ('A:H-alpha 6564', ['A:H' + alpha, 6564.61])
                ])
        


def contResiduals(p, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    v = p.valuesdict()
    intcp = v['intercept']
    try:
        slope = v['slope']
        model = intcp + slope * x
    except BaseException:
        model = intcp
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data) / eps


def linesResiduals(p, x, data=None, eps=None):
    # unpack parameters:
    # extract .value attribute for each parameter
    v = p.valuesdict()
    model = 0
    n = len(v) // 3
    for i in range(n):
        li = 'l' + str(i) + '_'
        lc = v[li + 'center']
        la = v[li + 'amplitude']
        ls = v[li + 'sigma']
        model += la / (np.sqrt(2 * np.pi) * ls) * np.exp(-(x - lc) * (x - lc) * 0.5 / (ls * ls))
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data) / eps


def fitContinuum(sp):
    """Fit the continuum defined in the guess."""
    from lmfit import Parameters, minimize
    slope = sp.guess.slope
    intcpt = sp.guess.intcpt
    xg, yg = zip(*sp.guess.xy)
    xg = np.array(xg)
    xg *= 1. + sp.gal.z
    wc = sp.gal.w
    fc = sp.gal.f
    ec = sp.gal.e
    c = sp.gal.c
    idx = np.where(((wc > xg[0]) & (wc < xg[1]) & c == 1) | ((wc > xg[2]) & (wc < xg[3]) & c == 1))
    x = wc[idx]
    y = fc[idx]
    e = ec[idx]
    # Definition of the model
    fit_params = Parameters()
    fit_params.add('intercept', value=intcpt)
    if slope != 0:
        fit_params.add('slope', value=slope)
    out = minimize(contResiduals, fit_params, args=(x,), kws={'data':y,'eps':e}, method='leastsq')
    # out = minimize(contResiduals, fit_params, args=(x,), kws={'data': y, 'eps': e}, method='Nelder')
    par = out.params#.valuesdict()
    p1 = par['intercept']
    ic = p1.value
    eic = p1.stderr
    try:
        p2 = par['slope']
        s = p2.value
        es = p2.stderr
    except:
        s = 0.0
        es = 0.0
    return ic, eic, s, es


def fitLines(sp, intercept, slope):
    """Fit the lines defined in the guess."""
    from lmfit import Parameters, minimize
    # Select the input values for the fit
    z = sp.gal.z
    wc = sp.gal.w
    fc = sp.gal.f
    ec = sp.gal.e
    c = sp.gal.c
    xg, yg = zip(*sp.guess.xy)
    xg = np.array(xg)
    xg *= (1. + z)  # Back to observed
    idx = np.where((wc > xg[0]) & (wc < xg[3]) & c == 1)
    x = wc[idx]
    y = fc[idx]
    e = ec[idx]
    continuum = intercept + slope * x
    y -= continuum
    # Normalization
    norm = np.abs(np.median(continuum))
    print('Normalization factor ', norm)
    y /= norm
    e /= norm
    # Define the model
    fit_params = Parameters()
    # Define lines
    for i, line in enumerate(sp.emlines + sp.ablines):
        li = 'l' + str(i) + '_'
        x0 = line.x0 * (1. + z)
        fit_params.add(li + 'center', value=x0, min=(x0 - 10), max=(x0 + 10))
        A = line.A / norm
        print('amplitude ', A)
        if A > 0:
            fit_params.add(li + 'amplitude', value=A, min=0.1 * A, max=A * 10)
        else:
            fit_params.add(li + 'amplitude', value=A, max=0.1 * A, min=A * 10)
        sigma = line.fwhm / 2.355 * (1. + z)
        fit_params.add(li + 'sigma', value=sigma, min=sigma / 2., max=sigma * 2)
    # Minimize
    out = minimize(linesResiduals, fit_params, args=(x,),
                   kws={'data': y, 'eps': e}, method='leastsq')
    #               kws={'data': y, 'eps': e}, method='Nelder')
    # Return lines fitted parameters
    pars = out.params#.valuesdict()
    nlines = len(pars) // 3
    print("Number of lines fitted: ", nlines)
    linepars = []
    for i in range(nlines):
        li = 'l' + str(i) + '_'
        center = pars[li + 'center'].value  # Observed
        centerErr = pars[li + 'center'].stderr  # Observed
        sigma = pars[li + 'sigma'] .value   # Observed
        sigmaErr = pars[li + 'sigma'] .stderr   # Observed
        A =  pars[li+'amplitude'].value
        Aerr = pars[li+'amplitude'].stderr
        amplitude = A * norm / (np.sqrt(2 * np.pi) * sigma)
        amplitudeErr = amplitude * (Aerr/A + sigmaErr/sigma)
        linepars.append([center, centerErr, amplitude, amplitudeErr, sigma, sigmaErr])
    return linepars
