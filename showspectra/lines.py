import numpy as np

def define_lines():
    import collections
    alpha = u'\u03B1'
    beta = u'\u03B2'
    delta = u'\u03B3'
    gamma = u'\u03B4'
    eps = u'\u03B5'
    #        alpha.encode('utf8')
    return collections.OrderedDict([
        ('OVI 1033', ['OVI', 1033.83]),
        ('Ly-alpha 1215', ['Ly' + alpha, 1215.67]),
        ('N-V 1240', ['N-V', 1240.14]),
        ('OI 1304', ['OI', 1304.35]),
        ('Si-IV 1396', ['Si-IV', 1396.76]),
        ('C-IV 1549', ['C-IV', 1549.06]),
        ('C-III 1908', ['C-III', 1908.73]),
        ('Mg-II 2798', ['Mg-II', 2798.75]),
        ('[OII] 3728', ['[OII]', 3728.48]),
        ('[NeIII] 3869', ['[NeIII]', 3869.85]),
        ('H-8 3890', ['H-8', 3890.15]),
        ('H-eps 3971', ['H' + eps, 3971.20]),
        ('H-delta 4102', ['H' + delta, 4102.89]),
        ('H-gamma 4341', ['H' + gamma, 4341.68]),
        ('H-beta 4862', ['H' + beta, 4862.68]),
        ('[OIII] 4960', ['[OIII]', 4960.30]),
        ('[OIII] 5008', ['[OIII]', 5008.24]),
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

def contResiduals(p,x,data=None,eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    v = p.valuesdict()
    intcp = v['intercept']
    slope = v['slope']    
    model = intcp+slope*x
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data)/eps

def linesResiduals(p,x,data=None,eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    v = p.valuesdict()
    model = 0
    n = len(v)//3
    for i in range(n):
        li = 'l'+str(i)+'_'
        lc = v[li+'center']
        la = v[li+'amplitude']
        ls = v[li+'sigma']
        model += la/(np.sqrt(2*np.pi)*ls)*np.exp(-(x-lc)*(x-lc)*0.5/(ls*ls))
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data)/eps

def fitContinuum(sp):
    """Fit the continuum defined in the guess."""
    from lmfit import Parameters, minimize
    slope = sp.guess.slope
    intcpt = sp.guess.intcpt
    xg,yg = zip(*sp.guess.xy)
    xg = np.array(xg)
    xg *= 1. + sp.gal.z
    wc = sp.gal.wc
    fc = sp.gal.fc
    ec = sp.gal.ec
    c = sp.gal.c
    idx = np.where(((wc > xg[0]) & (wc < xg[1]) & c==1) | ((wc > xg[2]) & (wc < xg[3]) & c==1))
    x = wc[idx]
    y = fc[idx]
    e = ec[idx]
    # Definition of the model
    fit_params = Parameters()
    fit_params.add('intercept', value=intcpt)
    fit_params.add('slope', value=slope)
    # out = minimize(contResiduals, fit_params, args=(x,), kws={'data':y,'eps':e},method='leastsq')
    out = minimize(contResiduals, fit_params, args=(x,), kws={'data':y,'eps':e},method='Nelder')
    par = out.params.valuesdict() 
    return par['intercept'], par['slope']

def fitLines(sp,intercept,slope):
    """Fit the lines defined in the guess."""
    from lmfit import Parameters, minimize
    # Select the input values for the fit
    z = sp.gal.z
    wc = sp.gal.wc
    fc = sp.gal.fc
    ec = sp.gal.ec
    c = sp.gal.c
    xg,yg = zip(*sp.guess.xy)
    xg = np.array(xg)
    xg *= (1. + sp.gal.z)  # Back to observed
    idx = np.where((wc > xg[0]) & (wc < xg[3]) & c==1)
    x = wc[idx]
    y = fc[idx]
    e = ec[idx]
    continuum = intercept + slope * x
    y -= continuum
    # Normalization
    norm = np.abs(np.median(continuum))
    print('Normalization factor ',norm)
    y /= norm
    e /= norm
    # Define the model
    fit_params = Parameters()
    # Define lines
    for i,line in enumerate(sp.emlines+sp.ablines):
        li = 'l'+str(i)+'_'
        x0 = line.x0 * (1. + z)
        fit_params.add(li+'center', value=x0,min=(x0-10), max=(x0+10))
        A = line.A/norm
        print('amplitude ', A)
        if A > 0:
            fit_params.add(li+'amplitude', value=A, min=0.1*A, max=A*10)
        else:
            fit_params.add(li+'amplitude', value=A, max=0.1*A, min=A*10)
        sigma = line.fwhm/2.355 * (1. + z)
        fit_params.add(li+'sigma', value=sigma, min=sigma/2., max=sigma*2)
    # Minimize
    # out = minimize(linesResiduals, fit_params, args=(x,), kws={'data':y,'eps':e},method='leastsq')
    out = minimize(linesResiduals, fit_params, args=(x,), kws={'data':y,'eps':e},method='Nelder')            
    # Return lines fitted parameters
    pars = out.params.valuesdict()
    nlines = len(pars)//3
    print ("Number of lines fitted: ", nlines)
    linepars = []
    for i in range(nlines):
        li = 'l'+str(i)+'_'
        center = pars[li+'center']
        sigma = pars[li+'sigma']
        amplitude = pars[li+'amplitude']*norm/(np.sqrt(2*np.pi)*sigma)
        linepars.append([center, amplitude, sigma])
    return linepars
    