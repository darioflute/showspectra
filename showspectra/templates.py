import numpy as np
from astropy.io.fits import getdata
import os, collections

class Template(object):
    """ SDSS template spectrum """
    def __init__(self,wave,flux,template):
        self.w = wave
        self.f = flux     # flux for cross-correlation
        self.t = template # flux for display
            
def readSdssTemplate(file):
    data, hdr = getdata(file, 0, header=True)
    cdelt=hdr['cd1_1']
    crval=hdr['crval1']
    wave=crval+np.arange(len(data[0]))*cdelt
    wave = np.exp(np.log(10.)*wave)
    return Template(wave,data[1],data[0])
        
def readTemplates(self):
    """ Read templates for cross correlation"""
    # spDR2-023.fit -- Early-type galaxy
    # spDR2-027.fit -- Late-type galaxy
    # spDR2-029.fit -- QSO
    # spDR2-032.fit -- High-luminosity QSO
    # extension 1 contains the spectrum after median filtering subtraction
    here = os.path.abspath(os.path.dirname(__file__))
    path = here + '/Templates/'
    self.EarlyTypeGalaxy = readSdssTemplate(path+'spDR2-023.fit')
    self.Galaxy1 = readSdssTemplate(path+'spDR2-024.fit')
    self.Galaxy2 = readSdssTemplate(path+'spDR2-025.fit')
    self.Galaxy3 = readSdssTemplate(path+'spDR2-026.fit')
    self.LateTypeGalaxy = readSdssTemplate(path+'spDR2-027.fit')
    self.LumRedGalaxy = readSdssTemplate(path+'spDR2-028.fit')
    self.QSO = readSdssTemplate(path+'spDR2-029.fit')
    #        self.BALQSO1 = readSdssTemplate(path+'spDR2-030.fit')
    #        self.BALQSO2 = readSdssTemplate(path+'spDR2-031.fit')
    #        self.HighLumQSO = readSdssTemplate(path+'spDR2-032.fit')
    
    self.templates  =  collections.OrderedDict()
    self.templates['EarlyTypeGalaxy'] = self.EarlyTypeGalaxy
    self.templates['LateTypeGalaxy'] = self.LateTypeGalaxy
    self.templates['LumRedGalaxy'] = self.LumRedGalaxy
    self.templates['Galaxy1'] = self.Galaxy1
    self.templates['Galaxy2'] = self.Galaxy2
    self.templates['Galaxy3'] = self.Galaxy3
    #        self.templates['QSO'] = self.QSO
    #        self.templates['BALQSO1'] = self.BALQSO1
    #        self.templates['BALQSO2'] = self.BALQSO2
    #        self.templates['HighLumQSO'] = self.HighLumQSO
        
        
