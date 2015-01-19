###############################################################################
# apogee.spec.plot: various way to plot APOGEE spectra
###############################################################################
from functools import wraps
import numpy
from galpy.util import bovy_plot
import apogee.tools.read as apread
_LOG10LAMBDA0= 4.179 
_DLOG10LAMBDA= 6.*10.**-6.
_NLAMBDA= 8575
def specPlotInputDecorator(func):
    """Decorator to parse input to spectral plotting"""
    @wraps(func)
    def input_wrapper(*args,**kwargs):
        if len(args) == 2 and isinstance(args[0],(list,numpy.ndarray)) \
                and isinstance(args[1],(list,numpy.ndarray)):
            # wavelength, spectrum
            return func(args[0],args[1],**kwargs)
        elif len(args) == 1 and isinstance(args[0],(list,numpy.ndarray)):
            # spectrum on standard re-sampled wavelength grid
            lam= 10.**numpy.arange(_LOG10LAMBDA0,
                                   _LOG10LAMBDA0+_NLAMBDA*_DLOG10LAMBDA,
                                   _DLOG10LAMBDA)
            return func(lam,args[0],**kwargs)
        elif isinstance(args[0],(int,str)) and isinstance(args[1],str):
            # location ID and APOGEE ID (loc ID can be string for 1m sample)
            if kwargs.get('apStar',False):
                spec, hdr= apread.apStar(args[0],args[1],header=True,
                                         ext=kwargs.get('ext',1))
                spec= spec[numpy.amin([kwargs.get('apStarIndx',1),
                                       len(spec)-1])]
            else: #aspcapStar
                spec, hdr= apread.aspcapStar(args[0],args[1],header=True,
                                             ext=kwargs.get('ext',1))
            lam= 10.**numpy.arange(hdr['CRVAL1'],
                                   hdr['CRVAL1']+len(spec)*hdr['CDELT1'],
                                   hdr['CDELT1'])
            return func(lam,spec,**kwargs)
    return input_wrapper

@specPlotInputDecorator
def chunks(*args,**kwargs):
    """
    NAME:
       chunks
    PURPOSE:
       plot chunks of the spectrum in one row
    INPUT:
       Either:
          (a) wavelength, spectrum (\AA,spectrum units)
          (b) spectrum (assumed on standard APOGEE re-sampled wavelength grid)
          (c) location ID, APOGEE ID (default loads aspcapStar, loads extension ext(=1); apStar=True loads apStar spectrum)
    KEYWORDS:
       ext= (1) extension to load
       apStar= (False) if True, load the apStar spectrum
       apStarIndx= (1) index in the apStar spectrum to load
       bovy_plot.bovy_plot kwargs
    OUTPUT:
       plot to output
    HISTORY:
       2015-01-18 - Written (based on older code) - Bovy (IAS)
    """
    bovy_plot.bovy_plot(args[0],args[1])
    return None
