from functools import wraps
import numpy
import apogee.tools.read as apread
from apogee.tools import toAspcapGrid
def specFitInput(func):
    """Decorator to parse the input for spectral fitting"""
    @wraps(func)
    def input_wrapper(*args,**kwargs):
        spec= args[0]
        specerr= args[1]
        if isinstance(specerr,str): # locID+APOGEE-ID; array
            ispec= apread.aspcapStar(spec,specerr,ext=1,header=False,
                                     aspcapWavegrid=True)
            ispecerr= apread.aspcapStar(spec,specerr,ext=2,header=False,
                                        aspcapWavegrid=True)
            spec= ispec
            specerr= ispecerr
        elif (isinstance(specerr,(list,numpy.ndarray)) \
                  and isinstance(specerr[0],str)): # locID+APOGEE-ID; array
            nspec= len(specerr)
            ispec= numpy.empty((nspec,7214))
            ispecerr= numpy.empty((nspec,7214))
            for ii in range(nspec):
                ispec[ii]= apread.aspcapStar(spec[ii],specerr[ii],ext=1,
                                             header=False,aspcapWavegrid=True)
                ispecerr[ii]= apread.aspcapStar(spec[ii],specerr[ii],ext=2,
                                                header=False,aspcapWavegrid=True)
            spec= ispec
            specerr= ispecerr
        elif isinstance(specerr,(list,numpy.ndarray)) \
                and isinstance(specerr[0],(float,numpy.float32,
                                           numpy.float64,numpy.ndarray)) \
            and ((len(specerr.shape) == 1 and len(specerr) == 8575)
                 or (len(specerr.shape) == 2 and specerr.shape[1] == 8575)): #array on apStar grid
            spec= toAspcapGrid(spec)
            specerr= toAspcapGrid(specerr)
        return func(spec,specerr,*args[2:],**kwargs)
    return input_wrapper

