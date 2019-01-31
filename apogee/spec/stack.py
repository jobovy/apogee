###############################################################################
# apogee.spec.stack: stack APOGEE spectra in various ways
###############################################################################
import numpy
_BIGERR= 10.**7.
def median(spec,mask=None):
    """
    NAME:
       median
    PURPOSE:
       median stack a set of spectra
    INPUT:
       spec - array of spectra (nspec,nwave)
       mask= (None) if set, use this mask (1/True for inclusion)
    OUTPUT:
       median spectrum
    HISTORY:
       2015-01-26 - Written - Bovy (IAS@KITP)
    """
    if mask is None:
        mask= True^numpy.isnan(spec)
    else:
        mask= mask.astype('bool')
        mask*= True^numpy.isnan(spec)
    out= numpy.zeros(spec.shape[1])+numpy.nan
    for ii in range(spec.shape[1]):
        out[ii]= numpy.median(spec[mask[:,ii],ii])
    return out

def invvar(spec,specerr=None,return_error=False):
    """
    NAME:
       invvar
    PURPOSE:
       Inverse-variance stack a set of spectra
    INPUT:
       spec - array of spectra (nspec,nwave)
       specerr - (None) if set, use these errors
       return_error= (False) if True, also return the error on the stack
    OUTPUT:
       (stacked spectrum,error)
    HISTORY:
       2015-02-05 - Written - Bovy (IAS@KITP)
    """
    if specerr is None:
        specerr= numpy.ones_like(spec)
    spec[numpy.isnan(spec)]= 2.
    specerr[numpy.isnan(spec)]= _BIGERR
    err2= 1./numpy.sum(1./specerr**2.,axis=0)
    out= numpy.sum(spec/specerr**2.,axis=0)*err2
    out[out == 2.]= numpy.nan
    if return_error:
        return (out,numpy.sqrt(err2))
    else:
        return out
