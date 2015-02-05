###############################################################################
# apogee.spec.stack: stack APOGEE spectra in various ways
###############################################################################
import numpy
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
        mask= True-numpy.isnan(spec)
    else:
        mask= mask.astype('bool')
        mask*= True-numpy.isnan(spec)
    out= numpy.zeros(spec.shape[1])+numpy.nan
    for ii in range(spec.shape[1]):
        out[ii]= numpy.median(spec[mask[:,ii],ii])
    return out
