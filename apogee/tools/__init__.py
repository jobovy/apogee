import os.path
import numpy
from scipy import optimize
import path as appath
import download as download
import fitsio
from periodictable import elements
try:
    # Need to have allStar
    filePath= appath.allStarPath()
    if not os.path.exists(filePath):
        download.allStar()
    indexArrays= fitsio.read(appath.allStarPath(),3)
except ValueError:
    _INDEX_ARRAYS_LOADED= False
else:
    _INDEX_ARRAYS_LOADED= True
    _PARAM_SYMBOL= [index.strip().lower() for index in indexArrays['PARAM_SYMBOL'].flatten()]
    _ELEM_SYMBOL= [index.strip().lower() for index in indexArrays['ELEM_SYMBOL'].flatten()]
    _ELEM_NUMBER_DICT= dict((elem,
                             elements.__dict__[elem.capitalize()].number)
                            for elem in _ELEM_SYMBOL)
def paramIndx(param):
    """
    NAME:
       paramIndx
    PURPOSE:
       return the index into the PARAM/FPARAM  arrays corresponding to a given stellar parameter 
    INPUT:
       param - the stellar parameter (one of TEFF,LOGG,LOG10VDOP,METALS,C,N,ALPHA)
    OUTPUT:
       index into PARAM/FPARAM array
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    if not _INDEX_ARRAYS_LOADED: raise ImportError("paramIndx function cannot be used, because the allStar file could not be properly loaded")
    if param.lower() == 'alpha': return _PARAM_SYMBOL.index('o mg si s ca ti')
    else: 
        try:
            return _PARAM_SYMBOL.index(param.lower())
        except ValueError:
            raise KeyError("Stellar parameter %s not recognized" % param)

def elemIndx(elem):
    """
    NAME:
       elemIndx
    PURPOSE:
       return the index into the ELEM/FELEM arrays corresponding to a given element
    INPUT:
       elem - the element (string like 'C')
    OUTPUT:
       index into ELEM/FELEM array
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    if not _INDEX_ARRAYS_LOADED: raise ImportError("elemIndx function cannot be used, because the allStar file could not be properly loaded")
    try:
        return _ELEM_SYMBOL.index(elem.lower())
    except ValueError:
        raise KeyError("Element %s is not part of the APOGEE elements (can't do everything!) or something went wrong)" % elem)

def vac2air(wave):
    """
    NAME:
       vac2air
    PURPOSE:
       Convert from vacuum to air wavelengths (Morton 1991; see http://classic.sdss.org/dr7/products/spectra/vacwavelength.html)
    INPUT:
       wave - vacuum wavelength in \AA
    OUTPUT:
       air wavelength in \AA
    HISTORY:
       2014-12-04 - Written - Bovy (IAS)
    """
    return wave/(1.+2.735182*10.**-4.+131.4182/wave**2.+2.76249*10.**8./wave**4.)

def air2vac(wave):
    """
    NAME:
       air2vac
    PURPOSE:
       Convert from air to vacuum wavelengths (Morton 1991; see http://classic.sdss.org/dr7/products/spectra/vacwavelength.html)
    INPUT:
       wave - air wavelength in \AA
    OUTPUT:
       vacuum wavelength in \AA
    HISTORY:
       2014-12-04 - Written - Bovy (IAS)
    """
    return optimize.brentq(lambda x: vac2air(x)-wave,wave-100,wave+100.)

def toAspcapGrid(spec):
    """
    NAME:
       toAspcapGrid
    PURPOSE:
       convert a spectrum from apStar grid to the ASPCAP grid (w/o the detector gaps)
    INPUT:
       spec - spectrum (or whatever) on the apStar grid; either (nwave) or (nspec,nwave)
    OUTPUT:
       spectrum (or whatever) on the ASPCAP grid
    HISTORY:
       2015-02-17 - Written - Bovy (IAS)
    """
    if len(spec.shape) == 2: # (nspec,nwave)
        out= numpy.zeros((spec.shape[0],7214))
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,7214))
        spec= numpy.reshape(spec,(1,len(spec)))
    out[:,:2920]= spec[:,322:3242]
    out[:,2920:5320]= spec[:,3648:6048]
    out[:,5320:]= spec[:,6412:8306]
    if oneSpec:
        return out[0]
    else:
        return out

def toApStarGrid(spec):
    """
    NAME:
       toApStarGrid
    PURPOSE:
       convert a spectrum from the ASPCAP grid (w/o the detector gaps) to the apStar grid
    INPUT:
       spec - spectrum (or whatever) on the ASPCAP grid; either (nwave) or (nspec,nwave)
    OUTPUT:
       spectrum (or whatever) on the apStar grid
    HISTORY:
       2015-02-17 - Written - Bovy (IAS)
    """
    if len(spec.shape) == 2: # (nspec,nwave)
        out= numpy.zeros((spec.shape[0],8575))
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,8575))
        spec= numpy.reshape(spec,(1,len(spec)))
    out[:,322:3242]= spec[:,:2920]
    out[:,3648:6048]= spec[:,2920:5320]
    out[:,6412:8306]= spec[:,5320:]
    if oneSpec:
        return out[0]
    else:
        return out
