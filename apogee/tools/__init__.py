import os.path
import numpy
from scipy import optimize, interpolate
from . import path as appath
from . import download as download
import fitsio
import warnings
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
    _PARAM_SYMBOL= [index.strip().lower().decode("utf-8")  for index in indexArrays['PARAM_SYMBOL'].flatten()]
    _ELEM_SYMBOL= [index.strip().lower().decode("utf-8")  for index in indexArrays['ELEM_SYMBOL'].flatten()]
    _ELEM_NUMBER_DICT= dict((elem,
                             elements.__dict__[elem.capitalize()].number)
                            for elem in _ELEM_SYMBOL 
                            if elem != 'ci' and elem != 'tiii')
    _ELEM_NUMBER_DICT['CI']= elements.__dict__['C'].number
    _ELEM_NUMBER_DICT['TiII']= elements.__dict__['Ti'].number

# Detector limits used in pix2wv and wv2pix
apStarBlu_lo = 322
apStarBlu_hi = 3242
apStarGre_lo = 3648
apStarGre_hi = 6048
apStarRed_lo = 6412
apStarRed_hi = 8306
aspcapBlu_start = 0
aspcapGre_start = apStarBlu_hi-apStarBlu_lo+aspcapBlu_start
aspcapRed_start = apStarGre_hi-apStarGre_lo+aspcapGre_start
aspcapTotal = apStarRed_hi-apStarRed_lo+aspcapRed_start

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
        return _ELEM_SYMBOL.index(elem.lower().decode("utf-8"))
    except ValueError:
        raise KeyError("Element %s is not part of the APOGEE elements (can't do everything!) or something went wrong)" % elem)

def atomic_number(elem):
    """
    NAME:
       atomic_number
    PURPOSE:
       return the atomic number of a given element
    INPUT:
       elem - element
    OUTPUT:
       atomic number
    HISTORY:
       2015-03-10 - Written - Bovy (IAS)
    """
    try:
        return _ELEM_NUMBER_DICT[elem.lower()]
    except (NameError,KeyError):
        return elements.__dict__[elem.lower().capitalize()].number

def vac2air(wave,sdssweb=False):
    """
    NAME:
       vac2air
    PURPOSE:
       Convert from vacuum to air wavelengths (See Allende Prieto technical note: http://hebe.as.utexas.edu/apogee/docs/air_vacuum.pdf)
    INPUT:
       wave - vacuum wavelength in \AA
       sdssweb= (False) if True, use the expression from the SDSS website (http://classic.sdss.org/dr7/products/spectra/vacwavelength.html)
    OUTPUT:
       air wavelength in \AA
    HISTORY:
       2014-12-04 - Written - Bovy (IAS)
       2015-04-27 - Updated to CAP note expression - Bovy (IAS)
    """
    if sdssweb:
        return wave/(1.+2.735182*10.**-4.+131.4182/wave**2.+2.76249*10.**8./wave**4.)
    else:
        return wave/(1.+0.05792105/(238.0185-(10000./wave)**2.)+0.00167917/(57.362-(10000./wave)**2.))

def air2vac(wave,sdssweb=False):
    """
    NAME:
       air2vac
    PURPOSE:
       Convert from air to vacuum wavelengths (See Allende Prieto technical note: http://hebe.as.utexas.edu/apogee/docs/air_vacuum.pdf)
    INPUT:
       wave - air wavelength in \AA
       sdssweb= (False) if True, use the expression from the SDSS website (http://classic.sdss.org/dr7/products/spectra/vacwavelength.html)
    OUTPUT:
       vacuum wavelength in \AA
    HISTORY:
       2014-12-04 - Written - Bovy (IAS)
       2015-04-27 - Updated to CAP note expression - Bovy (IAS)
    """
    return optimize.brentq(lambda x: vac2air(x,sdssweb=sdssweb)-wave,
                           wave-20,wave+20.)

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
        out= numpy.zeros((spec.shape[0],7214),dtype=spec.dtype)
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,7214),dtype=spec.dtype)
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
        out= numpy.zeros((spec.shape[0],8575),dtype=spec.dtype)
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,8575),dtype=spec.dtype)
        spec= numpy.reshape(spec,(1,len(spec)))
    out[:,322:3242]= spec[:,:2920]
    out[:,3648:6048]= spec[:,2920:5320]
    out[:,6412:8306]= spec[:,5320:]
    if oneSpec:
        return out[0]
    else:
        return out

totalpix = 8575 # total pixels on the apstar grid
wv0 = 10**4.179 # starting wavelength
wvs = numpy.zeros(totalpix) # create empty array for wavelengths
# Assign all wavelengths
wvs[0] = wv0
for i in range(1,totalpix):
    wvs[i] = 10**(6e-6 + numpy.log10(wvs[i-1]))
pixels = numpy.arange(0,totalpix)
aspcapwvs = toAspcapGrid(wvs)
apStar_pixel_interp = interpolate.interp1d(wvs,pixels,kind='linear',
                                           bounds_error=False)

def pix2wv(pix,apStarWavegrid=False):
    """
    NAME:
       pix2wv
    PURPOSE:
       convert pixel to wavelength
    INPUT:
       pix - pixel (int), range of pixels (tuple) or list of pixels (list/numpy array)
             float input will be converted to integers
       apStarWavegrid = (False) uses aspcapStarWavegrid by default
    OUTPUT:
       wavelength(s) in Angstroms corresponding to input pixel(s)
    HISTORY:
       2016-10-18 - Written - Price-Jones
    """
    # choose wavelength array to source from
    if apStarWavegrid:
        wvlist = wvs
        maxpix = totalpix
    elif not apStarWavegrid:
        wvlist = aspcapwvs
        maxpix = aspcapTotal
    # Check input cases
    if isinstance(pix,float):
        pix = int(pix)
    if isinstance(pix,int):
        if pix >= 0 and pix < maxpix:
            return wvlist[pix]
        else:
            warnings.warn("pixel outside allowed pixel range",RuntimeWarning)
            return numpy.nan
    elif isinstance(pix,tuple):
        if pix[0] >= 0 and pix[1] < maxpix:
            return wvlist[int(pix[0]):int(pix[1]):int(pix[2])]
        else:
            warnings.warn("pixel bounds outside allowed pixel range",RuntimeWarning)
            return numpy.nan
    elif isinstance(pix,(list,numpy.ndarray)):
        wavelengths = numpy.zeros(len(pix))
        for p in range(len(pix)):
            if pix[p] >= 0 and pix[p] < maxpix:
                wavelengths[p] = wvlist[int(pix[p])]
            else:
                warnings.warn("pixel outside allowed pixel range",RuntimeWarning)
                wavelength[p] = numpy.nan
        return wavelengths
    # If input not recognized inform the user
    elif not isinstance(wv,(int,float,tuple,list,numpy.ndarray)):
        warnings.warn("unrecognized pixel input",RuntimeWarning)
        return None

def wv2pix(wv,apStarWavegrid=False):
    """
    NAME:
       wv2pix
    PURPOSE:
       convert wavelength to pixel using interpolated function
    INPUT:
       wv - wavelength (int), range of wavelengths (tuple) or list of wavelengths
            (list/numpy array) in Angstroms
       apStarWavegrid = (False) uses aspcapStarWavegrid by default
    OUTPUT:
       array of pixel(s) corresponding to input wavelength(s)
       nan - indicates input wavelength(s) outside the range
       None - indicates input wavelength type not recognized
       0 - indicates the wavelength can be found but is outside the bounds of the a
           spcapStarWavegrid
    HISTORY:
       2016-10-18 - Written - Price-Jones
    """
    # Check input cases
    if isinstance(wv,(int,float)):
        if wv >= wvs[0] and wv <= wvs[-1]:
            pixels = apStar_pixel_interp(wv)
        else:
            warnings.warn("wavelength outside allowed wavelength range",RuntimeWarning)
            return numpy.nan 
    elif isinstance(wv,tuple):
        if wv[0] >= wvs[0] and wv[1] <= wvs[-1]:
            wvlist = numpy.arange(wv[0],wv[1],wv[2])
            pixels = apStar_pixel_interp(wvlist)
        else:
            warnings.warn("wavelength bounds outside allowed wavelength range",RuntimeWarning)
            return numpy.nan       
    elif isinstance(wv,(list,numpy.ndarray)):
        pixels = numpy.zeros(len(wv))
        allinside=True
        for w in range(len(wv)):
            if wv[w] >= wvs[0] and wv[w] <= wvs[-1]:
                pixels[w] = apStar_pixel_interp(wv[w])
            else:
                allinside=False
                pixels[w] = numpy.nan
        if not allinside:
            warnings.warn("wavelength outside allowed wavelength range",RuntimeWarning)
    # If input not recognized inform the user
    elif not isinstance(wv,(int,float,tuple,list,numpy.ndarray)):
        warnings.warn("unrecognized wavelength input",RuntimeWarning)
        return None

    if apStarWavegrid:
        return pixels.astype(int)

    # If on aspcapStarWavegrid, convert appropriately    
    elif not apStarWavegrid:        
        # find where pixel list matches detectors
        blue = numpy.where((pixels >= apStarBlu_lo) & (pixels < apStarBlu_hi))
        green = numpy.where((pixels >= apStarGre_lo) & (pixels < apStarGre_hi))
        red = numpy.where((pixels >= apStarRed_lo) & (pixels < apStarRed_hi))
        # find where pixel list does not match detectors
        if pixels.size > 1:
            nomatch = (numpy.array([i for i in range(len(pixels)) if i not in blue[0] and i not in green[0] and i not in red[0]]),)
             # adjust pixel values to match aspcap wavegrid
            pixels[blue] -= (apStarBlu_lo-aspcapBlu_start)
            pixels[green] -= (apStarGre_lo-aspcapGre_start)
            pixels[red] -= (apStarRed_lo-aspcapRed_start)
        # Case of single wavelength
        elif pixels.size == 1:
            if blue[0].size==1:
                pixels -= (apStarBlu_lo-aspcapBlu_start)
            elif green[0].size==1:
                pixels -= (apStarGre_lo-aspcapGre_start)
            elif red[0].size==1:
                pixels -= (apStarRed_lo-aspcapRed_start)
            elif blue[0].size==0 and green[0].size==0 and red[0].size==0:
                nomatch = ([0],)
                pixels = 0
        return numpy.floor(pixels).astype(int)
