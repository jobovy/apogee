import os.path
import numpy
from scipy import optimize, interpolate
from . import path as appath
from . import download as download
try:
    import fitsio
    fitsread = fitsio.read
except ImportError:
    import astropy.io.fits as pyfits
    fitsread= pyfits.getdata
import warnings
from periodictable import elements
try:
    # Need to have allStar
    filePath= appath.allStarPath()
    if not os.path.exists(filePath):
        download.allStar()
    indexArrays= fitsread(appath.allStarPath(),3)
except ValueError:
    _INDEX_ARRAYS_LOADED= False
else:
    _INDEX_ARRAYS_LOADED= True
    if type(indexArrays['PARAM_SYMBOL'][0,0]) == numpy.dtype(bytes):
        _PARAM_SYMBOL= [index.strip().lower().decode("utf-8")  
                        for index in indexArrays['PARAM_SYMBOL'].flatten()]
        _ELEM_SYMBOL= [index.strip().lower().decode("utf-8") 
                       for index in indexArrays['ELEM_SYMBOL'].flatten()]
    else:
        _PARAM_SYMBOL= [index.strip().lower()  
                        for index in indexArrays['PARAM_SYMBOL'].flatten()]
        _ELEM_SYMBOL= [index.strip().lower()
                       for index in indexArrays['ELEM_SYMBOL'].flatten()]
    _ELEM_NUMBER_DICT= dict((elem,
                             elements.__dict__[elem.capitalize()].number)
                            for elem in _ELEM_SYMBOL 
                            if elem != 'ci' and elem != 'tiii')
    _ELEM_NUMBER_DICT['CI']= elements.__dict__['C'].number
    _ELEM_NUMBER_DICT['TiII']= elements.__dict__['Ti'].number


# DR12 abundance uncertainty coefficients  as a function of Teff, [M/H], SNR
# from http://www.sdss.org/dr12/irspec/abundances/
# see also Holtzman et al 2015

_ch_12coeff=[-3.350,0.769,-0.919,-0.066]
_nh_12coeff=[-2.704,0.291,-0.591,-0.078]
_oh_12coeff=[-3.649,0.670,-0.614,-0.093]
_nah_12coeff=[-2.352,-0.002,-0.915,-0.263]
_mgh_12coeff=[-3.537,0.263,-0.825,-0.297]
_alh_12coeff=[-2.764,0.471,-0.868,-0.162]
_sih_12coeff=[-3.150,0.383,-0.224,-0.105]
_sh_12coeff=[-3.037,0.507,-0.625,-0.299]
_kh_12coeff=[-2.770,0.216,-0.667,-0.275]
_cah_12coeff=[-3.226,0.284,-0.879,-0.429]
_tih_12coeff=[-3.186,0.657,-0.819,-0.068]
_vh_12coeff=[-1.608,0.900,-0.400,-0.418]
_mnh_12coeff=[-3.031,0.639,-0.661,-0.326]
_feh_12coeff=[-3.357,0.098,-0.303,-0.071]
_nih_12coeff=[-3.153,0.135,-0.493,-0.185]
_mh_12coeff=[-3.603,0.109,-0.433,0.039]
_alpha_12coeff=[-4.360,0.060,-0.848,-0.096]

DR12_XH_coeff = {'C_H':_ch_12coeff,'N_H':_nh_12coeff,'O_H':_oh_12coeff,
                 'NA_H':_nah_12coeff,'MG_H':_mgh_12coeff,'AL_H':_alh_12coeff,
                 'SI_H':_sih_12coeff,'S_H':_sh_12coeff,'K_H':_kh_12coeff,
                 'CA_H':_cah_12coeff,'TI_H':_tih_12coeff,'V_H':_vh_12coeff,
                 'MN_H':_mnh_12coeff,'FE_H':_feh_12coeff,'NI_H':_nih_12coeff,
                 'METALS':_mh_12coeff,'ALPHAFE':_alpha_12coeff}


# DR13 abundance uncertainty coefficients  as a function of Teff, [M/H], SNR
# from http://www.sdss.org/dr13/irspec/abundances/

_cfe_13coeff=[-3.243,0.608,-0.757,-0.257]
_cIfe_13coeff=[-2.804,0.403,-0.743,-0.319]
_nfe_13coeff=[-2.671,0.373,-0.407,-0.192]
_ofe_13coeff=[-3.410,1.471,-0.778,-0.182]
_nafe_13coeff=[-2.389,0.140,-0.926,-0.323]
_mgfe_13coeff=[-3.980,0.284,-0.949,-0.115]
_alfe_13coeff=[-2.616,-0.192,-0.628,-0.399]
_sife_13coeff=[-3.464,0.548,-0.482,-0.212]
_pfe_13coeff=[-1.988,0.384,-0.568,-0.369]
_sfe_13coeff=[-2.199,-0.030,-0.402,-0.295]
_kfe_13coeff=[-3.098,0.208,-0.583,-0.496]
_cafe_13coeff=[-3.520,0.153,-0.895,-0.405]
_tife_13coeff=[-3.108,0.295,-0.741,-0.185]
_tiIIfe_13coeff=[-2.192,0.328,-0.538,-0.267]
_vfe_13coeff=[-2.447,1.030,-1.096,-0.519]
_crfe_13coeff=[-3.191,0.290,-0.775,-0.455]
_mnfe_13coeff=[-3.523,0.235,-0.614,-0.488]
_feh_13coeff=[-5.316,0.202,-0.874,0.019]
_cofe_13coeff=[-2.062,1.064,-0.656,-0.523]
_nife_13coeff=[-4.067,0.442,-0.816,-0.395]
_cufe_13coeff=[-2.140,-0.096,-0.559,-0.426]
_gefe_13coeff=[-1.893,0.258,-0.665,-0.395]
_rbfe_13coeff=[-2.325,0.466,-1.117,-0.360]
_mh_13coeff=[-3.730,0.232,-0.524,0.013]
_alpha_13coeff=[-4.219,0.053,-0.794,-0.127]

DR13_XH_coeff={'C_FE':_cfe_13coeff,'CI_FE':_cIfe_13coeff,'N_FE':_nfe_13coeff,
               'O_FE':_ofe_13coeff,'NA_FE':_nafe_13coeff,'MG_FE':_mgfe_13coeff,
               'AL_FE':_alfe_13coeff,'SI_FE':_sife_13coeff,'P_FE':_pfe_13coeff,
               'S_FE':_sfe_13coeff,'K_FE':_kfe_13coeff,'CA_FE':_cafe_13coeff,
               'TI_FE':_tife_13coeff,'TIII_FE':_tiIIfe_13coeff,
               'V_FE':_vfe_13coeff,'CR_FE':_crfe_13coeff,'MN_FE':_mnfe_13coeff,
               'FE_H':_feh_13coeff,'CO_FE':_cofe_13coeff,'NI_FE':_nife_13coeff,
               'CU_FE':_cufe_13coeff,'GE_FE':_gefe_13coeff,
               'RB_FE':_rbfe_13coeff,'M_H':_mh_13coeff,
               'ALPHA_M':_alpha_13coeff}

drcoeffs = {'12':DR12_XH_coeff,'13':DR13_XH_coeff}


# Detector limit by data release
apStarInds = {'10':{'blue':(322,3242),'green':(3648,6048),'red':(6412,8306)},
              '11':{'blue':(322,3242),'green':(3648,6048),'red':(6412,8306)},
              '12':{'blue':(322,3242),'green':(3648,6048),'red':(6412,8306)},
              '13':{'blue':(246,3274),'green':(3585,6080),'red':(6344,8335)},
              '14':{'blue':(246,3274),'green':(3585,6080),'red':(6344,8335)},
              '16':{'blue':(246,3274),'green':(3585,6080),'red':(6344,8335)},
              'current':{'blue':(246,3274),'green':(3585,6080),'red':(6344,8335)}
             }

def _apStarPixelLimits(dr=None):
  """
  NAME: 
      _apStarPixelLimits
  PURPOSE:
      return the apStar pixel bounds for each detector for the chosen data 
      release by unpacking apStarInds.
  INPUT
      dr - string referring to data release, e.g. '12' 
  OUTPUT:
      bounds of blue, green and red detectors.
  HISTORY:
      2018-02-05 - Written - Price-Jones (UofT)
  """
  if dr is None: 
    dr=appath._default_dr()
  inds = apStarInds[dr]
  apStarBlu_lo,apStarBlu_hi = inds['blue']
  apStarGre_lo,apStarGre_hi = inds['green']
  apStarRed_lo,apStarRed_hi = inds['red']
  return apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi

def _aspcapPixelLimits(dr=None):
  """
  NAME: 
      _aspcapPixelLimits
  PURPOSE:
      return the ASPCAP pixel bounds for each detector for the chosen data 
      release by unpacking apStarInds.
  INPUT
      dr - string referring to data release, e.g. '12' 
  OUTPUT:
      starting pixel of the blue, green and red detectors, as well as the 
      total spectrum length in pixels
  HISTORY:
      2018-02-05 - Written - Price-Jones (UofT)
  """
  if dr is None:
    dr=appath._default_dr()
  apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=dr)
  aspcapBlu_start = 0
  aspcapGre_start = apStarBlu_hi-apStarBlu_lo+aspcapBlu_start
  aspcapRed_start = apStarGre_hi-apStarGre_lo+aspcapGre_start
  aspcapTotal = apStarRed_hi-apStarRed_lo+aspcapRed_start
  return aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal

# Wavegrid parameters used in apStarWavegrid and pix2wv
_LOG10LAMBDA0= 4.179 
_DLOG10LAMBDA= 6.*10.**-6.
_NLAMBDA= 8575

def apStarWavegrid():
    return 10.**numpy.arange(_LOG10LAMBDA0,
                             _LOG10LAMBDA0+_NLAMBDA*_DLOG10LAMBDA,
                             _DLOG10LAMBDA)

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

def sigma_XH(elem,Teff=4500.,M_H=0.,SNR=100.,dr=None):
    """
    NAME:
       sigma_XH
    PURPOSE:
       return uncertainty in a given element at specified effective 
       temperature, metallicity and signal to noise ratio (functional form
       taken from Holtzman et al 2015)
    INPUT:
       elem - string element name following the ASPCAP star naming convention
              i.e. for DR12 carbon, string is 'C_H'
       Teff - effective temperature or array thereof  in K, defaults to 4500 K
       M_H  - metallicity or array thereof, defaults to 0
       SNR  - signal to noise ratio or array thereof, defaults to 100
       dr   - data release
    OUTPUT:
       float or array depending on shape of Teff, M_H and SNR input
    HISTORY:
       2017-07-24 - Written - Price-Jones (UofT)
    """
    if dr is None: dr=appath._default_dr()
    A,B,C,D = drcoeffs[dr][elem]
    logsig = A + B*((Teff-4500.)/1000.) + C*M_H + D*(SNR-100)
    return numpy.exp(logsig)


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

def toAspcapGrid(spec,dr=None):
    """
    NAME:
       toAspcapGrid
    PURPOSE:
       convert a spectrum from apStar grid to the ASPCAP grid (w/o the detector gaps)
    INPUT:
       spec - spectrum (or whatever) on the apStar grid; either (nwave) or (nspec,nwave)
       dr - data release of pixel bounds to use
    OUTPUT:
       spectrum (or whatever) on the ASPCAP grid
    HISTORY:
       2015-02-17 - Written - Bovy (IAS)
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones (UofT)
    """
    apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=dr)    
    aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=dr)
    if len(spec.shape) == 2: # (nspec,nwave)
        out= numpy.zeros((spec.shape[0],aspcapTotal),dtype=spec.dtype)
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,aspcapTotal),dtype=spec.dtype)
        spec= numpy.reshape(spec,(1,len(spec)))
    out[:,:aspcapGre_start]= spec[:,apStarBlu_lo:apStarBlu_hi]
    out[:,aspcapGre_start:aspcapRed_start]= spec[:,apStarGre_lo:apStarGre_hi]
    out[:,aspcapRed_start:]= spec[:,apStarRed_lo:apStarRed_hi]
    if oneSpec:
        return out[0]
    else:
        return out

def toApStarGrid(spec,dr=None):
    """
    NAME:
       toApStarGrid
    PURPOSE:
       convert a spectrum from the ASPCAP grid (w/o the detector gaps) to the apStar grid
    INPUT:
       spec - spectrum (or whatever) on the ASPCAP grid; either (nwave) or (nspec,nwave)
       dr - data release of pixel bounds to use
    OUTPUT:
       spectrum (or whatever) on the apStar grid
    HISTORY:
       2015-02-17 - Written - Bovy (IAS)
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones (UofT)
    """
    apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=dr)    
    aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=dr)
    if len(spec.shape) == 2: # (nspec,nwave)
        out= numpy.zeros((spec.shape[0],8575),dtype=spec.dtype)
        oneSpec= False
    else:
        oneSpec= True
        out= numpy.zeros((1,8575),dtype=spec.dtype)
        spec= numpy.reshape(spec,(1,len(spec)))
    out[:,apStarBlu_lo:apStarBlu_hi]= spec[:,:aspcapGre_start]
    out[:,apStarGre_lo:apStarGre_hi]= spec[:,aspcapGre_start:aspcapRed_start]
    out[:,apStarRed_lo:apStarRed_hi]= spec[:,aspcapRed_start:]
    if oneSpec:
        return out[0]
    else:
        return out

wvs = apStarWavegrid()
aspcapwvs = toAspcapGrid(wvs)
pixels = numpy.arange(0,_NLAMBDA)
apStar_pixel_interp = interpolate.interp1d(wvs,pixels,kind='linear',
                                           bounds_error=False)

def pix2wv(pix,apStarWavegrid=False,dr=None):
    """
    NAME:
       pix2wv
    PURPOSE:
       convert pixel to wavelength
    INPUT:
       pix - pixel (int), range of pixels (tuple) or list of pixels (list/numpy array)
             float input will be converted to integers
       apStarWavegrid = (False) uses aspcapStarWavegrid by default
       dr - data release of pixel bounds to use
    OUTPUT:
       wavelength(s) in Angstroms corresponding to input pixel(s)
    HISTORY:
       2016-10-18 - Written - Price-Jones
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones (UofT)
    """
    # choose wavelength array to source from
    if apStarWavegrid:
        wvlist = wvs
        maxpix = _NLAMBDA
    elif not apStarWavegrid:
        aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=dr)
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
        if isinstance(pix,list):
          pix = numpy.array(pix)
        wavelengths = numpy.zeros(len(pix))
        valid = (pix>=0) & (pix<maxpix)
        invalid = (pix<0) | (pix>maxpix)
        wavelengths[valid] = wvlist[pix[valid].astype(int)]
        wavelengths[invalid] = numpy.nan
        if sum(invalid)!=0:
            warnings.warn("pixel outside allowed pixel range",RuntimeWarning)
        return wavelengths
    # If input not recognized inform the user
    elif not isinstance(wv,(int,float,tuple,list,numpy.ndarray)):
        warnings.warn("unrecognized pixel input",RuntimeWarning)
        return None

def wv2pix(wv,apStarWavegrid=False,dr=None):
    """
    NAME:
       wv2pix
    PURPOSE:
       convert wavelength to pixel using interpolated function
    INPUT:
       wv - wavelength (int), range of wavelengths (tuple) or list of wavelengths
            (list/numpy array) in Angstroms
       apStarWavegrid = (False) uses aspcapStarWavegrid by default
       dr - data release of pixel bounds to use
    OUTPUT:
       array of pixel(s) corresponding to input wavelength(s)
       nan - indicates input wavelength(s) outside the range
       None - indicates input wavelength type not recognized
       0 - indicates the wavelength can be found but is outside the bounds of the a
           spcapStarWavegrid
    HISTORY:
       2016-10-18 - Written - Price-Jones
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones
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
        if isinstance(wv,list):
          wv = numpy.array(wv)
        pixels = numpy.zeros(len(wv))
        valid = (wv>=wvs[0]) & (wv<wvs[-1])
        invalid = (wv<wvs[0]) | (wv>wvs[-1])
        pixels[valid] = apStar_pixel_interp(wv[valid])
        pixels[invalid] = numpy.nan
        if sum(invalid)!=0:
            warnings.warn("wavelength outside allowed wavelength range",RuntimeWarning)
    # If input not recognized inform the user
    elif not isinstance(wv,(int,float,tuple,list,numpy.ndarray)):
        warnings.warn("unrecognized wavelength input",RuntimeWarning)
        return None

    if apStarWavegrid:
        return pixels.astype(int)

    # If on aspcapStarWavegrid, convert appropriately    
    elif not apStarWavegrid:
        apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=dr)    
        aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=dr)
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
