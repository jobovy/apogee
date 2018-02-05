from functools import wraps
import os, os.path
import shutil
import subprocess
import numpy
from scipy import special
import apogee.tools.read as apread
import apogee.tools.path as appath
from apogee.tools import toAspcapGrid,_aspcapPixelLimits
from apogee.spec.plot import apStarWavegrid
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
            aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
            nspec= len(specerr)
            ispec= numpy.empty((nspec,aspcapTotal))
            ispecerr= numpy.empty((nspec,aspcapTotal))
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

def convert_modelAtmosphere(**kwargs):
    """
    NAME:
       convert_modelAtmosphere
    PURPOSE:
       Convert a model atmosphere to MOOG format
    INPUT:
       Either:
          (a) modelatm= (None) can be set to the filename of a model atmosphere
          (b) specify the stellar parameters for a grid point in model atm by
              - lib= ('kurucz_filled') spectral library
              - teff= (4500) grid-point Teff
              - logg= (2.5) grid-point logg
              - metals= (0.) grid-point metallicity
              - cfe= (0.) grid-point carbon-enhancement
              - afe= (0.) grid-point alpha-enhancement
              - dr= return the path corresponding to this data release
       vmicro= (2.) microturbulence (km/s) (only used if the MOOG-formatted atmosphere file doesn't already exist)
    OUTPUT:
       (none; just converts and caches the model atmosphere
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
       2015-03-21 - Adjusted to also work for off-grid atmosphers - Bovy (IAS)
    """
    # Get the filename of the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    if not modelatm is None:
        if isinstance(modelatm,str) and os.path.exists(modelatm):
            modelfilename= modelatm
        elif isinstance(modelatm,str):
            raise ValueError('modelatm= input is a non-existing filename')
        else: # model atmosphere instance
            raise ValueError('modelatm= in moogsynth should be set to the name of a file')
    else:
        modelfilename= appath.modelAtmospherePath(**kwargs)
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    outname= modelbasename.replace('.mod','.org')
    if os.path.exists(os.path.join(modeldirname,outname)): return None
    shutil.copy(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'scripts/makemoogmodel.awk'),modeldirname)
    try:
        stdout= open(os.path.join(modeldirname,outname),'w')
        stderr= open('/dev/null','w')
        subprocess.check_call(['awk','-f','makemoogmodel.awk',
                               'vmicro=%.1f' % kwargs.get('vmicro',2.),
                               modelfilename],
                              cwd=modeldirname,
                              stdout=stdout,stderr=stderr)
        stdout.close()
        stderr.close()
    except: raise
    finally:
        os.remove(os.path.join(modeldirname,'makemoogmodel.awk'))
    return None

def vmacro(x,vmacro=6.,sparse=False,norm=True):
    """
    NAME:
       vmacro
    PURPOSE:
       compute the proper macroturbulence kernel
    INPUT:
       x - Array of X values for which to compute the macroturbulence kernel, in pixel offset relative to pixel centers; the kernel is calculated at the x offsets for each pixel center; x need to be 1/integer equally-spaced pixel offsets
       vmacro= (6.) macroturbulence in km/s (FWHM)
       sparse= (False) if True, return a sparse representation that can be passed to apogee.spec.lsf.convolve for easy convolution      
       norm= (True) if False, don't normalize to sum to 1 (useful to check whether the kernel actually integrates to 1)
    OUTPUT:
       LSF-like array of the macroturbulence
    HISTORY:
       2015-03-23 - Written - Bovy (IAS)
    """
    from apogee.spec.lsf import sparsify
    # Convert vmacro to Gaussian sigma / c
    sigvm= vmacro/3./10.**5./2./numpy.sqrt(2.*numpy.log(2.))
    # Are the x unit pixels or a fraction 1/hires thereof?
    hires= int(1./(x[1]-x[0]))
    # Setup output
    wav= apStarWavegrid()
    l10wav= numpy.log10(wav)
    dowav= l10wav[1]-l10wav[0]
    # Hi-res wavelength for output
    hireswav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    # Calculate kernel
    lam= numpy.tile(hireswav,(len(x),1)).T
    dlam= 10.**(numpy.tile(numpy.log10(hireswav),(len(x),1)).T\
                    +numpy.tile(x,(len(hireswav),1))*dowav)/lam-1.
    u= numpy.fabs(dlam/sigvm)
    out= 2./numpy.sqrt(numpy.pi)*u\
        *(numpy.exp(-u**2.)/u-numpy.sqrt(numpy.pi)*special.erfc(u))
    out[dlam == 0.]= 2./numpy.sqrt(numpy.pi)
    out*= (1.+dlam)*numpy.log(10.)/sigvm
    if norm: out/= numpy.tile(numpy.sum(out,axis=1),(len(x),1)).T
    if sparse: out= sparsify(out)
    return out

def _chi2(mspec,spec,specerr,weights=None):
    """Internal function that calculates the chi^2 for a given model,
     assumes that the wavelength axis==-1"""
    if not weights is None:
        return numpy.sum(weights*(mspec-spec)**2./specerr**2,axis=-1)
    else:
        return numpy.sum((mspec-spec)**2./specerr**2,axis=-1)

