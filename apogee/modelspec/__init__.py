from functools import wraps
import os, os.path
import shutil
import subprocess
import numpy
import apogee.tools.read as apread
import apogee.tools.path as appath
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

